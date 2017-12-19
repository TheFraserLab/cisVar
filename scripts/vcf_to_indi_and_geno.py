#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parse VCF files into genotypes and individuals files.

EX:
    chr   position    ref alt NA18486 NA18489 ...
    chr1  10583   g   a   0   1 ...
    chr1  10611   C   G   0   0 ...
    chr1  13302   c   t   1   0 …

and:
    NA18486
    NA18489
"""
import os as _os
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import subprocess as _sub
import argparse as _argparse

import fyrd


EXPECTED_HEADER = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
]


def geno_file(vcf, check_header=None, ind_list=None,
              get_header=False, get_inds=False, log=_sys.stderr):
    """Wrapper for genotype files, currently only vcf

    Parameters
    ----------
    vcf : str
        Path to a vcf file, can be gzipped
    check_header : list, optional
        If a list is provided, assert that header matches before iterating.
    ind_list : list of str, optional
        If list provided, filter output to only include these individuals.
    get_header : bool, optional
        Return a header instead of yielding a row
    get_inds : bool, optional
        Return a list of individuals instead of yielding a row

    Yields
    ------
    chr : str
    pos : str (base 1)
    ref : str
    alt : str
    inds : list of str
    """
    name = vcf.split('.')
    if 'vcf' in name[-2:]:
        return vcf_file(vcf, check_header, ind_list, get_header, get_inds, log)
    else:
        raise NotImplementedError('Currently only VCFs supported')


def vcf_file(vcf, check_header=None, ind_list=None,
             get_header=False, get_inds=False, log=_sys.stderr):
    """VCF file iterator with goodies

    Parameters
    ----------
    vcf : str
        Path to a vcf file, can be gzipped
    check_header : list, optional
        If a list is provided, assert that header matches before iterating.
    ind_list : list of str, optional
        If list provided, filter output to only include these individuals.
    get_header : bool, optional
        Return a header instead of yielding a row
    get_inds : bool, optional
        Return a list of individuals instead of yielding a row
    bad : file
        File or handle to write errors

    Yields
    ------
    chr : str
    pos : str (base 1)
    ref : str
    alt : str
    inds : list of str
    """
    # Get header at genotype loc
    with _open_zipped(vcf) as fin:
        line = fin.readline()
        while True:
            if not line.startswith('#'):
                break
            pos    = fin.tell()
            header = line
            line   = fin.readline()

    # Get genotype location
    try:
        gti = line.split('\t')[8].split(':').index('GT')
    except ValueError:
        raise ValueError(
            'GT column missing from VCF format string. In need '
            'the genotype info location to know where to look '
            'for the genotype (looks like 0|1)'
        )

    headers = header.strip().split('\t')
    inds = headers[9:]
    if check_header and header != check_header:
        sys.stderr.write('{} bad header\n'.format(vcf))
        sys.stderr.write('Should be:\n{}\nInd:\n{}\n'.format(check_header, header))
        raise Exception('{} file had an invalid header'.format(vcf))
    if get_header:
        return header
    if get_inds:
        return inds

    return _vcf_file(vcf, inds, ind_list, gti, pos, log)


def _vcf_file(vcf, inds, ind_list, gti, pos, log):
    """Iterator for VCF."""

    # Get individual locations for lookup/filter
    ind_pos = []
    if ind_list:
        if sorted(ind_list) != sorted(inds):
            for ind in ind_list:
                ind_pos.append(inds.index(ind))

    short_lines = 0
    bad_gts = 0
    count = 0

    # Actually parse file
    with _open_zipped(vcf) as fin:
        fin.seek(pos)
        count = 0
        for line in fin:
            f = line.strip().split('\t')
            out = [f[0], f[1], f[3], f[4]]
            gti = f[8].split(':').index('GT')
            ind_data = f[9:]
            if not len(ind_data) == len(inds):
                short_lines += 1
                continue
            if ind_pos:
                these_inds = []
                for loc in ind_pos:
                    these_inds.append(ind_data[loc])
                assert len(these_inds) == len(ind_list)
            else:
                these_inds = ind_data
            for ind in these_inds:
                gt = parse_gt(ind.split(':')[gti])
                if not gt:
                    #  _sys.stderr.write(
                        #  'Bad line: {}\t{}\n'.format(count, line)
                    #  )
                    bad_gts += 1
                    continue
                out.append(gt)
            yield out
    log.write('Total: {}\nBad Genotypes: {}\nShort Lines (not enough inds): {}\n'
              .format(count, bad_gts, short_lines))
    return


def parse_gt(gt):
    """Convert common genotypes into 0,1,2."""
    if '|' in gt:
        if gt == '0|0':
            gt = '0'
        elif gt in ['0|1', '1|0']:
            gt = '1'
        elif gt == '1|1':
            gt = '2'
        else:
            gt = None
    else:
        gt = None
    return gt


def parse_files(geno_files, geno_output, individuals, skip_indels=True,
                log_dir='.'):
    """Loop through all geno files in parallel and then recombine.

    Parameters
    ----------
    geno_files : list of str
        Genotype file locations, currently only VCF supported, gzipped fine.
    geno_output : str
        file to write genotypes
    individuals : str
        Path to a file of individuals, one per row. If file does not exist,
        will be written and all will be kept.
    skip_indels : bool, optional
        Don't include rows with indels
    log_dir : str, optional
        Directory to write parse statistics to
    """
    if not isinstance(geno_files, list):
        raise ValueError('genotype files must be a list')

    log_dir = _os.path.abspath(_os.path.expandvars(log_dir))
    if not _os.path.isdir(log_dir):
        _os.makedirs(log_dir)

    # Individuals
    if _os.path.isfile(individuals):
        print('Using existing individuals file.')
        with open(individuals) as fin:
            inds = fin.read().strip().split('\n')
        final_inds = inds
    else:
        print('Creating new individuals file with all individuals.')
        inds = geno_file(
            geno_files[0], get_inds=True
        )
        with open(individuals, 'w') as fout:
            fout.write('\n'.join(inds))
        final_inds = inds
        print('Individuals written')
        inds = None

    print('\nUsing {} individuals.\n'.format(len(final_inds)))

    # Header
    header = geno_file(
        geno_files[0], get_header=True
    )

    # Parse files
    count = 0
    print('Looping through VCFs')
    jobs = []
    ourdir, us = _os.path.split(__file__)
    us = us.split('.')[0]
    for geno in geno_files:
        log_fl = _os.path.join(log_dir, _os.path.basename(geno) + '.log')
        count += 1
        ofl = '{}.{}'.format(geno_output, count)
        if geno_output.endswith('gz'):
            ofl = ofl + '.gz'
        jobs.append(
            (
                ofl,
                fyrd.submit(
                    parse_file, (geno, ofl, header, inds, skip_indels, log_fl),
                    time="04:00:00", mem='12GB', syspaths=[ourdir],
                    scriptpath=log_dir, outpath=log_dir, clean_files=False,
                    clean_outputs=False, partition='high_priority',
                    imports=[
                        'from {} import *'.format(us),
                        'from {} import _open_zipped'.format(us),
                        'from {} import _vcf_file'.format(us)
                    ]
                )
            )
        )

    j = [i[1] for i in jobs]
    fyrd.wait(j)

    print('Checking jobs')

    # Check jobs
    failed = []
    for ofl, job in jobs:
        if not job.state == 'completed':
            failed.append((ofl, job))
    if failed:
        print('Failed: {}'.format(failed))
        raise Exception('Failed')

    fls = [i[0] for i in jobs]
    for fl in fls:
        if not _os.path.isfile(fl) or not _os.path.isfile(fl + '.done'):
            failed.append(fl)
    if failed:
        raise Exception('The following files failed: {}'.format(failed))

    # Merge
    print('Merging files')
    with open(geno_output, 'w') as fout:
        fout.write(
            'chrom\tposition\tref\talt\t{}\n'.format('\t'.join(final_inds))
        )
    _sub.check_call(
        'cat {} >> {}'.format(
            ' '.join(fls), geno_output
        ),
        shell=True
    )

    print('Removing tmp files')
    for temp_file in fls:
        if _os.path.isfile(temp_file):
            _os.remove(temp_file)


def parse_file(geno, outfile, header, inds=None, skip_indels=False, log=None):
    """File parser, see parse_files for information."""
    if not log:
        log = geno + '.parse.log'
    indels = 0
    kept = 0
    superbad = []
    print('Creating', outfile)
    with _open_zipped(outfile, 'w') as gn, open(log, 'a') as lg:
        lg.write('Parse initiated\n')
        for f in geno_file(geno, check_header=header, ind_list=inds, log=lg):
            # Write core records
            if skip_indels and (len(f[2]) > 1 or len(f[3]) > 1):
                indels += 1
                continue
            chrom = f[0] if f[0].startswith('chr') else 'chr' + f[0]
            outstr = (
                '{chr}\t{pos}\t{ref}\t{alt}'.format(
                    chr=chrom, pos=f[1], ref=f[2], alt=f[3]
                )
            )
            # Write genotypes
            if inds:
                if len(f[4:]) != len(inds):
                    superbad.append(f)
                    continue
            outstr += '\t{}'.format('\t'.join(f[4:]))
            gn.write(outstr + '\n')
            kept += 1
        lg.write('Indels: {}\nKept: {}\nParse Complete\n'.format(indels, kept))
        if superbad:
            lg.write('Length failure despite length filtering:\n{}\n'
                     .format(superbad))
    with open(outfile + '.done') as fout:
        fout.write('Done\n')
    return


def _open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Also returns already opened files unchanged, text mode automatic for
    compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile
    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'
    if not isinstance(infile, str):
        raise ValueError("I cannot open a filename that isn't a string.")
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        else:
            return _bz2.BZ2File(infile, mode)
    return open(infile, mode)


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('genotype_file', help="Genotype output file")
    parser.add_argument('individual_file',
                        help=("Individual file, create if does not exist, "
                              "if exists, filter by these inividuals"))
    parser.add_argument('vcf_files', nargs='+', help="VCF Files")

    parser.add_argument('-s', '--skip_indels', action='store_true',
                        help="Don't include indels in output")
    parser.add_argument('-l', '--log-dir', default='.',
                        help="Log directory")

    args = parser.parse_args(argv)

    parse_files(
        args.vcf_files, args.genotype_file,
        args.individual_file, skip_indels=args.skip_indels,
        log_dir=args.log_dir
    )


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
