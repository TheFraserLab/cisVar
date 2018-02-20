#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert a bed file into an alleles file.
"""
import os as _os
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse
import multiprocessing as mp

try:
    import fyrd
except ImportError:
    fyrd = False


def make_alleles(mpileup_file, vcf_files, allele_file, use_fyrd=True, cores=1):
    """Use VCF files plus the mpileup file to make alleles file.

    Note: We only use the mpileup file to limit the data in the vcf files.
    Arguments should be writable files or file names
    """
    if isinstance(vcf_files, str):
        vcf_files = [vcf_files]
    else:
        vcf_files = list(vcf_files)
    jobs = {}
    if use_fyrd and fyrd:
        print('Using fyrd, cores ignored')
    else:
        print('Using multiprocessing')
        pool = mp.Pool(int(cores))
    for vcf in vcf_files:
        if use_fyrd and fyrd:
            print('Using fyrd, cores ignored')
            jobs[vcf] = fyrd.submit(
                get_vcf, (vcf,),
                mem='12GB', time='04:00:00',
                syspaths=[_os.path.dirname(_os.path.abspath(__file__))],
                imports=[
                    'from make_alleles import *',
                    'from make_alleles import open_zipped'
                ]
            )
        else:
            jobs[vcf] = pool.apply_async(get_vcf, (vcf,))

    all_loci = {}
    if use_fyrd and fyrd:
        fyrd.wait(jobs.values())

    all_loci = {}
    for vcf, job in jobs.items():
        print('Getting', vcf)
        out = job.get()
        all_loci.update(out)

    other = {}
    for k, v in all_loci.items():
        other[flip_chrom(k)] = v
    all_loci.update(other)

    with open_zipped(mpileup_file) as fin, open_zipped(allele_file, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            ref, alt = all_loci[fields[0]][fields[1]]
            if len(ref) > 1 or len(alt) > 1:
                continue
            fout.write(
                '{}\t{}\t{}\t{}\n'.format(
                    fields[0],
                    fields[1],
                    ref,
                    alt
                )
            )
    return 0


def get_vcf(vcf):
    """Convert a vcf file to a dict of 'chr'->'pos'->(ref, alt)
    """
    all_loci = {}
    with open_zipped(vcf) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[0] not in all_loci:
                all_loci[fields[0]] = {}
            all_loci[fields[0]][fields[1]] = (fields[3], fields[4])
    return all_loci


def flip_chrom(x):
    """Return alternate chromosome formatting."""
    if x.startswith('chr'):
        return x[3:]
    else:
        return 'chr{}'.format(x)


def open_zipped(infile, mode='r'):
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

    parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    # Optional arguments
    parser.add_argument('-c', '--cores', default=1, type=int,
                        help='Number of cores to use')
    parser.add_argument('--skip-fyrd', action='store_false',
                        help='Do not use fyrd even if available')

    parser.add_argument('mpileup',
                        help='The M-pileup file to match')

    # Optional files
    parser.add_argument('allele_file', nargs='?',
                        help="Output file (Default: STDOUT)")
    parser.add_argument('vcf_files', nargs='+',
                        help="VCF Files")

    args = parser.parse_args(argv)

    ofl = args.allele_file if args.allele_file else _sys.stdout

    return make_alleles(
        args.mpileup, args.vcf_files, ofl,
        use_fyrd=args.skip_fyrd, cores=args.cores
    )

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
