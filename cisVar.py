#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cisVar: Find cis QTLs based on an experimental selection method

Ashley Tehranchi <tehranchi576@gmail.com>

Stanford University

Version: 2.0.0b1
Created: 2015-12-12
Updated: 2018-04-24

Example usage:
cisVar mpileup -F <SampleName> -f <fastaFile> -p <mpileupBEDfile> -B <sortedBam>
cisVar post -F <SampleName> -r <readDepth> -a <allelesFile>
cisVar geno -F <SampleName> -r <readDepth> -i <individualsFile> -g <genotypesFile>
cisVar qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>

Note:
The qtls regression step will use approximately 32GB of memory on an averaged-
sized dataset.

The geno step will use approximately 20GB of memory on the same dataset.

"""
from __future__ import print_function
import os
import sys
import argparse
import operator
import subprocess

import bz2 as _bz2
import gzip as _gzip
from datetime import datetime

# psutil is used to print memory usage if installed, otherwise ignored
try:
    import psutil
except ImportError:
    psutil = None


############################################################################
# mpileup BAM files with hg19 masked genome with blacklist regions removed #
############################################################################


def mpileup(allCHRfasta, mpileupBEDfile, sortedBam, prefix_name):
    """Run mpileup on the provided files."""
    # outputs prefix.mpileup.txt
    subprocess.check_call(
        "samtools mpileup -Q 25 -f " + allCHRfasta + " -l " + mpileupBEDfile \
        + " " + sortedBam + " > " + prefix_name + ".mpileup.txt",
        shell=True
    )


###########################################################################
#                 POST-ChIP Calculation and File Trimming                 #
###########################################################################


def postcalc(prefix_name, trial_depths, allelesFileName):
    """Calculate POST frequencies from mpileup results."""
    print("Post-ChIP calculation for min reads {}".format(trial_depths))

    infilename  = prefix_name + ".mpileup.txt"
    inalleles   = allelesFileName
    outfilename = prefix_name + "." + str(trial_depths) + ".POST.txt"
    depth = int(trial_depths)

    pileup_dictionary = {}
    count = 1
    bad   = 0
    # mpileup format:
    # chromosome, 1-based coordinate, reference base, # reads covering the
    # site, read bases, and base qualities
    with open_zipped(infilename, 'r') as pileup:
        for line in pileup:
            count = count + 1
            rows = line.rstrip().split("\t")
            if not len(rows) == 6:
                bad += 1
                continue

            # SNP name is just chr.position
            snp = rows[0] + "." + rows[1]
            # cov is the total coverage for the SNP
            cov = int(rows[3])
            # Normalize base (called read here) to lowercase to make counting
            # easier
            read = str(rows[4]).lower()

            # Only work on those with coverage exceeding the minimum depth
            if cov < depth:
                Acount = 'NA'
                Ccount = 'NA'
                Gcount = 'NA'
                Tcount = 'NA'
                ALTallele = 'NA'
                ALTdepth = 'NA'
                REFdepth = 'NA'
                RefFreq = 'NA'
                AltFreq = 'NA'
            else:
                # Count all A's, A is plus strand, a in negative strand,
                # count both
                Acount = read.count('a')
                Ccount = read.count('c')
                Gcount = read.count('g')
                Tcount = read.count('t')

                # Make list of all counts
                countslist = [Acount, Ccount, Gcount, Tcount]
                # Find index of largest number and the value
                # ALT is assumed to be most common of the 4 bases
                index, ALTdepth = max(
                    enumerate(countslist), key=operator.itemgetter(1)
                )

                indels = rows[4].count('-') + rows[4].count('+')
                REFdepth = cov - Acount - Ccount - Gcount - Tcount - indels

                # Adjusts total_depth/cov to remove third alleles/sequencing
                # errors
                snp_depth = REFdepth + ALTdepth
                rows[3] = str(snp_depth)

                if index == 0:
                    ALTallele = 'A'
                elif index == 1:
                    ALTallele = 'C'
                elif index == 2:
                    ALTallele = 'G'
                elif index == 3:
                    ALTallele = 'T'
                else:
                    ALTallele = "NA"

                RefFreq = "NA"
                AltFreq = "NA"

                if snp_depth == 0:
                    RefFreq = str(0)
                    print(
                        "Reference frequency is 0 for {}, probably a mistake"
                        .format(snp)
                    )
                else:
                    postfreq = float(REFdepth) / float(snp_depth)
                    if snp_depth > depth and postfreq >= 0:
                        # To get ref freq:
                        # totaldepth - altdepth = refdepth / totaldepth
                        RefFreq = str(postfreq)
                        AltFreq = str(1 - postfreq)

            pileup_dictionary[snp] = (
                rows[2:] + [
                    Acount, Ccount, Gcount, Tcount, ALTdepth,
                    REFdepth, ALTallele, RefFreq, AltFreq
                ]
            )
            assert len(pileup_dictionary[snp]) == 13

    print('{} mpileup rows were too short and thus skipped'.format(bad))

    linecounter = 1
    with open_zipped(inalleles) as alleles, open_zipped(outfilename, 'w') as outfile:
        for line in alleles:
            row = line.rstrip().split('\t')
            snp = row[0] + '.' + row[1]

            if snp in pileup_dictionary:
                matchingVector = [str(x) for x in pileup_dictionary[snp]]
            else:
                continue

            genoRefAllele = row[2]
            genoAltAllele = row[3]

            # For alt alleles with 0 reads, force alt allele to be imputed geno
            # alt allele
            if matchingVector[8] == str(0):
                matchingVector[10] = genoAltAllele

            # For ref alleles with 0 reads, force ref allele to be imputed geno
            # ref allele
            if matchingVector[9] == str(0):
                matchingVector[0] = genoRefAllele

            postChipRefAllele = matchingVector[0].upper()
            postChipRefFreq   = matchingVector[11]
            postChipAltAllele = matchingVector[10].upper()
            postChipAltFreq   = matchingVector[12]

            if postChipRefAllele == genoRefAllele and postChipAltAllele == genoAltAllele:
                outAllele = postChipRefAllele
                outfreq = postChipRefFreq
            elif postChipRefAllele == genoAltAllele and postChipAltAllele == genoRefAllele:
                outAllele = postChipAltAllele
                outfreq = postChipAltFreq
                linecounter = linecounter +1
            else:
                outAllele = 'NA'
                outfreq = 'NA'

            outfile.write(
                '\t'.join(
                    [row[0], row[1], '\t'.join(matchingVector),
                     str(outAllele), str(outfreq)]
                ) + '\n'
            )

    print("post-frequency calculation DONE")


###########################################################################
#      POST.txt file triming, addition of header, BED file creation       #
###########################################################################


def postTrim(prefix_name, trial_depths):
    """Tidy the POST frequency outputs."""

    with open_zipped("mpileup.header.txt", 'w') as headerOUT:
        headerOUT.write(
            '\t'.join(
                ['Chr', 'position', 'REFallele', 'Depth', 'Acount', 'Ccount',
                 'Gcount', 'Tcount', 'ALTdepth', 'REFDepth', 'ALTallele',
                 'REFfreq', 'ALTfreq', 'POSTallele', 'POSTfreq']
            ) + '\n'
        )

    # File names
    prfx = "{}.{}".format(prefix_name, trial_depths)
    post_fl = "{}.POSTth.txt".format(prfx)
    bed_fl = "{}.bed".format(prfx)

    print("File trimming")

    # Filter out unneeded columns
    script = (
        "cat mpileup.header.txt <("
        "sed '/NA/d' {prefix}.POST.txt | cut -f-4,7- "
        ") > {prefix}.POSTth.txt"
    ).format(prefix=prfx)
    subprocess.check_call('bash -c "{}"'.format(script), shell=True)

    print("File trimming DONE")
    print("Writing bed file")

    # Make bed file
    with open_zipped(post_fl) as openfile, open_zipped(bed_fl, 'w') as bedOUT:
        openfile.readline()  # Discard header
        for line in openfile:
            ln = line.rstrip().split('\t')
            snp = ln[1]
            bedOUT.write('\t'.join([ln[0], str(int(snp) - 1), snp]) + '\n')

    print("Bed file writing complete")


###########################################################################
#                 Genotype Extraction for POST-calc SNPs                  #
###########################################################################


def genoExtract(prefix_name, trial_depths, genosFile, individualslist=None):
    """Filter genotypes by SNP and individual.

    First loops through genotypes and filters to a temp file, then transposes,
    then filters transposed lines by individual.

    Takes approximately 10 minutes and uses 10GB of memory on a genotype file
    of 77 million SNPs and an individuals files of 66 individuals.

    Params
    ------
    prefix_name : str
    trial_depths : int
    genosFile : str
        Path to genotype file to parse. Must be in the format:
            chrom\\tpos\\tref\\talt\\tgenotypes...
        Genotypes must be in 0/1/2 format, and the file must contain a header
        with the individual names.
        Note: Genotype file *must* be sorted by position.
    individualslist : str, optional
        Path to file of individuals to include separated by newlines
    """
    new_prefix = prefix_name + "." + str(trial_depths)
    postdata = new_prefix + ".POSTth.txt"
    out_name = new_prefix + ".genotypes.txt"
    print("File prefix used:", new_prefix)

    # Track memory usage
    process = show_mem()

    # Get list of sites with POSTfreq data
    print("Getting SNPs from POST data")
    with open_zipped(postdata, 'r') as postin:
        # Drop first line if it is a header (position is not numberic)
        header = postin.readline().strip().split('\t')
        if header[1].isdigit():
            postin.seek(0)
        postlocs = [
            rows[0] + '.' + rows[1] for rows in [
                line.split('\t') for line in postin.readlines()
            ]
        ]
    post_snps = frozenset(postlocs)
    l = len(postlocs)
    # The POST data must be unique per SNP, if it isn't we have a big problem
    # somewhere in the earlier steps
    assert l == len(post_snps)

    print("Got {} SNPs".format(l))
    show_mem(proc=process)

    done_snps = set()
    add_to_done = done_snps.add
    extras = 0
    indels = 0
    dups   = 0

    print("Filtering genotypes by SNP")
    print("Time: {}".format(datetime.now()))
    # Geno will have the same number of columns as the POST data has rows but
    # will contain numeric matrix only, one column per SNP, and one row per
    # individual *sorted to match individuals.txt*
    with open_zipped(genosFile) as genos_in:
        # Parse header
        header = genos_in.readline().strip().split('\t')
        if header[1].isdigit():
            raise Exception(
                'Genotype file {} has no header, cannot parse individuals'
                .format(genosFile)
            )
        # Extract the individuals from the header
        inds = header[4:]

        # Parse the file
        gt_dict = {}
        for line in genos_in:
            row = line.strip().split('\t')
            # Filter by SNP
            snp = row[0] + '.' + row[1]
            if snp not in post_snps:
                extras += 1
                continue
            if snp in done_snps:
                dups += 1
                continue
            # Filter out INDELS
            if len(row[2]) > 1 or len(row[3]) > 1:
                if len(row[2]) > 1:
                    raise Exception(
                        'Genotype file has more than one char for ref at '
                        'line:\n{}'.format(line)
                    )
                for allele in row[3].split(','):
                    if len(allele) > 1:
                        indels += 1
                        continue
            # Keep only the genotype data, not the position data
            gt_dict[snp] = row[4:]
            # Done by adding this at the end, we can still get duplicates
            # where the first entry was an indel or something else.
            add_to_done(snp)

    print("Genotype filtering complete")
    show_mem(proc=process)

    print("Time: {}".format(datetime.now()))
    print("Dropped {} unneeded SNPs".format(extras))
    print("Dropped {} indels".format(indels))
    print("Dropped {} duplicates".format(dups))

    print("Checking all POST SNPs matched")
    if len(postlocs) != len(done_snps):
        err_file = new_prefix + ".missing.snps"
        with open_zipped(err_file, "w") as fout:
            fout.write(
                '\n'.join(sorted([i for i in postlocs if i not in done_snps]))
            )
        raise Exception(
            "{} SNPs in POST were not in the genotype file, written to {}"
            .format(len(postlocs)-len(done_snps), err_file)
        )
    print("Done")

    # Get individual locations
    if individualslist:
        print(
            "Getting individuals from individuals file ({})."
            .format(individualslist)
        )
        with open_zipped(individualslist) as fin:
            individuals = [i.strip() for i in fin.read().strip().split('\n')]
        ind_locs = []
        missing = []
        for ind in individuals:
            if ind in inds:
                ind_locs.append(inds.index(ind))
            else:
                missing.append(ind)
        if len(missing) > 0:
            sys.stderr.write(
                'The following individuals are missing '
                'from the genotype file:\n'
            )
            sys.stderr.write('\n'.join(missing))
            raise IndexError('Missing some individuals from genotypes')
        assert len(ind_locs) == len(individuals)
    else:
        print('No individuals file provided, keeping all individuals')
        ind_locs = None

    print("Creating sorted matrix of only required individuals")
    mat = []
    for snp in postlocs:
        if ind_locs:
            gt = gt_dict[snp]
            mat.append([gt[i] for i in ind_locs])
        else:
            mat.append(gt_dict[snp])
    print("Done")
    show_mem(proc=process)
    del gt_dict
    print("Number of individuals included: {}".format(len(mat[0])))

    print("Transposing and writing")
    with open_zipped(out_name, 'w') as fout:
        fout.write(
            '\n'.join(
                ['\t'.join(i) for i in zip(*mat)]
            )
        )
    print("Done")
    show_mem(proc=process)
    del mat

    print("Genotype filtering complete")
    show_mem(proc=process)


###########################################################################
#                  Regressuion with Z-score and P-value                   #
###########################################################################


def regPVAL(prefix_name, trial_depths, numIndv):
    """Run the actual regression R code."""
    new_prefix = prefix_name + "." + str(trial_depths)
    postdata = new_prefix + ".POSTth.txt"
    genosout_name = new_prefix + ".genotypes.txt"

    ### make sure location of R script is correct
    r_script = os.path.join(os.path.dirname(__file__), 'regression_qtls.R')
    if not os.path.isfile(r_script):
        raise OSError('File not found: {}'.format(r_script))

    subprocess.check_call(
        ['Rscript', r_script, postdata, genosout_name, new_prefix, numIndv]
    )

    print("regression DONE")


###########################################################################
#                              Tidy Results                               #
###########################################################################

orig_headers = [
    'Chr', 'position', 'REFallele', 'Depth', 'Acount', 'Ccount', 'Gcount',
    'Tcount', 'ALTdepth', 'REFDepth', 'ALTallele', 'REFfreq', 'ALTfreq',
    'POSTallele', 'POSTfreq', 'prechipfreq', 'pvalue', 'zvalue',
    'prevar', 'postvar', 'SNPpostfreq', 'SNPprefreq'
]
new_headers  = [
    'chrom', 'position', 'ref', 'depth', 'a_count', 'c_count', 'g_count',
    't_count', 'alt_depth', 'ref_depth', 'alt', 'ref_freq', 'alt_freq',
    'post_allele', 'post_freq', 'pre_freq', 'p_value', 'z_value',
    'pre_variance', 'post_variance', 'snp_postfreq', 'snp_prefreq'
]
final_headers = [
    'chrom', 'position', 'rsid', 'open_allele', 'closed_allele',
    'pre_freq', 'post_freq', 'pre_variance', 'post_variance', 'beta',
    'p_value', 'z_value', 'ref', 'alt', 'depth', 'ref_depth',
    'alt_depth', 'ref_freq', 'alt_freq', 'snp_postfreq', 'snp_prefreq',
    'a_count', 'c_count', 'g_count', 't_count'
]
float_dtypes = {
    'position': 'object',
    'REFfreq': 'float128',
    'ALTfreq': 'float128',
    'POSTfreq': 'float128',
    'prechipfreq': 'float128',
    'pvalue': 'float128',
    'zvalue': 'float128',
    'prevar': 'float128',
    'postvar': 'float128',
    'SNPpostfreq': 'float128',
    'SNPprefreq': 'float128',
}


def parse_regression(regfile, textfile=None, pandasfile=None, bed=None):
    """Parse a regression file into a pandas datafram.

    Optionally output either pickled panda or textfile.

    For all files, gzipped files or open file handle are permissable.

    Parameters
    ----------
    regfile : file name
        Tab delimited output from regression, header strictly enforced
    textfile : file name, optional
        Path to write tab delimited formatted output
    pandasfile : file name, optional
        Path to write pickled pandas dataframe
    bed : file name, optional
        Path to a bed file that contains rsids

    Returns
    -------
    pandas.DataFrame
        Parsed DataFrame
    """
    import numpy as np
    import pandas as pd
    print('Loading regression file')
    df = pd.read_csv(regfile, sep='\t', low_memory=False, dtype=float_dtypes,
                     float_precision=128)
    assert df.columns.tolist() == orig_headers
    df.columns = new_headers
    # Set the index to chr.position
    df['idx'] = df.chrom.astype(str) + '.' + df.position.astype(str)
    df = df.set_index('idx', drop=True)

    # Add rsID if possible
    if bed:
        print('Loading rsids')
        df['rsid'] = bed_to_series(bed)
    else:
        print('No bed file, all rsids will be "Unknown"')
        df['rsid'] = 'Unknown'

    # Sort the data by position (default is by best pvalue)
    print('Sorting')
    df = df.sort_values(['chrom', 'position'])

    # Force datatypes and uppercase
    df['position'] = df.position.astype(int)
    for i in ['ref', 'alt', 'post_allele']:
        df[i] = df[i].str.upper()
    # Add the beta
    df['beta'] = df.post_freq - df.pre_freq

    # Add open and closed
    print('Calculating open/closed')
    df['open_allele'] = np.where(df.post_freq > df.pre_freq, df.post_allele, df.alt)
    df['closed_allele'] = np.where(df.post_freq > df.pre_freq, df.alt, df.post_allele)
    df['open_allele'] = np.where(df.post_freq == df.pre_freq, 'NA', df.open_allele)
    df['closed_allele'] = np.where(df.post_freq == df.pre_freq, 'NA', df.closed_allele)

    # Sort columns
    print('Adjusting columns')
    df = df[final_headers]

    # Write outputs
    if pandasfile:
        print('Writing pickled panda to', pandasfile)
        df.to_pickle(pandasfile)
    if textfile:
        print('Writing parsed dataframe to', textfile)
        df.to_csv(textfile, sep='\t')

    return df


###########################################################################
#                                 Helpers                                 #
###########################################################################


def bed_to_series(bedfile):
    """Convert a bed file to a series of chr.pos->name."""
    import pandas as pd
    args = dict(sep='\t', usecols=[0, 2, 3])
    with open_zipped(bedfile) as fin:
        if not fin.readline().startswith('#'):
            args['header'] = None
    print('Reading in bed file')
    bed = pd.read_csv(bedfile, **args)
    bed.columns = ['chrom', 'pos', 'name']
    bed['idx'] = bed.chrom.astype(str) + '.' + bed.pos.astype(str)
    bed = bed.drop_duplicates('idx', keep=False)
    bed = bed.set_index('idx', drop=True)
    print('Parsed %s SNPs', len(bed))
    return bed.name


def show_mem(proc=None):
    """Print memory usage."""
    if not psutil:
        return None
    if not proc:
        proc = psutil.Process(os.getpid())
    print(
        "Memory: {:.2f}GB"
        .format(proc.memory_info().rss/1024/1024/1024)
    )
    return proc


def open_zipped(infile, mode='r'):
    """Return file handle of file regardless of compressed or not.

    Returns already opened files unchanged
    Text mode automatic for compatibility with python2.
    """
    # Return already open files
    if hasattr(infile, 'write'):
        return infile

    # Make text mode automatic
    if len(mode) == 1:
        mode = mode + 't'

    # Refuse to handle non-strings that aren't files.
    if not isinstance(infile, str):
        raise ValueError("i cannot open a filename that isn't a string.")

    # Treat '-' appropriately
    if infile == '-':
        if 'w' in mode:
            return sys.stdout
        return sys.stdin

    # If possible open zipped files
    if infile.endswith('.gz'):
        return _gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(_bz2, 'open'):
            return _bz2.open(infile, mode)
        return _bz2.bz2file(infile, mode)

    # Fall back on regular open
    return open(infile, mode)


###########################################################################
#                          Command Line Parsing                           #
###########################################################################


def main():
    """Run as a script."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Subcommands
    modes = parser.add_subparsers(
        dest='inputCommand',
        metavar='{mpileup,post,geno,qtls}'
    )

    # Shared commands
    shared_cmds = argparse.ArgumentParser(add_help=False)
    run_opts = shared_cmds.add_argument_group("Run Options")
    run_opts.add_argument(
        '-F', '--SampleName', dest='prefix_name', default='cis_var',
        help='sample/population name'
    )
    run_opts.add_argument(
        '-r', '--readDepth', type=int, dest='trial_depths', default=20,
        help='minimum read depth per variant'
    )


    # mpileup
    mpileup_cmd_grp = modes.add_parser(
        'mpileup', description="Run mpileup",
        parents=[shared_cmds], help="Run mpileup",
        aliases=['m'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mpileup_cmd = mpileup_cmd_grp.add_argument_group(
        "mpileup Options (Required)"
    )
    mpileup_cmd.add_argument(
        '-f', '--fasta', required=True, dest='allCHRfasta',
        help='fasta file with all chromosomes'
    )
    mpileup_cmd.add_argument(
        '-B', '--BAMfile', required=True, dest='sortedBam',
        help='sorted BAM file'
    )
    mpileup_cmd.add_argument(
        '-p', '--mpileupBEDfile', required=True, dest='mpileupBEDfile',
        help='The  mpileup BED file'
    )

    # POST
    post_cmd_grp = modes.add_parser(
        'post', description="Run POST frequency calculation",
        parents=[shared_cmds], help="Run POST frequency calculation",
        aliases=['p'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    post_cmd = post_cmd_grp.add_argument_group("POST Options (Required)")
    post_cmd.add_argument(
        '-a', '--allelesFile', required=True, dest='allelesFileName',
        help='format: chr1    10583   G   A'
    )

    # geno
    geno_cmd_grp = modes.add_parser(
        'geno', description="Munge genotypes to prepare for regression",
        parents=[shared_cmds],
        help="Munge genotypes to prepare for regression", aliases=['g'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    geno_cmd = geno_cmd_grp.add_argument_group("Genotype Options (Required)")
    geno_cmd.add_argument(
        '-g', '--genoFile', required=True, dest='genosFile',
        help='The genotypes file'
    )
    geno_cmd.add_argument(
        '-i', '--individualsFile', required=False, dest='individualslist',
        help='list of individuals matching genotype matrix; one indv per line'
    )

    # Regression
    qtls_cmd_grp = modes.add_parser(
        'qtls', description="Run the regression",
        parents=[shared_cmds], aliases=['q', 'regression', 'r'],
        help="Run the regression",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    qtls_cmd = qtls_cmd_grp.add_argument_group("Regression Options (Required)")
    qtls_cmd.add_argument(
        '-n', '--numberIndividuals', required=True, dest='numIndv',
        help='The number of individuals in the pool'
    )

    # Tidy
    tidy_cmd = modes.add_parser(
        'tidy', description="Tidy up regression, call open/closed",
        parents=[shared_cmds],
        help="Tidy up regression, call open/clsoed", aliases=['t'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    inputs  = tidy_cmd.add_argument_group('inputs')
    outputs = tidy_cmd.add_argument_group('outputs')
    inputs.add_argument('-b', '--bedfile',
                        help="BED file to extract rsIDs from (optional)")
    outputs.add_argument('-t', '--textfile',
                         help="Parsed output (Default: STDOUT)")
    outputs.add_argument('-p', '--pandasfile',
                         help="Parsed dataframe")


    args = parser.parse_args()

    if args.inputCommand in ['mpileup', 'm']:
        mpileup(
            args.allCHRfasta, args.mpileupBEDfile,
            args.sortedBam, args.prefix_name
        )

    elif args.inputCommand in ['post', 'p']:
        postcalc(args.prefix_name, args.trial_depths, args.allelesFileName)
        postTrim(args.prefix_name, args.trial_depths)

    elif args.inputCommand in ['geno', 'g']:
        genoExtract(
            args.prefix_name, args.trial_depths,
            args.genosFile, args.individualslist
        )

    elif args.inputCommand in ['qtls', 'q', 'regression', 'r']:
        regPVAL(args.prefix_name, args.trial_depths, args.numIndv)

    elif args.inputCommand in ['tidy', 't']:
        prefix = args.prefix_name + '.' + str(args.trial_depths)
        ifl = prefix + '.total.txt'
        ofl = args.textfile if args.textfile else sys.stdout
        parse_regression(
            ifl, textfile=ofl, pandasfile=args.pandasfile,
            bed=args.bedfile
        )

    else:
        sys.stderr.write('Unrecognized command {}'.format(args.inputCommand))
        parser.print_help(file=sys.stderr)
        return 2


if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
