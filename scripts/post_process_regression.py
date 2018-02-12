#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Take a completed regression from cisVar.py and create a dataframe, call
open/closed alleles, sort, and tidy.
"""
import sys as _sys
import bz2 as _bz2
import gzip as _gzip
import argparse as _argparse
import logging as _log

import numpy as np
import pandas as pd

# Formatted logging
_log_formatter = _log.Formatter(
    "%(asctime)s [%(name)s] %(levelname)s – %(message)s"
)
_log_formatter_simp = _log.Formatter(
    "%(asctime)s %(levelname)s – %(message)s"
)
log = _log.getLogger('process_regression')
_log_out = _log.StreamHandler()
_log_out.setFormatter(_log_formatter_simp)
_log_out.setLevel(_log.INFO)
log.addHandler(_log_out)
log.setLevel(_log.INFO)

###############
#  Constants  #
###############

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


####################
#  Core Functions  #
####################


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
    log.info('Loading regression file')
    df = pd.read_csv(regfile, sep='\t', low_memory=False, dtype=float_dtypes,
                     float_precision=128)
    assert df.columns.tolist() == orig_headers
    log.info('%s regression records', len(df))
    log.debug('Changing headers')
    df.columns = new_headers
    # Set the index to chr.position
    log.debug('Indexing')
    df['idx'] = df.chrom.astype(str) + '.' + df.position.astype(str)
    df = df.set_index('idx', drop=True)

    # Add rsID if possible
    if bed:
        log.info('Loading rsids')
        df['rsid'] = bed_to_series(bed)
    else:
        log.warn('No bed file, all rsids will be "Unknown"')
        df['rsid'] = 'Unknown'

    # Sort the data by position (default is by best pvalue)
    log.info('Sorting')
    df = df.sort_values(['chrom', 'position'])

    # Force datatypes and uppercase
    df['position'] = df.position.astype(int)
    for i in ['ref', 'alt', 'post_allele']:
        df[i] = df[i].str.upper()
    # Add the beta
    df['beta'] = df.post_freq - df.pre_freq

    # Add open and closed
    log.info('Calculating open/closed')
    df['open_allele'] = np.where(df.post_freq > df.pre_freq, df.post_allele, df.alt)
    df['closed_allele'] = np.where(df.post_freq > df.pre_freq, df.alt, df.post_allele)
    df['open_allele'] = np.where(df.post_freq == df.pre_freq, 'NA', df.open_allele)
    df['closed_allele'] = np.where(df.post_freq == df.pre_freq, 'NA', df.closed_allele)

    # Sort columns
    log.debug('Adjusting columns')
    df = df[final_headers]

    # Write outputs
    if pandasfile:
        log.info('Writing pickled panda to %s', pandasfile)
        df.to_pickle(pandasfile)
    if textfile:
        log.debug('Writing parsed dataframe to %s', textfile)
        df.to_csv(textfile, sep='\t')

    return df


def bed_to_series(bedfile):
    """Convert a bed file to a series of chr.pos->name."""
    args = dict(sep='\t', usecols=[0, 2, 3])
    log.debug('Checking for bed head')
    with _open_zipped(bedfile) as fin:
        if not fin.readline().startswith('#'):
            args['header'] = None
    log.debug('Reading in bed file')
    bed = pd.read_csv(bedfile, **args)
    bed.columns = ['chrom', 'pos', 'name']
    log.debug('Indexing bed')
    bed['idx'] = bed.chrom.astype(str) + '.' + bed.pos.astype(str)
    bed = bed.drop_duplicates('idx', keep=False)
    bed = bed.set_index('idx', drop=True)
    log.debug('Returning rsid series')
    log.debug('Parsed %s SNPs', len(bed))
    return bed.name


#######################
#  Support Functions  #
#######################


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

    inputs = parser.add_argument_group('inputs')
    outputs = parser.add_argument_group('outputs')
    inputs.add_argument('regfile',
                        help="Regression matrix (Default: STDIN)")
    inputs.add_argument('-b', '--bedfile',
                        help="BED file to extract rsIDs from")
    outputs.add_argument('-t', '--textfile',
                         help="Parsed output (Default: STDOUT)")
    outputs.add_argument('-p', '--pandasfile',
                         help="Parsed dataframe")
    outputs.add_argument('-l', '--logfile',
                         help="Log file (Default: STDERR)")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Verbose logging")
    parser.add_argument('-n', '--notext', action='store_true',
                        help="Suppress text output, output only dataframe")

    args = parser.parse_args(argv)

    ifl = args.regfile if args.regfile else _sys.stdin
    ofl = args.textfile if args.textfile else _sys.stdout

    # Adjust logging
    if args.verbose:
        log.setLevel(_log.DEBUG)
    loghandle = None
    if args.logfile:
        # Setup file
        loghandle = _open_zipped(args.logfile, 'a')
        handler = _log.StreamHandler(loghandle)
        handler.setFormatter(_log_formatter)
        if args.verbose:
            handler.setLevel(_log.DEBUG)
        log.addHandler(handler)
        log.info('%s Starting Run', __file__)
    elif args.verbose:
        # Set STDERR logger to debug if no logfile
        log.handlers = []
        _log_out.setLevel(_log.DEBUG)
        log.addHandler(_log_out)

    if args.notext:
        if not args.pandasfile:
            _sys.stderr.write('Pandas file required in notext mode.\n')
            _sys.stderr.write(parser.format_help() + '\n')
            return 1
        parse_regression(ifl, pandasfile=args.pandasfile)
    else:
        parse_regression(ifl, textfile=ofl, pandasfile=args.pandasfile,
                         bed=args.bedfile)

    if loghandle:
        loghandle.close()

if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
