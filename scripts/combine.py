#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Combine a bunch of cisVar pandas files by sample (e.g. population).

Requires a search string such as prefix.{sample}.regression.pd.

Writes
------
prefix.combined.regression.pd
    A simple merger of all DataFrames
prefix.merged.regression.pd
    A per-snp merger based on p-value
"""
import os as _os
import re as _re
import sys as _sys
import argparse as _argparse
from glob import glob as _glob

import numpy as np
import pandas as pd
from scipy import stats as sts


def find_dfs(search_str):
    """Use search_str to find files to merge.

    Params
    ------
    search_str : str
        String in format 'prefix.{sample}.suffix', {sample} required

    Returns
    -------
    dfs : dict
        Dictionary of {name: path} for all dataframes
    """
    query = '{sample}'
    if query not in search_str:
        raise ValueError('Search string must contain {sample}')
    ss = search_str.replace(query, '*')
    q  = _re.compile(search_str.replace(query, '(.*)'))
    c  = _re.compile(search_str.replace(query, 'combined'))
    return {
        q.findall(i)[0]: i for i in _glob(_os.path.expandvars(ss))
        if not c.search(i)
    }


def combine_dfs(dfs, column_name='population'):
    """Combine pandas dfs.

    Params
    ------
    dfs : dict
        {name: path}
    column_name : str, optional
        A name for the column to merge on

    Returns
    -------
    DataFrame
    """
    cdfs = []
    for name, path in dfs.items():
        print('Working on {}'.format(path))
        df = pd.read_pickle(path)
        df[column_name] = name
        df['snp'] = df.index.to_series()
        df['idx'] = df[column_name].astype(str) + '.' + df['snp'].astype(str)
        cdfs.append(df)
    print('Combining')
    return pd.concat(cdfs)


def combine_and_write(search_str, column_name='population'):
    """Combine pandas dfs by search string.

    Params
    ------
    search_str : str
        String in format 'prefix.{sample}.suffix', {sample} required
    column_name : str, optional
        A name for the column to merge on

    Writes
    ------
    DataFrame
        Combined DataFrame at search_string.replace('{sample}', 'combined')

    Returns
    -------
    DataFrame
    """
    comb = combine_dfs(find_dfs(search_str), column_name=column_name)
    ofl = search_str.replace('{sample}', 'combined')
    print('Writing to {}'.format(ofl))
    comb.to_pickle(ofl)
    return comb


def merge_dataframe(df, column_name='population'):
    """Merge the result of combine_dfs into per-SNP dataframe."""
    df['snp'] = df.chrom.astype(str) + '.' + df.position.astype(str)
    print('Getting open and closed combined data, this might take a while.')
    # Takes about 20 minutes
    open_closed = get_open_closed_counts(df)
    keep_cols = [
        'snp', 'chrom', 'position', 'rsid', 'pre_freq', 'post_freq',
        'pre_variance', 'post_variance', 'beta', 'p_value', 'z_value', 'ref',
        'alt', 'depth', 'ref_depth', 'alt_depth', 'snp_postfreq', 'snp_prefreq',
        'population'
    ]
    comb_functions = {
        'chrom': first, 'position': first, 'rsid': first,
        'ref': first, 'alt': first,
        'pre_freq': np.median, 'post_freq': np.median,
        'pre_variance': np.median, 'post_variance': np.median,
        'beta': np.median, 'z_value': np.median, 'depth': np.median,
        'ref_depth': np.median, 'alt_depth': np.median,
        'snp_postfreq': np.median, 'snp_prefreq': np.median,
        column_name: ','.join, 'p_value': fishers
    }
    print('Collapsing dataframe, this also takes a while.')
    # Takes about 10 minutes
    df = df[keep_cols].groupby('snp').agg(comb_functions)
    # Split the fishers statistics
    df['fishers_value'], df['fishers_p'] = df.p_value.str
    print('Merging in open/closed data')
    df = pd.merge(df, open_closed, left_index=True, right_index=True)
    print('Organizing')
    final_cols = [
        'chrom', 'position', 'rsid', 'fishers_p', 'ref', 'alt', 'open_best',
        'closed_best', 'pop_count', column_name, 'pops_agree',
        'most_pops_agree', 'all_open_closed', 'open_closed_freqs',
        'fishers_value', 'pre_freq', 'post_freq', 'pre_variance',
        'post_variance', 'beta', 'z_value', 'depth', 'ref_depth', 'alt_depth',
        'snp_postfreq', 'snp_prefreq'
    ]
    final_names = [
        'chrom', 'position', 'rsid', 'fishers_p', 'ref', 'alt', 'open_best',
        'closed_best', column_name + '_count', column_name + 's', column_name +
        's_agree', 'most_{}s_agree'.format(column_name), 'all_open_closed',
        'open_closed_freqs', 'fishers_value', 'pre_freq_median',
        'post_freq_median', 'pre_variance_median', 'post_variance_median',
        'beta_median', 'z_value_median', 'depth_median', 'ref_depth_median',
        'alt_depth_median', 'snp_postfreq_median', 'snp_prefreq_median'
    ]
    df = df[final_cols]
    df.columns = final_names
    print('Done')
    return df


def merge_and_write(df, ofl, column_name='population'):
    """Merge the dataframe and write."""
    df = merge_dataframe(df, column_name=column_name)
    df.to_pickle(ofl)
    return df


###########################################################################
#                          Combination Functions                          #
###########################################################################

def first(x):
    """Return the first member of a list."""
    try:
        x = x.to_series()
    except AttributeError:
        pass
    return list(x)[0]


def fishers(x):
    """Return fishers, fishers_p for pvalues."""
    return sts.combine_pvalues(x.astype(np.float), method='fisher')



def summarize(x):
    """Create a set of open, closed allele pairs.

    Args:
        DataFrame: Should be provided by apply, expects only one SNP

    Returns:
        tuple: tuple, int, tuple: open_closed_alleles, total populations, open_close_frequencies
    """
    pop_count = len(x.open_allele)
    # Get the counts for every open, closed combination
    oc = pd.Series([(i, j) for i, j in zip(x.open_allele.tolist(), x.closed_allele.tolist())]).value_counts()
    total = oc.sum()
    ocs = []
    freqs = []
    for oc, cnt in zip(oc.index.to_series().tolist(), oc.tolist()):
        ocs.append(oc)
        freqs.append(round(cnt/total, 3))
    return tuple(ocs), total, tuple(freqs)


def pops_agree(x):
    """True only if there is only one open-closed pair in the row, use with apply, expects a row of a DF."""
    return len(x.all_open_closed) == 1


def pops_agree_50(x):
    """True if more than 50% of populations agree on open-closed pair, False if 50% or less agree."""
    return x.open_closed_freqs[0] > .5


def get_open_closed_counts(x):
    """Take a dataframe (from apply), tally the open/closed SNPs, and return some stats.

    Returns:
        DataFrame: open and closed _best represent the most common open/closed alleles, pop count
                   is all populations with this SNP, all_open_closed is a tuple of tuples of all open-closed pairs,
                   open_closed_freqs tells you how many populations had each of the open-closed pairs. Other
                   two columns do what they say on the tin.
    """
    open_closed_data = pd.concat(x[['snp', 'open_allele', 'closed_allele']].groupby('snp').apply(summarize).str, axis=1)
    open_closed_data.columns = ['all_open_closed', 'pop_count', 'open_closed_freqs']
    open_closed_data.index.name = None
    # Explicitly set the best alleles, assuming sorted by frequency (which is done by summarize())
    open_closed_data['open_best'], open_closed_data['closed_best'] = open_closed_data['all_open_closed'].str[0].str
    open_closed_data['pops_agree'] = open_closed_data.apply(pops_agree, axis=1)
    open_closed_data['most_pops_agree'] = open_closed_data.apply(pops_agree_50, axis=1)
    # Set best alleles to NA if the top allele is not better than next best allele
    top_freq = open_closed_data.open_closed_freqs.str[0]
    nxt_freq = np.where(open_closed_data.pops_agree, 0, open_closed_data.open_closed_freqs.str[1])
    for i in ['open_best', 'closed_best']:
        open_closed_data[i] = np.where(top_freq > nxt_freq, open_closed_data[i], 'NA')
    # Return organized DF
    return open_closed_data[['open_best', 'closed_best', 'pop_count', 'all_open_closed', 'pops_agree', 'most_pops_agree', 'open_closed_freqs']]


def main(argv=None):
    """Run as a script."""
    if not argv:
        argv = _sys.argv[1:]

    parser  = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter
    )

    # Positional arguments
    parser.add_argument(
        'search_str',
        help='e.g. name.{sample}.regression.pd, used to find files'
    )
    parser.add_argument(
        '-c', '--column-name', default='population',
        help='Name for combine column, e.g. population'
    )
    parser.add_argument(
        '--no-merge', action='store_false', dest='merge',
        help='Produce only a combined dataframe, not a merged dataframe' +
        '. merging can add half an hour over combination, which takes seconds'
    )

    args = parser.parse_args(argv)

    comb = combine_and_write(args.search_str, column_name=args.column_name)
    if len(comb) < 10:
        return 1
    if args.merge:
        ofl = args.search_str.replace('{sample}', 'merged')
        merged = merge_and_write(comb, ofl, column_name=args.column_name)
    return 0


if __name__ == '__main__' and '__file__' in globals():
    _sys.exit(main())
