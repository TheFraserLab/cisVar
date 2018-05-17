#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cisVar: Find cis QTLs based on an experimental selection method

Ashley Tehranchi <tehranchi576@gmail.com>

Stanford University

Version: 2.0.0b1
Created: 2015-12-12
Updated: 2018-05-16

Example usage:
cisVar prep -F test_new -i individuals.txt.gz --chrom-format chr
cisVar mpileup -F <SampleName> -f <fastaFile> -B <sortedBam>
cisVar post -F <SampleName> -r <readDepth>
cisVar geno -F <SampleName> -r <readDepth> -i <individualsFile>
cisVar qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>
cisVar tidy -F <SampleName> -r <readDepth> -p out.dataframe.pandas -t out.dataframe.txt

Note:
The qtls regression step will use approximately 32GB of memory on an averaged-
sized dataset.

The geno step will use approximately 20GB of memory on the same dataset.

"""
from __future__ import print_function
import os
import sys
import bz2
import gzip
import argparse
import operator
import subprocess
import multiprocessing as mp

from time import sleep
from datetime import datetime

from collections import defaultdict

# psutil is used to print memory usage if installed, otherwise ignored
try:
    import psutil
except ImportError:
    psutil = None

CORES = mp.cpu_count()

__version__ = '2.0.0b2'

###########################################################################
#                 Create Necessary Input Files from VCFs                  #
###########################################################################


def prep_files(vcf_files, prefix_name, inds=None, limit_file=None,
               alleles_file=None, chrom_format=0, skip_indels=True, cores=4):
    """Create genotype, bed, and allele files for pipeline.

    Bed is used by mpileup to define SNPs to pileup.
    Alleles file is used to set default ref and alt in POST.

    Note, this only needs to be created once per experiment, not once per
    group/pool. The same files can be used for each population/group/pool,
    they are further filtered in the genotype step.

    This code first loads all locations from the limit file (if present) and
    any alleles to override. These two lists are kept in memory. It then loops
    through each vcf_file in parallel, line-by-line, writing temp files, before
    concatenating the results into a single file.

    Params
    ------
    vcf_files : list
        A list of VCF files to parse, should contain genotypes, ref and alt.
    prefix_name : str
        A name for the output files.
    inds : list or path, optional
        Either a list of individuals or a list of paths to individuals files.
        Should include all individuals in this experiment, in all pools.
    limit_file : path, optional
        Path to a bed/vcf/simple file of SNPs to consider. Should be a subset
        of SNPs in the VCF files. BED3 format is fine, simple format is just
        chromosome and position (1-based) to exclude, tab separated.
    alleles_file : path, optional
        Path to a VCF/bed/txt file with alleles, it text, should be
        chr position (base-1). Can be same file as limit file if desired.
    chrom_format : {0, 1, 2}, optional
        0 = leave
        1 = force chr#
        2 = force #
    skip_indels : bool, optional
        Don't include rows with indels
    cores : int, optional
        Number of parallel cores to use on multiple vcf files.

    Writes
    ------
    prefix.genotypes.txt.gz
        A tab separated file, including header with individuals, of:
        chrom, position (base-1), ref (upper), alt (upper), genotypes (0/1/2)
    prefix.locations.bed.gz
        A BED file of all of the locations to be tested, including names
    prefix.individuals.txt.gz
        A newline separated file of individuals, replicates inds if provided,
        overwrites even if already exists.
    """
    if isinstance(vcf_files, str):
        vcf_files = [vcf_files]
    for fl in vcf_files:
        if not os.path.isfile(fl):
            raise OSError('Could not find file {}'.format(fl))
    assert isinstance(prefix_name, str)
    par = False if len(vcf_files) == 1 or cores < 2 else True
    if limit_file:
        sys.stderr.write(
            "Loading test locations from {}\n".format(limit_file)
        )
        limit_locs = parse_loc_file(limit_file)
    else:
        limit_locs = None
    if alleles_file:
        sys.stderr.write(
            "Loading alleles from {}\n".format(alleles_file)
        )
        alleles = parse_loc_file(alleles_file, mode='alleles')
    else:
        alleles = None
    if par:
        pool = mp.Pool(cores)

    # Individuals
    if os.path.isfile(inds):
        with open_zipped(inds) as fl:
            inds = fl.read().strip().split('\n')
    if not inds:
        inds = geno_file(vcf_files[0], get_inds=True)
    final_inds = list(inds)

    # Write the final individuals, order preserved
    ind_output = '{}.individuals.txt.gz'.format(prefix_name)
    with open_zipped(ind_output, 'w') as fout:
        fout.write('\n'.join(final_inds))

    # Primary parsing
    sys.stderr.write("Begining genotype file parse\n")

    # Header
    header = geno_file(vcf_files[0], get_header=True)

    count = 0
    jobs = []
    geno_output = '{}.genotypes.txt.gz'.format(prefix_name)
    bed_output = '{}.locations.bed.gz'.format(prefix_name)
    for fl in vcf_files:
        ofl = '{}.{}.gz'.format(geno_output, count)
        bfl = '{}.{}.gz'.format(bed_output, count)
        args = (
            fl, ofl, header, bfl, final_inds, limit_locs,
            alleles, skip_indels, chrom_format
        )
        if par:
            jobs.append([ofl, bfl, pool.apply_async(parse_file, args)])
        else:
            parse_file(*args)
            jobs.append((ofl, bfl, None))
        count += 1

    ofls = [i[0] for i in jobs]
    bfls = [i[1] for i in jobs]
    jobs = [i[2] for i in jobs]

    if par:
        for job in jobs:
            job.wait()
        sleep(2)

    for ofl in ofls:
        done_file = ofl + '.done'
        if not os.path.isfile(done_file):
            sleep(1)
            if not os.path.isfile(done_file):
                print('Cannot find {}, assuming failure'.format(done_file))
                raise Exception('Missing done file')
        os.remove(done_file)

    for bfl in bfls:
        done_file = bfl + '.done'
        if not os.path.isfile(done_file):
            sleep(1)
            if not os.path.isfile(done_file):
                print('Cannot find {}, assuming failure'.format(done_file))
                raise Exception('Missing done file')
        os.remove(done_file)

    # Merge
    print('Merging files')
    with open_zipped(geno_output, 'w') as fout:
        fout.write(
            'chrom\tposition\tref\talt\t{}\n'.format('\t'.join(final_inds))
        )
    subprocess.check_call(
        'zcat {} | gzip >> {}'.format(' '.join(ofls), geno_output),
        shell=True
    )
    subprocess.check_call(
        'zcat {} | gzip > {}'.format(' '.join(bfls), bed_output), shell=True
    )

    print('Removing tmp files')
    for temp_file in ofls + bfls:
        if os.path.isfile(temp_file):
            os.remove(temp_file)


def parse_file(geno, outfile, header, bedfile, inds=None, limit_locs=None,
               alleles=None, skip_indels=False, chrom_format=0, log=None,
               fail_if_too_few_inds=True):
    """File parser, see parse_files for information."""
    if not log:
        log = outfile + '.parse.log'
    kept = 0
    indels = 0
    skipped = 0
    superbad = []
    print('Creating {} and {} from {}'.format(outfile, bedfile, geno))
    if chrom_format == 0:
        def chromit(x):
            """Return x."""
            return x
    elif chrom_format == 1:
        def chromit(x):
            """Force x to start with chr."""
            return x if x.startswith('chr') else 'chr' + x
    elif chrom_format == 2:
        def chromit(x):
            """Force x to not start with chr."""
            return x if not x.startswith('chr') else x[3:]
    else:
        raise ValueError('chrom_format must be one of 0, 1, 2')

    if not limit_locs:
        limit_locs = dict()
    if not alleles:
        alleles = dict()

    lg = open_zipped(log, 'a')
    with open_zipped(outfile, 'w') as gn, open_zipped(bedfile, 'w') as bd:
        lg.write('Parse initiated\n')
        for f in geno_file(geno, check_header=header, ind_list=inds, log=lg):
            # Write core records
            if skip_indels and (len(f[2]) > 1 or len(f[3]) > 1):
                indels += 1
                continue
            chrom = chromit(f[0])
            pos = int(f[1])
            # Drop excluded positions
            try:
                if pos not in limit_locs[chrom]:
                    skipped += 1
                    continue
            except KeyError:
                pass

            # Override alleles
            try:
                if pos in alleles[chrom]:
                    ref, alt = alleles[chrom][pos]
                else:
                    ref = f[2]
                    alt = f[3]
            except KeyError:
                ref = f[2]
                alt = f[3]

            name = f[4]
            gt_inds = f[5:]

            outstr = (
                '{chr}\t{pos}\t{ref}\t{alt}'.format(
                    chr=chrom, pos=pos, ref=ref, alt=alt
                )
            )
            # Write genotypes
            if inds:
                gt_ind_l = len(gt_inds)
                ind_l = len(inds)
                if gt_ind_l < ind_l:
                    lg.write(
                        '{chr}.{pos} had too few ({ind}) individuals\n'
                        .format(chr=chrom, pos=f[1], ind=len(gt_inds))
                    )
                    if fail_if_too_few_inds:
                        superbad.append(f)
                    continue
                if gt_ind_l > ind_l:
                    lg.write(
                        '{chr}.{pos} had too many ({ind}) individuals\n'
                        .format(chr=chrom, pos=f[1], ind=len(gt_inds))
                    )
                    superbad.append(f)
                    continue
            gn.write('{}\t{}\n'.format(outstr, '\t'.join(gt_inds)))
            bd.write('{}\t{}\t{}\t{}\n'.format(chrom, pos-1, pos, name))
            kept += 1
    lg.write(
        'Skipped: {}\nIndels: {}\nKept: {}\nParse Complete\n'
        .format(skipped, indels, kept)
    )
    print(
        'Finished {}. Kept {} SNPs. Dropped {} indels and skipped {}'
        .format(outfile, kept, indels, skipped)
    )
    if superbad:
        msg = 'Length failure despite length filtering:\n{}\n'.format(superbad)
        lg.write(msg)
        print(msg)

    print('Writing done files')
    with open(outfile + '.done', 'w') as fout:
        fout.write('Done\n')
    with open(bedfile + '.done', 'w') as fout:
        fout.write('Done\n')

    lg.close()
    return


############################################################################
# mpileup BAM files with hg19 masked genome with blacklist regions removed #
############################################################################


def mpileup(allCHRfasta, sortedBam, prefix_name, mpileupBEDfile=None):
    """Run mpileup on the provided files.

    Params
    ------
    allCHRfasta : str
        Path to genome wide fasta file
    sortedBam : str
        Path to BAM file sorted by position
    prefix_name : str
        Prefix to use for default file names
    mpileupBEDfile : str, optional
        Path to BED file for mpileup, defaults to file output by prep.

    Writes
    ------
    prefix_name.mpileup.txt.gz
    """
    if not mpileupBEDfile:
        mpileupBEDfile = '{}.locations.bed.gz'.format(prefix_name)
    cmnd = (
        "samtools mpileup -Q 25 -f {} -l {} {} | gzip > {}.mpileup.txt.gz"
        .format(allCHRfasta, mpileupBEDfile, sortedBam, prefix_name)
    )
    sys.stderr.write('Running {}\n'.format(cmnd))
    subprocess.check_call(cmnd, shell=True)


###########################################################################
#                 POST-ChIP Calculation and File Trimming                 #
###########################################################################


def postcalc(prefix_name, inalleles, trial_depths):
    """Calculate POST frequencies from mpileup results."""
    print("Post-ChIP calculation for min reads {}".format(trial_depths))

    infilename  = prefix_name + ".mpileup.txt.gz"
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
            row = line.rstrip().split('\t')[:4]
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


def regPVAL(prefix_name, trial_depths, numIndv=None):
    """Run the actual regression R code."""
    new_prefix = prefix_name + "." + str(trial_depths)
    postdata = new_prefix + ".POSTth.txt"
    genosout_name = new_prefix + ".genotypes.txt"
    if not numIndv:
        wc = subprocess.check_output(['wc', '-l', genosout_name])
        numIndv = wc.decode().split(' ')[0]
    numIndv = str(int(numIndv))

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
        return gzip.open(infile, mode)
    if infile.endswith('.bz2'):
        if hasattr(bz2, 'open'):
            return bz2.open(infile, mode)
        return bz2.bz2file(infile, mode)

    # Fall back on regular open
    return open(infile, mode)


# Simple location file handling

def parse_loc_file(infile, mode='locations'):
    """Return a dict of SNPs from a bed/vcf/text of chr postition (base-1)."""
    if mode not in ['locations', 'alleles', 'names']:
        raise ValueError(
            'mode must be one of locations, alleles, names. Is {}'
            .format(mode)
        )
    names = infile.split('.')
    if 'vcf' in names:
        fl_tp = 'vcf'
    elif 'bed' in names:
        fl_tp = 'bed'
    else:
        fl_tp = 'txt'
    splt_char = '\t ' if fl_tp == 'txt' else '\t'
    with open_zipped(infile) as fl:
        for line in fl:
            loc = fl.tell()
            # Comment lines
            if line.startswith('#'):
                loc = fl.tell()
            # First non-comment
            else:
                header = line.strip().split(splt_char)
                if line[1].isdigit():
                    # Not a header
                    fl.seek(loc)
                break
        # All non-head lines
        if mode == 'locations':
            res = _get_locs(fl, fl_tp, splt_char)
        elif mode == 'alleles':
            res = _get_alleles(fl, fl_tp, splt_char)
        elif mode == 'names':
            res = _get_names(fl, fl_tp, splt_char)
        else:
            raise ValueError('Invalid mode')
        # Generalize chromosome name
        for k, v in res.items():
            res[flip_chrom(k)] = v
        return res


def _get_locs(fl, tp, sp):
    result = defaultdict(set)
    for line in fl:
        f = line.rstrip().split(sp)
        if tp == 'bed':
            loc = int(f[2])
        else:
            loc = int(f[1])
        result[f[0]].add(loc)
    return result


def _get_alleles(fl, tp, sp):
    if tp == 'bed':
        raise ValueError("Can't get alleles from a bed file")
    result = defaultdict(dict)
    for line in fl:
        f = line.rstrip().split(sp)
        loc = int(f[1])
        if tp == 'vcf':
            ref = f[3]
            alt = f[4]
        else:
            ref = f[2]
            alt = f[3]
        result[f[0]][loc] = (ref, alt)
    return result


def _get_names(fl, tp, sp):
    result = defaultdict(dict)
    for line in fl:
        f = line.rstrip().split(sp)
        if tp == 'bed':
            name = f[3]
            loc = int(f[2])
        else:
            name = f[2]
            loc = int(f[1])
        result[f[0]][loc] = name
    return result


def flip_chrom(chrom):
    """Return opposite of chr# or # for chromosome name."""
    if isinstance(chrom, int):
        chrom = str(chrom)
    if chrom.startswith('chr'):
        return chrom[3:]
    return 'chr' + chrom


# Genotype file handling


def geno_file(infile, check_header=None, ind_list=None,
              get_header=False, get_inds=False, log=sys.stderr):
    """Wrapper for genotype files, currently only vcf

    Parameters
    ----------
    infile : str
        Path to a vcf file, can be zipped.
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
    name : str
    inds : list of str
    """
    name = infile.split('.')
    if 'vcf' in name:
        return vcf_file_gt(
            infile, check_header, ind_list, get_header, get_inds, log
        )
    else:
        raise NotImplementedError(
            'Currently only VCFs supported, filetype required in name'
        )


def vcf_file_gt(vcf, check_header=None, ind_list=None,
                get_header=False, get_inds=False, log=sys.stderr):
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
    name : str
    inds : list of str
    """
    # Get header at genotype loc
    with open_zipped(vcf) as fin:
        line = fin.readline()
        while True:
            if not line.startswith('#'):
                break
            header = line
            pos    = fin.tell()
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
    total_inds = len(inds)
    if ind_list:
        ind_count = len(ind_list)
        if sorted(ind_list) != sorted(inds):
            for ind in ind_list:
                ind_pos.append(inds.index(ind))
    else:
        ind_count = len(inds)

    short_lines = 0
    bad_gts = 0
    count = 0

    # Actually parse file
    with open_zipped(vcf) as fin:
        fin.seek(pos)
        count = 0
        for line in fin:
            count += 1
            f = line.strip().split('\t')
            out = [f[0], int(f[1]), f[3], f[4], f[2]]
            gti = f[8].split(':').index('GT')
            ind_data = f[9:]
            if not len(ind_data) == total_inds:
                short_lines += 1
                continue
            if ind_pos:
                these_inds = []
                for loc in ind_pos:
                    these_inds.append(ind_data[loc])
                if len(these_inds) != ind_count:
                    short_lines += 1
                    continue
            else:
                these_inds = ind_data
            for ind in these_inds:
                gt = parse_gt(ind.split(':')[gti])
                if not gt:
                    break
                out.append(gt)
            if (len(out) - 5) != ind_count:
                #  sys.stderr.write(
                    #  'Bad line: {}\t{}\t{}\nLen: {}\tExpected: {}\n'.format(
                        #  count, f[0], f[1], len(out)-4, ind_count
                    #  )
                #  )
                #  raise Exception('Bugger')
                bad_gts += 1
                continue
            yield out
    log.write(
        'Total: {}\nBad Genotypes: {}\nShort Lines (not enough inds): {}\n'
        .format(count, bad_gts, short_lines)
    )
    print(
        '{}: Total: {}, Bad Genotypes: {}, Short Lines (not enough inds): {}'
        .format(vcf, count, bad_gts, short_lines)
    )
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
    elif '/' in gt:
        if gt == '0/0':
            gt = '0'
        elif gt in ['0/1', '1/0']:
            gt = '1'
        elif gt == '1/1':
            gt = '2'
        else:
            gt = None
    elif gt.isdigit():
        gt = gt
    else:
        gt = None
    return gt



###########################################################################
#                          Command Line Parsing                           #
###########################################################################


def main(argv=None):
    """Run as a script, parse command line args."""
    if not argv:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Subcommands
    modes = parser.add_subparsers(
        dest='inputCommand',
        metavar='{prep,mpileup,post,geno,qtls,tidy,get_snake}'
    )

    # Shared commands
    shared_cmds = argparse.ArgumentParser(add_help=False)
    run_opts = shared_cmds.add_argument_group("Run Options")
    run_opts.add_argument(
        '-F', '--SampleName', dest='prefix_name', default='cis_var',
        help='sample/population name (default: cis_var)'
    )
    run_opts.add_argument(
        '-r', '--readDepth', type=int, dest='trial_depths', default=20,
        help='minimum read depth per variant (default: 20)'
    )

    # prep
    prep_cmd_grp = modes.add_parser(
        'prep', description="Prepare genotype files",
        parents=[shared_cmds], help="Prepare genotype files",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    prep_cmd = prep_cmd_grp.add_argument_group(
        "Prep Options"
    )
    prep_cmd.add_argument(
        '-i', '--all-inds',
        help='File of individuals in all groups, one per line'
    )
    prep_cmd.add_argument(
        '-l', '--limit-file',
        help='BED/VCF/txt file of SNPs to consider'
    )
    prep_cmd.add_argument(
        '-a', '--allele-file',
        help='BED/VCF/txt file of alleles to override VCF allels (subset of vcf)'
    )
    prep_cmd.add_argument(
        '--chrom-format', choices={'chr', 'num', 'ignore'}, default='ignore',
        help='chr: make format "chr#", num: make format "#", ' +
        'ignore: do nothing (default: ignore)'
    )
    prep_cmd.add_argument(
        '--include-indels', action='store_false',
        help='Do not skip indels'
    )
    prep_cmd.add_argument(
        '-c', '--cores', type=int, default=CORES,
        help='Number of cores to use (default: all)'
    )
    prep_cmd.add_argument(
        'vcf_files', nargs='+', help='VCF files with genotypes'
    )


    # mpileup
    mpileup_cmd_grp = modes.add_parser(
        'mpileup', description="Run mpileup",
        parents=[shared_cmds], help="Run mpileup",
        aliases=['m'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mpileup_cmd = mpileup_cmd_grp.add_argument_group(
        "mpileup Options"
    )
    mpileup_cmd.add_argument(
        '-f', '--fasta', required=True, dest='allCHRfasta',
        help='fasta file with all chromosomes (Required)'
    )
    mpileup_cmd.add_argument(
        '-B', '--BAMfile', required=True, dest='sortedBam',
        help='sorted BAM file (Required)'
    )
    mpileup_cmd.add_argument(
        '-p', '--mpileupBEDfile', required=False, dest='mpileupBEDfile',
        help='BED to use instead of the BED generated in the prep phase ' +
        '(Do not use if possible, use prep with limit instead)'
    )

    # POST
    post_cmd_grp = modes.add_parser(
        'post', description="Run POST frequency calculation",
        parents=[shared_cmds], help="Run POST frequency calculation",
        aliases=['p'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    post_cmd = post_cmd_grp.add_argument_group("POST Options (Deprecated)")
    post_cmd.add_argument(
        '-a', '--allelesFile', required=False, dest='genosFile',
        help='The genotypes file, ' +
        '(Optional, default is file created in prep)'
    )

    # geno
    geno_cmd_grp = modes.add_parser(
        'geno', description="Munge genotypes to prepare for regression",
        parents=[shared_cmds],
        help="Munge genotypes to prepare for regression", aliases=['g'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    geno_cmd = geno_cmd_grp.add_argument_group("Genotype Options")
    geno_cmd.add_argument(
        '-g', '--genoFile', required=False, dest='genosFile',
        help='The genotypes file, ' +
        '(Optional, default is file created in prep)'
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
    qtls_cmd = qtls_cmd_grp.add_argument_group("Regression Options")
    qtls_cmd.add_argument(
        '-n', '--numberIndividuals', dest='numIndv',
        help='The number of individuals in the pool ' +
        '(if omitted, calculated from genotype file length)'
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
    outputs.add_argument('-t', '--textfile', help="Parsed output")
    outputs.add_argument('-p', '--pandasfile', help="Parsed dataframe")

    # Download pipeline
    snk = modes.add_parser(
        'get_snake',
        description="Download Snakefile and config to current dir",
        help="Download Snakefile and config", aliases=['download'],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    snk.add_argument('-x', '--extra', action='store_true',
                     help='Get additional sample and cluster configs')

    args = parser.parse_args(argv)

    if args.inputCommand in ['prep']:
        if args.chrom_format == 'ignore':
            chrom_format = 0
        elif args.chrom_format == 'chr':
            chrom_format = 1
        elif args.chrom_format == 'num':
            chrom_format = 2
        prep_files(
            args.vcf_files, args.prefix_name, inds=args.all_inds,
            limit_file=args.limit_file, alleles_file=args.allele_file,
            chrom_format=chrom_format, skip_indels=args.include_indels,
            cores=args.cores
        )

    elif args.inputCommand in ['mpileup', 'm']:
        mpileup(
            args.allCHRfasta, args.sortedBam, args.prefix_name,
            mpileupBEDfile=args.mpileupBEDfile,
        )

    elif args.inputCommand in ['post', 'p']:
        alleles_file = args.genosFile
        if not alleles_file or not os.path.isfile(alleles_file):
            alleles_file = '{}.genotypes.txt.gz'.format(args.prefix_name)
        postcalc(args.prefix_name, alleles_file, args.trial_depths)
        postTrim(args.prefix_name, args.trial_depths)

    elif args.inputCommand in ['geno', 'g']:
        if args.genosFile:
            genosFile = args.genosFile
        else:
            genosFile = '{}.genotypes.txt.gz'.format(args.prefix_name)
        genoExtract(
            args.prefix_name, args.trial_depths,
            genosFile, args.individualslist
        )

    elif args.inputCommand in ['qtls', 'q', 'regression', 'r']:
        regPVAL(args.prefix_name, args.trial_depths, args.numIndv)

    elif args.inputCommand in ['tidy', 't']:
        prefix = args.prefix_name + '.' + str(args.trial_depths)
        ifl = prefix + '.total.txt'
        ofl = args.textfile if args.textfile else '{}.final.txt'.format(prefix)
        pfl = args.pandasfile if args.pandasfile else '{}.final.pd'.format(prefix)
        bedfile = args.bedfile
        if not args.bedfile:
            bfl = '{}.locations.bed.gz'.format(prefix)
            if os.path.isfile(bfl):
                bedfile = bfl
        parse_regression(
            ifl, textfile=ofl, pandasfile=pfl, bed=bedfile
        )

    elif args.inputCommand in ['get_snake', 'download']:
        try:
            from wget import download
        except ImportError:
            sys.stderr.write(
                'Cannot get files with wget installed, '
                'try pip install --user wget\n'
            )
            return 1
        url = 'https://raw.githubusercontent.com/TheFraserLab/cisVar/master/pipeline/{0}'
        files = ['Snakefile', 'cisvar_config.json']
        if args.extra:
            files += ['cisvar_samples.json', 'cluster.json']
        for fl in files:
            download(url.format(fl))
        print(
            'Done, edit cisvar_config.json for your needs, use README on '
            'github repo for help'
        )

    else:
        sys.stderr.write('Unrecognized command {}'.format(args.inputCommand))
        parser.print_help(file=sys.stderr)
        return 2


if __name__ == '__main__' and '__file__' in globals():
    sys.exit(main())
