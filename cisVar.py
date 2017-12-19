#!/usr/bin/python3

"""
#==============================================================================
#
#          FILE: cisVar (python 3)
#        AUTHOR: Ashley Tehranchi, tehranchi576@gmail.com
#  ORGANIZATION: Stanford University
#       CREATED: 2015-12-12 11:33
# Last modified: 2017-11-06 12:31
#
#
#==============================================================================


cisVar mpileup -F <SampleName> -f <fastaFile> -p <mpileupBEDfile> -B <sortedBam>
cisVar post -F <SampleName> -r <readDepth> -a <allelesFile>
cisVar geno -F <SampleName> -r <readDepth> -i <individualsFile> -g <genotypesFile>
cisVar qtls -F <SampleName> -r <readDepth> -n <numberIndividuals>

"""

######################################################################################
# ARGUEMENTS
######################################################################################


import argparse
import sys
import os
import operator
import subprocess
import scipy
import numpy as np
from scipy import stats
import pandas as pd


parser = argparse.ArgumentParser(description='Find cis QTLs based on an experimental selection method')
parser.add_argument('inputCommand', help='<mpileup>  <post>  <geno>  <qtls> ')
parser.add_argument('-F', '--SampleName', required=False, dest='prefix_name', help='sample/population name')
parser.add_argument('-f', '--fasta', required=False, dest='allCHRfasta', help='fasta file with all chromosomes')
parser.add_argument('-B', '--BAMfile', required=False, dest='sortedBam', help='sorted BAM file')
parser.add_argument('-r', '--readDepth', required=False, type=int, dest='trial_depths', help='minimum read depth per variant')
parser.add_argument('-i', '--individualsFile', required=False, dest='individualslist', help='list of individuals matching genotype matrix; one indv per line')
parser.add_argument('-a', '--allelesFile', required=False, dest='allelesFileName', help='format: chr1    10583   G   A')
parser.add_argument('-n', '--numberIndividuals', required=False, dest='numIndv', help='The number of individuals in the pool')
parser.add_argument('-p', '--mpileupBEDfile', required=False, dest='mpileupBEDfile', help='The  mpileup BED file')
parser.add_argument('-g', '--genoFile', required=False, dest='genosFile', help='The  genotypes file')
args = parser.parse_args()


WorkingDirectory = os.getcwd()

######################################################################################
# MPILEUP FROM BAM FILES USING HG19 MASKED GENOME WITH BLACKLISTED REGIONS REMOVED
######################################################################################

def mpileup(allCHRfasta, mpileupBEDfile, sortedBam, prefix_name):
    os.system("samtools mpileup -Q 25 -f " + allCHRfasta +" -l "+ mpileupBEDfile +" " + sortedBam + " > " + prefix_name + ".mpileup.txt")
    ## outputs prefix.mpileup.txt
    return()


######################################################################################
#POST-CHIP CALCULATION AND FILE TRIMMIMG PLUS HEADER ADDITION
######################################################################################

def postcalc(prefix_name, trial_depths, allelesFileName):
    os.system("module load samtools")
    print(prefix_name)
    print("Post-ChIP calculation for min reads ",trial_depths)
    infilename = prefix_name + ".mpileup.txt"
    inalleles = allelesFileName
    outfilename = prefix_name + "." + str(trial_depths) + ".POST.txt"
    depth = int(trial_depths)

    pileup_dictionary = {}
    count = 1
    bad = 0
    outfile = open(outfilename, 'w')
    with open(infilename, 'r') as pileup:
        for line in pileup:
            #print count
            count = count + 1
            rows = line.rstrip('\n').split("\t")
            if (int(rows[3]) < depth):
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
                if (int(rows[3]) >= depth):
                    read = str(rows[4]).lower() # isolate overlapping positions from reads, make all letters lowercase to count them properly
                    Acount = read.count('a')    # count all A's, A is plus strand, a in negative strand, count both
                    Ccount = read.count('c')
                    Gcount = read.count('g')
                    Tcount = read.count('t')
                    countslist = [Acount, Ccount, Gcount, Tcount]  #make list of all counts
                    index, ALTdepth = max(enumerate(countslist), key=operator.itemgetter(1)) #find index of largest number and the value
                    #print("index ", index)
                    indels = rows[4].count('-') + rows[4].count('+')
                    REFdepth = int(rows[3]) - int(Acount) - int(Ccount) - int(Gcount) - int(Tcount) - int(indels)
                    #print("indels, REFdepth ", indels, REFdepth)
                    ##added 9/2/14
                    #print(rows)
                    rows[3] = str(int(REFdepth) + int(ALTdepth)) ## Adjusts totaldepth to remove third alleles/sequencing errors
                    #print(rows)
                    ##
                    ALTallele = "NA"
                    if (index == 0):
                        ALTallele = 'A'
                    if (index == 1):
                        ALTallele = 'C'
                    if (index == 2):
                        ALTallele = 'G'
                    if (index == 3):
                        ALTallele = 'T'
                    RefFreq = "NA"
                    AltFreq = "NA"
                    if (rows[3] == 0):
                        RefFreq = str(0)
                        print("mistake")
                    else:
                        if (float(rows[3]) > depth) and ((float(REFdepth) / float(rows[3])) >= 0) :    ## depth greater than 0 and postfreq > 0
                            RefFreq = str(REFdepth / float(rows[3]))            # to get ref freq, totaldepth - altdepth = refdepth / totaldepth
                            AltFreq = str(1 - float(RefFreq))
            if not len(rows) == 6:
                bad += 1
                continue

            pileup_dictionary[rows[0] + "." + rows[1]] = (
                rows[2:] + [
                    Acount, Ccount, Gcount, Tcount, ALTdepth,
                    REFdepth, ALTallele, RefFreq, AltFreq
                ]
            )
            #print("ALTallele, row ",ALTallele, rows)
            assert len(pileup_dictionary[rows[0] + "." + rows[1]]) == 13
    linecounter =1
    print('{} mpileup rows were too short and thus skipped'.format(bad))
    with open(inalleles) as alleles:
        for line in alleles:
            row = line.rstrip().split('\t')
            genoRefAllele = row[2]
            genoAltAllele = row[3]
            try: matchingVector = pileup_dictionary[row[0] + '.' + row[1]]
            except: continue
            matchingVector = [str(x) for x in matchingVector]
            if (matchingVector[8] == str(0)): # for alt alleles with 0 reads, force alt allele to be imputed geno alt allele
                matchingVector[10] = genoAltAllele
            if (matchingVector[9] == str(0)): # for ref alleles with 0 reads, force ref allele to be imputed geno ref allele
                matchingVector[0] = genoRefAllele
            postChipRefAllele = matchingVector[0].upper()
            postChipRefFreq = matchingVector[11]
            postChipAltAllele = matchingVector[10].upper()
            postChipAltFreq = matchingVector[12]
            if (postChipRefAllele == genoRefAllele) and (postChipAltAllele == genoAltAllele):
                outAllele = postChipRefAllele
                outfreq = postChipRefFreq
            elif postChipRefAllele == genoAltAllele and (postChipAltAllele == genoRefAllele):
                outAllele = postChipAltAllele
                outfreq = postChipAltFreq
                linecounter = linecounter +1
            else:
                outAllele = 'NA'
                outfreq = 'NA'
            finalout = row[0] + "\t" + row[1] + "\t" + '\t'.join(matchingVector) + "\t" + str(outAllele) + "\t" + str(outfreq) + "\n"
            outfile.write(finalout)
    outfile.close()
    print("post-frequency calculation DONE")
    return ()


######################################################################################
# POST.TXT FILE TRIMMING, ADDITION OF HEADER, BED FILE CREATION
######################################################################################

def postTrim(prefix_name, trial_depths):
    headerOUT = open("mpileup.header.txt", 'w')
    headertemp = 'Chr'+'\t'+'position'+'\t'+'REFallele'+'\t'+'Depth'+'\t'+'Acount'+'\t'+'Ccount'+'\t'+'Gcount'+'\t'+'Tcount'+'\t'+'ALTdepth'+'\t'+'REFDepth'+'\t'+'ALTallele'+'\t'+'REFfreq'+'\t'+'ALTfreq'+'\t'+'POSTallele'+'\t'+'POSTfreq'+'\n'
    headerOUT.write(headertemp)
    headerOUT.close()

    print("File trimming ")
    os.system("sed '/NA/d' " + prefix_name + "." + str(trial_depths) +  ".POST.txt | cut -f-4,7- > " + prefix_name + "." + str(trial_depths) + ".POSTt.txt")
    os.system("cat mpileup.header.txt " + prefix_name + "."  + str(trial_depths) + ".POSTt.txt > " + prefix_name + "." + str(trial_depths) + ".POSTth.txt")
    print(os.system("head " + prefix_name + "." + str(trial_depths) + ".POSTth.txt"))
    postin = prefix_name + "." + str(trial_depths) + ".POSTt.txt"
    bedoutfile = prefix_name + "." + str(trial_depths) + ".bed"
    with open(postin) as openfile, open(bedoutfile, 'w') as bedOUT:
        for line in openfile:
            rowtmp = line.rstrip().split('\t')
            chrom = rowtmp[0]
            snp = rowtmp[1]
            startpos = int(snp) - 1
            OUTrow = str(chrom) + '\t' + str(startpos) + '\t' + str(snp) + '\n'
            bedOUT.write(OUTrow)
    print("File trimming DONE")
    return()

######################################################################################
# GENOTYPE EXTRACTION OF POSTCALC SNPS
######################################################################################

def genoExtract(prefix_name, trial_depths, individualslist, genosFile):
    new_prefix = prefix_name + "." + str(trial_depths)
    postdataIN = new_prefix + ".POSTth.txt"
    outName = new_prefix + ".genotypes.txt"
    print("file: ", new_prefix)

    # Get list of sites with POSTfreq data
    print("Getting SNPs")
    postlocs = pd.read_csv(postdataIN, sep='\t')
    postlocs['idx'] = postlocs[postlocs.columns[0]].astype(str) + '.' + postlocs[postlocs.columns[1]].astype(str)
    print("Got SNPs")

    print("Reading file Locations")
    with open(individualslist , 'r') as indv:
        indv_list = frozenset(indv.read().strip().split('\n'))
    dtypes = {i: int for i in indv_list}
    sidx = {}
    c = 0
    for i in indv_list:
        sidx[i] = c
        c += 1
    print("Individual Locations Read Complete")

    # Load genotypes
    print("Loading genotypes")
    genos = pd.read_csv(genosFile, sep='\t', dtype=dtypes)
    genos['idx'] = genos[genos.columns[0]].astype(str) + '.' + genos[genos.columns[1]].astype(str)
    print("Genotypes loaded")
    print("Dropping indels")
    l1 = len(genos)
    genos = genos[(genos.alt.str.len() == 1) & (genos.ref.str.len() == 1)]
    l2 = len(genos)
    print("Dropped {} indels".format(l1-l2))
    print("Dropping duplicates")
    genos = genos.drop_duplicates('idx', keep=False)
    l3 = len(genos)
    print("Dropped {} duplicates".format(l2-l3))
    print("Filtering by SNP")
    genos = genos[genos.idx.isin(postlocs.idx)]
    postlocs = postlocs[postlocs.idx.isin(genos.idx)]
    print("Sorting")
    genos = genos.sort_values(list(genos.columns[:2]))
    postlocs = postlocs.sort_values(list(postlocs.columns[:2]))
    genos = genos.set_index('idx', drop=True)
    postlocs = postlocs.set_index('idx', drop=True)
    print('Len {}: {}'.format(postdataIN, len(postlocs)))
    print('Len genotypes: {}'.format(len(genos)))
    print("Writing", postdataIN)
    postlocs.to_csv(postdataIN, sep='\t', index=False)
    print("Done")
    print("Writing", prefix_name+"geno.txt")
    genos.to_csv(prefix_name+"geno.txt", sep='\t')
    print("Done")
    print("Filtering by individual")
    geno_len = len(genos)
    genos = genos[list(indv_list)]
    assert len(genos.columns) == len(indv_list)
    genos = genos.dropna()
    final_len = len(genos)
    print('Dropped {} duplicates, writing {} lines to transpose'
          .format(geno_len-final_len, final_len))
    print("Transposing")
    trans = genos.T
    print("Writing", prefix_name+"transpose.txt")
    trans.to_csv(prefix_name+"transpose.txt", sep='\t')
    print("Done")

    print('Final filter')
    trans = trans[trans.index.to_series().isin(indv_list)]
    trans['sort'] = trans.index.to_series().map(sidx)
    trans = trans.sort_values('sort')
    trans = trans.drop('sort', axis=1)
    trans = trans.astype(int)
    print("Writing numeric matrix:", outName)
    trans.to_csv(outName, sep='\t', header=False, index=False)
    print("Done")
    return


######################################################################################
# REGRESSION WITH Z-SCORE, P-VALUE
######################################################################################

def regPVAL(prefix_name, trial_depths, numIndv):
    new_prefix = prefix_name + "." + str(trial_depths)
    postdataIN = new_prefix + ".POSTth.txt"
    genosoutName = new_prefix + ".genotypes.txt"
    ### make sure location of R script is correct
    orig_r = os.path.join(os.path.dirname(__file__), 'regression_qtls.R')
    if not os.path.isfile(orig_r):
        raise OSError('File not found: {}'.format(orig_r))
    commandLine = "sed -e 's!postdataIN!" + postdataIN + "!g' " + orig_r + " | sed -e 's!genosoutName!" + genosoutName + "!g' | sed -e 's!prefix!" + new_prefix + "!g' | sed -e 's!numIndiv!" + numIndv + "!g'> " + prefix_name + "1.R"
    print(commandLine)
    subprocess.check_call(commandLine, shell=True)
    commandout = "Rscript " +prefix_name+ "1.R"
    subprocess.check_call("chmod u+x " +prefix_name+ "1.R", shell=True)
    subprocess.check_call(commandout, shell=True)
    os.system("rm " +prefix_name+ "1.R")
    print("regression DONE")
    return


##################
# COMMANDS #
##################

if args.inputCommand == 'mpileup':
    mpileup(args.allCHRfasta, args.mpileupBEDfile, args.sortedBam, args.prefix_name)
elif args.inputCommand == 'post':
    postcalc(args.prefix_name, args.trial_depths, args.allelesFileName)
    postTrim(args.prefix_name, args.trial_depths)
elif args.inputCommand == 'geno':
    genoExtract(args.prefix_name, args.trial_depths, args.individualslist, args.genosFile)
elif args.inputCommand == 'qtls':
    regPVAL(args.prefix_name, args.trial_depths, args.numIndv)
