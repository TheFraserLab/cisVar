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
import os
import argparse
import operator
import subprocess

import psutil
from datetime import datetime

######################################################################################
# ARGUEMENTS
######################################################################################


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
    """Filter genotypes by SNP and individual.

    First loops through genotypes and filters to a temp file, then transposes,
    then filters transposed lines by individual.

    Takes approximately 10 minutes and uses 10GB of memory on a genotype file
    of 77 million SNPs and an individuals files of 66 individuals.

    Params
    ------
    prefix_name : str
    trial_depths : int
    individualslist : str
        Path to file of individuals separated by line
    genosFile : str
        Path to genotype file to parse. Must be in the format:
            chrom\\tpos\\tref\\talt\\tgenotypes...
        Genotypes must be in 0/1/2 format, and the file must contain a header
        with the individual names.
        Note: Genotype file *must* be sorted by position.
    """
    new_prefix = prefix_name + "." + str(trial_depths)
    postdata = new_prefix + ".POSTth.txt"
    out_name = new_prefix + ".genotypes.txt"
    print("File: ", new_prefix)
    # Track memory usage
    process = psutil.Process(os.getpid())

    print("Memory start: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))
    # Get list of sites with POSTfreq data
    print("Getting SNPs from POST data")
    with open(postdata, 'r') as postin:
        # Drop first line if it is a header (position is not numberic)
        header = postin.readline().strip().split('\t')
        if header[1].isdigit():
            postin.seek(0)
        postlocs = [
            rows[0] + '.' + rows[1] for rows in [
                line.split('\t') for line in postin.readlines()
            ]
        ]
    l = len(postlocs)
    postlocs = set(postlocs)
    # The POST data must be unique per SNP, if it isn't we have a big problem
    # somewhere in the earlier steps
    assert l == len(postlocs)
    print("Got {} SNPs".format(l))
    print("Memory SNPs: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))

    print("Filtering genotypes")
    extras = 0
    indels = 0
    dups   = 0
    print("Time: {}".format(datetime.now()))
    # Geno will have the same number of columns as the POST data has rows but
    # will contain numeric matrix only, one column per SNP, and one row per
    # individual *sorted to match individuals.txt*
    mat = []
    with open(genosFile) as genos_in:
        # Parse header
        header = genos_in.readline().strip().split('\t')
        if header[1].isdigit():
            raise Exception(
                'Genotype file {} has no header, cannot parse individuals'
                .format(genosFile)
            )

        # Check file is sorted
        pos = genos_in.tell()
        c = 1000
        print(
            "Genotype file must be sorted by position,",
            "checking first {0} lines".format(c)
        )
        lns = []
        while c:
            lns.append(int(genos_in.readline().split('\t')[1]))
            c -= 1
        assert lns == sorted(lns)
        print("Genotype file appears properly sorted, parsing.")
        genos_in.seek(pos)

        # Point a local variable at the method to save lookup time
        done_snps = set()
        add_to_done = done_snps.add
        # Parse the file
        mat  = [header[4:]]
        for line in genos_in:
            row = line.strip().split('\t')
            # Filter by SNP
            snp = row[0] + '.' + row[1]
            if snp not in postlocs:
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
            mat.append(row[4:])
            # Done
            add_to_done(snp)
    print("Genotype filtering complete")
    print("Memory genos: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))
    print("Time: {}".format(datetime.now()))
    print("Dropped {} unneeded SNPs".format(extras))
    print("Dropped {} indels".format(indels))
    print("Dropped {} duplicates".format(dups))

    print("Checking all POST SNPs matched")
    if len(postlocs) != len(done_snps):
        err_file = new_prefix + ".missing.snps"
        with open(err_file, "w") as fout:
            fout.write(
                '\n'.join(sorted([i for i in postlocs if i not in done]))
            )
        raise Exception(
            "{} SNPs in POST were not in the genotype file, written to {}"
            .format(len(postlocs)-len(done), err_file)
        )
    print("Done")

    print("Transposing to dictionary")
    trans = {i[0]: i[1:] for i in zip(*mat)}
    with open(out_name, 'w') as fout:
        fout.write(
            '\n'.join(
                ['\t'.join(i) for i in zip(*mat)]
            )
        )
    print("Memory transpose: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))
    del mat
    print("Memory transpose tidy: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))
    print("Done")

    print("Filtering by individual and writing")
    with open(individualslist , 'r') as indv, open(out_name, 'w') as fout:
        for ind in indv:
            ind = ind.strip()
            if ind in trans:
                fout.write('\t'.join(trans[ind]) + '\n')
            else:
                print('Missing:', ind)
    print("Memory End: {:.2f}GB".format(process.memory_info().rss/1024/1024/1024))


######################################################################################
# REGRESSION WITH Z-SCORE, P-VALUE
######################################################################################

def regPVAL(prefix_name, trial_depths, numIndv):
    new_prefix = prefix_name + "." + str(trial_depths)
    postdata = new_prefix + ".POSTth.txt"
    genosout_name = new_prefix + ".genotypes.txt"
    ### make sure location of R script is correct
    r_script = os.path.join(os.path.dirname(__file__), 'regression_qtls.R')
    if not os.path.isfile(r_script):
        raise OSError('File not found: {}'.format(r_script))
    subprocess.check_call(['Rscript', r_script, postdata, genosout_name, new_prefix, numIndv])
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
