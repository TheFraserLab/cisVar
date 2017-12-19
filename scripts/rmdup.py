#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script removes duplicates and filters out unmapped reads

Written by Eilon Sharon:
    `https://github.com/eilon-s/bioinfo_scripts/blob/master/rmdup.py`_
"""

import os
import sys
import random
import argparse

import pysam
import psutil


class SlimRead:

    """
    Simple class to hold read information.
    """

    def __init__(self, read, line_num):
        self.line_num = line_num
        self.tid = read.tid
        self.reference_start = read.reference_start
        self.reference_end = read.reference_end


def memory_usage_psutil():
    """Return the memory usage in MB."""
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0]/float(2 ** 20)
    return mem


def rmdup_pe(in_bam_file_name, out_bam_file_name):
    """Paired end duplicate removal."""

    sys.stderr.write("Running in paired-ends mode\n")
    in_bamfile = pysam.AlignmentFile(in_bam_file_name,"rb")

    sys.stderr.write("#refrences:%d\n" % in_bamfile.nreferences)

    # build hash of read pairs and uniq position reads pairs
    reads_dict = {}
    uniq_reads_corr_dict = {}
    r=0
    for read in in_bamfile:
        r += 1
        if (r % 1000000 == 0):
        #sys.stderr.write("parsing read %d\n" %r)
            sys.stderr.write("parsing read %d, memory: %dMB\n" % (r,memory_usage_psutil()))

        if read.is_paired:
        #sys.stderr.write("read id: %s " % read.qname)
            cur_read_key = read.qname.__hash__()

            # skip over unmapped reads using FLAG
            mapped_flags = [99,147,83,163,67,131,115,179,81,161,97,145,65,129,113,177]
            if read.flag in mapped_flags:

            # pair already in reads_dict
                if cur_read_key in reads_dict:
                    # already one read in the list
                    if len(reads_dict[cur_read_key]) == 1:
                    # saving the line of the read for memory efficiency
                        reads_dict[cur_read_key].append(SlimRead(read,r))

                        # found a pair -> add id to the uniq set
                        cur_first_read = reads_dict[cur_read_key][0]
                        cur_second_read = reads_dict[cur_read_key][1]

                        #print("found a pair\n%s\n%s\n" % (cur_first_read, cur_second_read) )
                        #sys.stderr.write("found a pair: %s" % read.qname)

                        if (cur_first_read.tid != cur_first_read.tid):
                            sys.stderr.write(
                                "Removing pairs that are mapped to two "
                                "different reference sequences (chromosomes)"
                            )
                            continue

                        # converting reference id to chr for the output to be sorted
                        cur_pair_key_str = (str(chr(65+cur_first_read.tid)) + "_" + str(chr(65+cur_second_read.tid)) + "_" +
                        str(min(cur_first_read.reference_start,cur_first_read.reference_end, cur_second_read.reference_start,cur_second_read.reference_end)) + "_" +
                        str(max(cur_first_read.reference_start,cur_first_read.reference_end, cur_second_read.reference_start,cur_second_read.reference_end)))

                        cur_pair_key = cur_pair_key_str.__hash__()

                        #sys.stderr.write("%s\n" % cur_pair_key)

                        # replacing the dictionary entry by just the line numbers
                        reads_dict[cur_read_key] = [cur_first_read.line_num, cur_second_read.line_num]

                        #sys.stderr.write("Read pair position key: %s" % cur_pair_key)

                        #if (cur_pair_key == "0_0_1472172_1472342"):
                        #  print("found a pair\n%s\n%s\n" % (cur_first_read, cur_second_read) )

                        # adding to the dictionary of unique pairs
                        if cur_pair_key in uniq_reads_corr_dict:
                        #sys.stderr.write("-----in dict: %s, %s" % (cur_pair_key,read.qname) )
                            uniq_reads_corr_dict[cur_pair_key].append(cur_read_key)

                        else:
                        #sys.stderr.write("not in dict: %s, %s" % (cur_pair_key,read.qname))
                            uniq_reads_corr_dict[cur_pair_key] = [cur_read_key]

                        #
                        #uniq_reads_corr_dict[]
                    else: # more than 2 rreads with the same id? error!
                        raise Exception("More than two reads with the same id: %s\n" % read)
                else: #adding the first read in a pair
                    #sys.stderr.write("Adding first read: %s\n" % read.qname)
                    # saving the line of the read for memory efficiency
                    reads_dict[cur_read_key] = [SlimRead(read,r)]

            else:
                continue

        else:
            sys.stderr.write("Found a single end read: %s\n" % read)

    in_bamfile.close()

    # implemented for pair ends
    sys.stderr.write(
        "Found {} uniq reads pairs ({} reads), memory: {}MB\n".format(
            len(uniq_reads_corr_dict), len(uniq_reads_corr_dict)*2,
            memory_usage_psutil()
        )
    )
    sys.stderr.write("Looking for duplicates...\n")

    max_num_of_dupliactes = 0
    cnt_dup_pos = 0

    out_reads_line_numbers = []

    for uniq_pos, pairs_ids in sorted(uniq_reads_corr_dict.items()):
        cur_num_pairs = len(pairs_ids)
        if (cur_num_pairs>1):
            cnt_dup_pos += 1
            max_num_of_dupliactes = max(max_num_of_dupliactes,cur_num_pairs)
            #print("duplicated position %s, #pairs: %d\n" % (uniq_pos, cur_num_pairs))
        #print("position %s, #pairs: %d\n" % (uniq_pos, cur_num_pairs))
        cur_pair_id = pairs_ids[random.randint(0,cur_num_pairs-1)]
        cur_pair = reads_dict[cur_pair_id]

        # writing to the output file
        #out_bamfile.write(cur_pair[0])
        #out_bamfile.write(cur_pair[1])

        # choosing which lines to print
        out_reads_line_numbers.append(cur_pair[0])
        out_reads_line_numbers.append(cur_pair[1])

    # sorting the reads
    out_reads_line_numbers.sort()

    sys.stderr.write(
        "Found {} duplicated segments in {} uniq segments ({:.2f}), max  #duplictaes: {}"
        .format(cnt_dup_pos, len(uniq_reads_corr_dict),
               (1.0*cnt_dup_pos)/len(uniq_reads_corr_dict),
               max_num_of_dupliactes)
    )

    # printing the output bam
    in_bamfile  = pysam.AlignmentFile(in_bam_file_name,"rb")
    out_bamfile = pysam.AlignmentFile(out_bam_file_name, "wb", template=in_bamfile)

    sys.stderr.write("Writing output bam file...\n")

    out_ind=0
    r=0
    for read in in_bamfile:
        r += 1
        if (r<out_reads_line_numbers[out_ind]):
            next
        elif (r==out_reads_line_numbers[out_ind]):
            out_bamfile.write(read)
            out_ind += 1
            if (out_ind >= len(out_reads_line_numbers)):
                break
        else:
            sys.stderr.write(
                "Error bug in printing the out file {},{},{}\n".format(
                    r, out_ind, out_reads_line_numbers[out_ind]
                )
            )
            sys.exit(1)

        if (r % 1000000 == 0):
            sys.stderr.write(
                "writing. going over read {}, printed {} reads memory: {}MB\n"
                .format(
                    r, out_ind, memory_usage_psutil()
                )
            )

    out_bamfile.close()
    in_bamfile.close()

    # writing a single random pair from each unique position
    sys.stderr.write("Finish writing to output file: %s" % (out_bam_file_name) )


def rmdup_se(in_bam_file_name, out_bam_file_name):
    """Single eng duplicate removal."""

    sys.stderr.write("Running in single end mode\n")
    in_bamfile  = pysam.AlignmentFile(in_bam_file_name,"rb")
    out_bamfile = pysam.AlignmentFile(out_bam_file_name, "wb", template=in_bamfile)

    # build hash of read pairs and uniq position reads pairs
    reads_dict = {}
    uniq_reads_corr_dict = {}
    r=0
    for read in in_bamfile:
        r += 1
        if (r % 10000 == 0):
            sys.stderr.write("parsing read %d, memroy: %dMB\n" % (r,memory_usage_psutil()))

        # adding to dict
        reads_dict[read.qname] = read

        # converting reference id to chr for the ourput to be sorted
        cur_pair_key = str(chr(65+read.tid)) + "_" + str(read.tid) + "_" + str(min(read.reference_start,read.reference_end)) + "_" + str(max(read.reference_start,read.reference_end))

        # adding to the dictionary of unique pairs
        if cur_pair_key in uniq_reads_corr_dict:
            #sys.stderr.write("-----in dict: %s, %s" % (cur_pair_key,read.qname) )
            uniq_reads_corr_dict[cur_pair_key].append(read.qname)
        else:
            #sys.stderr.write("not in dict: %s, %s" % (cur_pair_key,read.qname))
            uniq_reads_corr_dict[cur_pair_key] = [read.qname]

    sys.stderr.write("Found %d uniq positions for %d reads\n" % (len(uniq_reads_corr_dict),r))
    sys.stderr.write("Current memory: %d MB\n" % memory_usage_psutil())

    max_num_of_dupliactes = 0
    cnt_dup_pos = 0
    for uniq_pos,read_ids in sorted(uniq_reads_corr_dict.items()):
        cur_num_reads = len(read_ids)
        if (cur_num_reads>1):
            cnt_dup_pos += 1
            max_num_of_dupliactes = max(max_num_of_dupliactes,cur_num_reads)
            #print("duplicated position %s, #pairs: %d\n" % (uniq_pos, cur_num_pairs))
        #print("position %s, #pairs: %d\n" % (uniq_pos, cur_num_pairs))
        cur_read_id = read_ids[random.randint(0,cur_num_reads-1)]
        cur_read = reads_dict[cur_read_id]

        # writing to the output file
        out_bamfile.write(cur_read)

        # writing a single random pair from each unique position
        sys.stderr.write(
            "Finished! found {} duplicated segments in {} uniq segments ({:.2f}), max #duplictaes: {}"
            .format(
                cnt_dup_pos, len(uniq_reads_corr_dict),
                (1.0*cnt_dup_pos)/len(uniq_reads_corr_dict),
                max_num_of_dupliactes
            )
        )

        out_bamfile.close()
        in_bamfile.close()


def main():
    """Run as a script."""
    parser = argparse.ArgumentParser(
        "Identifies duplicated reads and select one (pair) randomly"
    )
    parser.add_argument("-p", action='store_true', dest='is_paired_ends',
                        help="Treat reads as paired")
    parser.add_argument("in_bam_file")
    parser.add_argument("out_bam_file")

    options = parser.parse_args()
    in_bam_file_name  = options.in_bam_file
    out_bam_file_name = options.out_bam_file

    # If randomize quality scores such that samtools will select a random pair
    # and not the highest score notice that this will make the quality scores
    # in downstream analysis random

    #  sys.stderr.write(
        #  "Warning: the -r flag is not implemented, the selection is "
        #  "random - i.e. quality independent\n"
    #  )

    if options.is_paired_ends:
        rmdup_pe(
            in_bam_file_name, out_bam_file_name
        )
    else:
        rmdup_se(
            in_bam_file_name, out_bam_file_name
        )

    sys.stderr.write("\nDone!\n")

if __name__ == '__main__':
    main()
