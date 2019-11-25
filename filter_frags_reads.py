#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script saves the classified transcripts in a file "out_frags_good" and the unclassified transcripts in "out_frags_fail"
# the files that did not classify with CCMetagen will get a "unk_taxid|unclassified" name
Difference - with filter_frags_transcripts - it get the different reads (rather than teh basename of the contig)
@ V.R.Marcelino
Created on Thu Mar  7 17:12:42 2019
"""
# load modules
import csv
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-if', '--input_frags', help='The path to the file containing frags', required=True)
parser.add_argument('-iccm', '--input_ccmetagen', help='The path to the ccmetagen result csv file', required=True)
parser.add_argument('-o', '--output', help='output prefix', required=True)

args = parser.parse_args()
frags_fp = args.input_frags
res_fp = args.input_ccmetagen
out_frags_good = args.output + "pass.txt"
out_frags_fail = args.output + "fail.txt"


# Read and store species in res file
res_match = []
with open(res_fp) as res:
    next (res) # skip first line
    for line in csv.reader(res):
        match = line[0]
#        species = line[18]
        res_match.append(match)


# Read frags file and save the filtered one
ofile_good = open(out_frags_good, "w")
ofile_fail = open(out_frags_fail, "w")

cg = 0 #count good reads
cb = 0 #count failed reads
with open(frags_fp) as frags:
    for line in csv.reader(frags, delimiter = '\t'):
        match = line[5]
        transcript_id = line[6].split(" ")[0]
        transcript_read = line[6].split("/")[-1]
        transcript = transcript_id + "/" + transcript_read

        if match in res_match:
            ofile_good.write(match + "\t" + transcript + "\n")
            cg += 1
            
        else:
            ofile_fail.write("unk_taxid|unclassified" + "\t" + transcript + "\n")
            cb += 1
            
ofile_good.close()
ofile_fail.close()

print ("")
print ("Saved %i sequences in the frags.pass.txt file." %(cg))
print ("%i reads did not meet the coverage QC threshold" %(cb))
print ("")


