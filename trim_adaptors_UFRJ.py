#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Trim tail and adaptors from seq sequences.
The 16S is a short amplicon, so the adaptor and tail has been sequenced in the final end of these sequences. 
We need therefore to trim them in order to merge the paired end reads.
As the adpators used by Lais are different from the ones I used, I'll use the primer sequences to trim them


Created on Jul 30 2014
modified: April 2019
"""
from Bio import SeqIO
import regex
import sys

# Help
if len(sys.argv) == 1:
    print ""
    print "Script to trimm the tail and adaptors from short freads"
    print ""
    print "Usage: supply the file name and inform if its forward or reverse"
    print "ex: python trim_adaptors.py my_fastq_file.fastq forward"
    print ""
    print "If you want to run it for multiple files, use the shell:"
    print "for file in *_R1_001.fastq; do trim_adaptors.py $file forward; done >> screen.out 2>> screen.err &"
    print ""
    sys.exit()
    
    
# input files and arguments:
input_file = str(sys.argv[1])
output_file = input_file + "_trimmed.fastq"


# Tail in forward reads:
# remove reverse primer (reverse complement of 16S R806)
# full primer: GGACTACHVGGGTWTCTAAT
# Reverse complement (the one that will match to the forward read): ATTAGAWACCCBDGTAGTCC
primer_tail_f = '(ATTAGA\wACCC\w\wGTAGTCC){e<=1}' # allowing 1 mismatch (\w = ambiguous positions)


# Tail in reverse reads:
# full primer 16S F515: GTGCCAGCMGCCGCGGTAA
# reverse complement o the forward primer: TTACCGCGGCKGCTGGCAC
primer_tail_r = '(TTACCGCGGC\wGCTGGCAC){e<=1}' # allowing 1 mismatch (\w = ambiguous positions)

#define trimmer function
def trimmer (records, tail):
    "Trims the rc primer and all the bases after it (adaptors and etc...)"
    for record in records:
        sequence = str(record.seq)        
#        index = record.seq.find(tail)
        index = regex.search((tail), sequence)       
        if index == None:
            #tail not found
            yield record
        else:
            # Trim it!
            cut_off = (int(index.span()[0]))
            # Check if it is in the end of the sequence (>100bp)
            if cut_off > 100:
                yield record [:cut_off]
            else:
                yield record
        

original_reads = SeqIO.parse(input_file, "fastq")

if str(sys.argv[2]) == "forward":
    trimmed_reads = trimmer(original_reads, primer_tail_f)
    
elif str(sys.argv[2]) == "reverse":
    trimmed_reads = trimmer(original_reads, primer_tail_r)

#Else return error:
else:
    print ("")
    print ("Please define if the file is forward or reverse!!")
   

count = SeqIO.write(trimmed_reads, output_file, "fastq")
print""
print "Saved %i reads in the %s file." % (count, output_file)


