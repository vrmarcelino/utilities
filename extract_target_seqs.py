#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract sequences with a good blast hit.
Takes as input the blastn output and a Trinity contig in fasta format.
Python version 3.5

Created on Fri Sep  8 14:02:44 2017
@author: VanessaRM
"""

# Load modules
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser
import sys

# Help text
if len(sys.argv) == 1:
    print ()
    print ("Script to extract sequences from a Trinity contig file that has a blast hit")
    print ()
    print ("USAGE: extract_target_seqs.py -bi <blast_input_file.txt> -ci <contigs_input.fas>")

    sys.exit()

# input files
parser = ArgumentParser()
parser.add_argument('-bi', '--blast_input', help='The path to the blast result', required=True)
parser.add_argument('-ci', '--contigs_input', help='Trinty contigs in fasta format', required=True)

args = parser.parse_args()
input_blast_res = args.blast_input
input_seqs = args.contigs_input
output_fp = "ITS_contigs.fas"
#input_blast_res = "bird1_fragment.txt"
#input_seqs = "Bird.01.Trinity.Trinity.fasta"



### OBTAIN THE ID OF CONTIGS MATCHING WITH ITS #### 

# Function that returns whether a given query (contig) had an ITS match or not
pattern_match = re.compile("Sequences producing significant alignments")
pattern_no_match = re.compile("No hits found")
def get_matches (n, all_the_lines):
    # n is the line to start search
    for i in all_the_lines[n:]:
        if pattern_match.search(i) != None:
            return "ITS"
        if pattern_no_match.search(i) != None:
            return "Something else"

# species name:
pattern_name = re.compile("Query=")

# loop through queries and find the ones that match with fungal ITS
wanted_contig_ids = []
with open(input_blast_res) as file:
  all_lines = file.readlines()
  count_lines = 0
  
  for line in all_lines:
      count_lines +=1
      
      query_line = pattern_name.search(line)
      if query_line != None: # Found a name, now check if it matches ITS:
          check = get_matches (count_lines,all_lines)
          
          if check == "ITS":
              splitted_line = line.split()
              wanted_contig_ids.append(splitted_line[1])
              
              
#### EXTRACT THE CONTIGS BASED ON THE IDs AND CLEAN SEQ NAME ####

extracted_contigs = []
for seq_record in SeqIO.parse(input_seqs, "fasta"):
    if seq_record.id in wanted_contig_ids:
        
        full_id = seq_record.id
        new_id = full_id.split()[0]
        seq = str(seq_record.seq)
        new_record = SeqRecord(Seq(seq), id= new_id, description='')
        
        extracted_contigs.append(new_record)

count_seqs = SeqIO.write(extracted_contigs, output_fp, "fasta")
print ()
print ("Saved %i reads in the %s file." %(count_seqs, output_fp))
print ()
