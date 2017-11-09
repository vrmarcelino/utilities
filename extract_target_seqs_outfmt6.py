#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract sequences with a good blast hit.
Takes as input the blastn output in tab format and a fasta file.
Python version 3.5

Created on Fri Sep  8 14:02:44 2017
@author: VanessaRM
"""

# Load modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser
import sys
import pandas as pd


# Help text
if len(sys.argv) == 1:
    print ()
    print ("Script to extract sequences from a fasta file with a good blast hit.")
    print ()
    print ("USAGE: extract_target_seqs.py -bi <blast_input_file.txt> -fi <fasta_input.fas> -o <output>")

    sys.exit()

# input files
parser = ArgumentParser()
parser.add_argument('-bi', '--blast_input', help='The path to the blast result', required=True)
parser.add_argument('-fi', '--fasta_input', help='Fasta file', required=True)
parser.add_argument('-o', '--output', help='output file', required=True)

args = parser.parse_args()
input_blast_res = args.blast_input
input_seqs = args.fasta_input
output_fp = args.output

#input_blast_res = "one_hit_all_fungi.txt"
#input_seqs = "all_fung_refseqs.fna"
#output_fp = "putative_rpb1.fas"


# Read blast output:
blast_tb = pd.read_table(input_blast_res, header = None)
wanted_contig_ids = pd.Series.tolist(blast_tb.iloc[:,0])

              
#### Extract contigs based on their seq name ####
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






