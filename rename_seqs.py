#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename fasta files based on a csv file with old and new names
@ V.R.Marcelino
Created on Mon Mar 12 15:36:17 2018
"""

import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = ArgumentParser()
parser.add_argument('-i', '--old_new_file', help='The path to the file containing old and new names', required=True)
parser.add_argument('-f', '--fasta', help='The path to the fasta file', required=True)
parser.add_argument('-o', '--output', help='The renamed file name', required=True)


args = parser.parse_args()
rename_file = args.old_new_file
fasta_fp = args.fasta
output_fp = args.output

#rename_file = "Rename_SRA_Runs.csv"
#fasta_fp = "Mito.mapped.fasta"
#output_fp = "Mito_mapped_renamed.fasta"

# store old and new names in a dict:
renames_df = pd.read_csv(rename_file,header=None,index_col=0)
names_dict = renames_df.to_dict('dict')[1]


# new fasta file with renamed seqs:
new_fasta = []
for seq_record in SeqIO.parse(fasta_fp, "fasta"):
    old_id = seq_record.id
    if old_id in names_dict:
        new_id = names_dict[old_id]        
        seq = str(seq_record.seq)
        new_record = SeqRecord(Seq(seq), id=new_id, description='')
        new_fasta.append(new_record)
    else:
        new_fasta.append(seq_record)


# Save renamed seqs
count = SeqIO.write(new_fasta,output_fp, "fasta")

print ("")
print ("Done. Saved %i sequences in the output file." % (count))
print ("")


