#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Clean reference dataset in fasta format
Keeps only one seq per genus.

Created on Thu Sep 14 09:18:56 2017
@author: VanessaRM
# -*- coding: utf-8 -*-
"""

from Bio import SeqIO
import sys

#help
if len(sys.argv) == 1:
    print ("")
    print ("Script to refine the reference datasets.")
    print ("")
    print ("Usage: del_repetitive_species_fas.py dataset.fas")
    print ("")
    print ("")
    sys.exit()

all_species = []
unique_genera = []
unique_records = []


input_dataset = str(sys.argv[1])
#input_dataset = "ITSSeq_renamed_no_dup.fas"


print ("")
print ("cleaning dataset...")

# All species
for seq_record in SeqIO.parse(input_dataset, "fasta"):
    all_species.append(seq_record.id)

# All seqs
#for seq_record in SeqIO.parse(input_dataset, "fasta"):
#    all_seqs.append(seq_record.seq)

len(unique_genera)
# Find the unique genera
for seq_record in SeqIO.parse(input_dataset, "fasta"):
    seq_name = seq_record.id
    genus = seq_name.split("_")[1]
    
    if genus not in unique_genera:
        unique_genera.append(genus)
        unique_records.append(seq_record)

#print results
deleted_records = (len(all_species) - len(unique_records))
print ("")

print ("number of deleted records is %i" %(deleted_records))
print ("%i species remained in the dataset" %(len(unique_records)))

# write to file
SeqIO.write(unique_records, "ITSSeq_one_per_genus.fas", "fasta")
print ("")
print ("Done!")

