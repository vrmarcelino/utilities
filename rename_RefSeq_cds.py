#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to rename fasta files that only ocntain an UID - put a species name.
Also removes identical species (gene duplications) and keeps only the longer seq
@ V.R.Marcelino
Created on Wed Nov  8 14:06:43 2017
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Entrez
import re

Entrez.email = "vrmarcelino@gmail.com"
input_seqs = "putative_rpb1.fas"
output_fp = "putative_rpb1s_renamed.fas"

# Function that returns a species name from a accession ID:
def taxonomer (UID):  
    handle = Entrez.efetch("nucleotide",rettype="gb", id=UID)
    results = SeqIO.read(handle, 'genbank')
    handle.close()
    desc_split = results.description.split(" ")
    species = desc_split[0] + "_" + desc_split[1]
    return species


# storage for the renamed seqs
new_fasta = []

# loop through sequences and rename them
# Also check for duplicates and chose the longets seq
new_fasta = {}
for seq_record in SeqIO.parse(input_seqs, "fasta"):
    get_gi = re.split("_|\|", seq_record.id)
    accession = get_gi[1] + "_" + get_gi[2]
    species_name = taxonomer(accession)
    seq = str(seq_record.seq)
    new_record = SeqRecord(Seq(seq), id=species_name, description='')
    
    if species_name in new_fasta.keys():
        # compare lenghts
        l_seq1 = len(new_fasta[species_name])
        l_seq2 = len(seq)
        print (species_name)
        print (l_seq1)
        print (l_seq2)
        if l_seq2 > l_seq1:
            new_fasta[species_name] = seq
    else:
        new_fasta[species_name] = seq
        print ("new")
        
# Save renamed seqs
count = 0
ofile = open(output_fp, "w")
for i in range(len(new_fasta)):
    count +=1
    ofile.write(">" + list(new_fasta.keys())[i] + "\n" + list(new_fasta.values())[i] + "\n")
ofile.close()


print ("")
print ("Done. Saved %i sequences in the output file." % (count))
print ("")

