#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Put species names on reference dataset (greengenes and silva reference databases downloaded from Qiime)
Created on Sat Aug 29 08:58:12 2015
Modified in Oct 2017

@author: VanessaRM
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-i', '--input_fasta', help='The path to the reference fasta file', required=True)
parser.add_argument('-t', '--taxa_name', help='id_to_taxonomy_map.txt', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
input_fasta = args.input_fasta
taxa_name = args.taxa_name
output_fp = args.output


# comment later
#input_fasta = "97_otus.fasta"
#taxa_name = "97_otu_taxonomy.txt"
#output_fp = "97_otus_tax.fas"



new_aln = []
# Create a dictionary with sequence IDs:
dict_taxa = {}
tax = open(taxa_name, "r" )
for line in tax:
    identifier_split = line.split('\t')
    seq_id = str(identifier_split[0])
    
    lineage = identifier_split[1]
    lineage = lineage.replace(" ", "")
    
    dict_taxa ['%s' %seq_id] = [lineage]
    
    
for seq_record in SeqIO.parse(input_fasta, "fasta"):
    seq = str(seq_record.seq)
    seq_id_number = seq_record.id
    new_id = str(dict_taxa[seq_id_number])
    new_id = new_id.replace("['k__","")
    new_id = new_id.replace("n']","")
    new_id = new_id.replace("\\","")
    new_id = new_id.replace(";","")
    print (new_id)
    new_record = SeqRecord(Seq(seq), id= new_id, description='')
    new_aln.append(new_record)
    
    

count = SeqIO.write(new_aln,output_fp, "fasta")

print ("")
print ("Done. Saved %i sequences in the output file." % (count))
print ("")

