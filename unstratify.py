#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to convert the output from humann2 
(stratified merged_paths.tsv or merged_genes.tsv)
into unstratified files containing the pathways only or the pathways per species only (without their sum).

@ V.R.Marcelino
Created on Tue Nov 28 15:40:43 2017
"""

import csv

input_paths = "merged_path_abundance.tsv"
input_genes = "merged_gene_families.tsv"

only_paths = []
only_species = []
with open(input_paths) as tsv:
    for line in csv.reader(tsv):
        if "|" in str(line):
            only_species.append(line)
        else:
            only_paths.append(line)
            

with open('paths.tsv', 'w') as tsvfile:
    tsvfile.writelines('\t'.join(i) + '\n' for i in only_paths)
with open('sp_paths.tsv', 'w') as tsvfile:
    tsvfile.writelines('\t'.join(i) + '\n' for i in only_species)



only_genes = []
only_species_genes = []
with open(input_genes) as tsv:
    for line in csv.reader(tsv):
        if "|" in str(line):
            only_species_genes.append(line)
        else:
            only_genes.append(line)
            

with open('genes.tsv', 'w') as tsvfile:
    tsvfile.writelines('\t'.join(i) + '\n' for i in only_genes)
with open('sp_genes.tsv', 'w') as tsvfile:
    tsvfile.writelines('\t'.join(i) + '\n' for i in only_species_genes)
