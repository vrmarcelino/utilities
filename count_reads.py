#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Count number of reads aligned to fungi in a bed file.
Normalizes the counts sequencing depth (i.e. gets relative abundances)

@ V.R.Marcelino
Created on Mon Nov 13 12:21:42 2017
"""

import pandas as pd

input_folder = "05_fungal_mRNA"
input_species = "Fungal_species.txt"
reads_fp = "n_reads.txt"

bed_files = !ls {input_folder}/*bed

species = pd.read_table(input_species, index_col="Species")
reads = pd.read_table(reads_fp, index_col = "Sample")

for f in bed_files:
    one_bed = open(f, 'r').read()
       
    sample_name = f.split("/")[1]
    species[sample_name] = "NA" # create a new column
    n_reads = reads.loc[sample_name] [0] # get the number of reads in this sample
    
    for ind, row in species.iterrows():
        n_raw = one_bed.count(ind)
        
        # Calculate the relative abundance of reads (%)
        n = n_raw * 100 / n_reads
        
        species.loc[ind,sample_name] = n


pd.DataFrame.to_csv(species, "reads_count_norm.csv")

