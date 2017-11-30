#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to count and normalize reads in a SAM file.

* Exclude PE reads that match to different species

* Normalizses as fragments (i.e. PE are considered 1 hit) per million reads

* Counts are a generic analog to the CPM in Humann2, or TPM (Transcipts per Million), except that gene length is not taken into consideration.

* Gene length are not taken into consideration because we are only using one gene

* WGD not taken into consideration because I didn't find evidence for a RPB1 paralog

@ V.R.Marcelino
Created on Wed Nov 29 16:24:18 2017

"""
# import stuff
import HTSeq
import pandas as pd

reads_fp = "reads.txt"
input_species = "Fungal_species.txt"
input_folder = "SAM"

# Open files
sam_files = !ls {input_folder}/*.sam
species = pd.read_table(input_species, index_col="Species")
reads = pd.read_table(reads_fp, index_col = "Sample")



## Function that takes a sam file and return a dict with reads and species
def get_quality_counts(sam):
    quality_hits = {}
    for alnmt in sam:
        if alnmt.read.name in quality_hits.keys():
            if quality_hits[alnmt.read.name] == alnmt.iv.chrom:
                print ("PE matches, hit already counted")
            else:
                print ("ambiguous, deleting hit")
                quality_hits.pop(alnmt.read.name, None)
        else:
            print("new match, adding hit")
            quality_hits[alnmt.read.name] = alnmt.iv.chrom
            
    # invert dict to have species as keys and reads as values:
    inv_quality_hits = {}
    for k, v in quality_hits.items():
        inv_quality_hits[v] = inv_quality_hits.get(v, [])
        inv_quality_hits[v].append(k)

    return inv_quality_hits


#Loop trhough the SAM files and count:
for sam in sam_files:
    one_sam = HTSeq.SAM_Reader(sam)
    sample_name = sam.split("/")[1]
    species[sample_name] = "NA" # create a new column
    n_reads = reads.loc[sample_name] [0] # get the number of reads in this sample

    counts_dict = get_quality_counts(one_sam) # get a dict with counts

    # loop trhough the species and get their counts:
    for ind, row in species.iterrows():
        n_raw = (len(counts_dict.get(ind,())))
    
        # Calculate the number of fragments per million reads (CHECK!!)
        n = n_raw * 1000000 / n_reads

        # append to the data_frame
        species.loc[ind,sample_name] = n



pd.DataFrame.to_csv(species, "reads_count_norm.csv")



    
