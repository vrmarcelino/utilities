#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate TPM per taxa
Takes as input a file with taxids and TPM counts
@ V.R.Marcelino
Created on 14 June 2019
"""

import pandas as pd 
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import cTaxInfo # where we define classes used here

from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-i', '--input_file', help='The path to the file containing taxids and transcripts TPM', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
rsem_res = args.input_file
output_fp = args.output + ".csv"

#rsem_res="ALG1_taxidTPM.tsv"


# function to get taxa from taxid:
def lineage_extractor(taxid,TaxInfo_object):
    list_of_taxa_ranks = ['superkingdom', 'kingdom', 'phylum', 'subphylum','class', 'order', 'family','genus', 'species']
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    names = ncbi.get_taxid_translator(lineage)
    
    for key, val in ranks.items():
        
        if val == list_of_taxa_ranks[0]:
            TaxInfo_object.Superkingdom = names[key]
            TaxInfo_object.Superkingdom_TaxId = key            
  
        elif val == list_of_taxa_ranks[1]:
            TaxInfo_object.Kingdom = names[key]
            TaxInfo_object.Kingdom_TaxId = key

        elif val == list_of_taxa_ranks[2]:
            TaxInfo_object.Phylum = names[key]
            TaxInfo_object.Phylum_TaxId = key
            
        elif val == list_of_taxa_ranks[3]:
            TaxInfo_object.SubPhylum = names[key]
            TaxInfo_object.SubPhylum_TaxId = key
            
        elif val == list_of_taxa_ranks[4]:
            TaxInfo_object.Class = names[key]
            TaxInfo_object.Class_TaxId = key
    
        elif val == list_of_taxa_ranks[5]:
            TaxInfo_object.Order = names[key]
            TaxInfo_object.Order_TaxId = key
          
        elif val == list_of_taxa_ranks[6]:
            TaxInfo_object.Family = names[key]
            TaxInfo_object.Family_TaxId = key
   
        elif val == list_of_taxa_ranks[7]:
            TaxInfo_object.Genus = names[key]
            TaxInfo_object.Genus_TaxId = key

        elif val == list_of_taxa_ranks[8]:
            TaxInfo_object.Species = names[key]
            TaxInfo_object.Species_TaxId = key

    return TaxInfo_object


# read rsem output
rsem_df = pd.read_csv(rsem_res, sep='\t', index_col=0)

new_df = pd.DataFrame(columns=['TPM','Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family','Genus', 'Species'])



# add a column with species name
for index, row in rsem_df.iterrows():
    
    TPM_number = row[0]
    new_df.at[index, 'TPM'] = TPM_number

    taxid = index
    
    if taxid == 'unk_taxid':
        pass
    
    else:
        match_info = cTaxInfo.TaxInfo()
        match_info = lineage_extractor(int(taxid), match_info)
        new_df.at[index, 'Superkingdom'] = match_info.Superkingdom
        new_df.at[index, 'Kingdom'] = match_info.Kingdom
        new_df.at[index, 'Phylum'] = match_info.Phylum
        new_df.at[index, 'Class'] = match_info.Class
        new_df.at[index, 'Order'] = match_info.Order
        new_df.at[index, 'Family'] = match_info.Family
        new_df.at[index, 'Genus'] = match_info.Genus
        new_df.at[index, 'Species'] = match_info.Species 

### Save
pd.DataFrame.to_csv(new_df, output_fp)

print ("Done!")



