#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Replace gene ID with a taxID
@ V.R.Marcelino
Created on Thu Jun 13 17:05:30 2019
"""
import pandas as pd
import re
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-if', '--input_frags', help='The path to the good frags file', required=True)
parser.add_argument('-ir', '--input_rsem', help='path to the RSEM file', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
frag_good = args.input_frags
RSEM_res = args.input_rsem
output_fp = args.output + ".tsv"

#frag_good = "ALG_1pass.txt"
#RSEM_res = "ALG_1.RSEM.isoforms.results.tsv"


# Make a dictionary of transcritp to taxID:

trans_dict = {}
frags = open(frag_good,'r')
for line in frags:
    splited = re.split("\t|\|",line)
    taxid = splited[0]
    transcript_full = splited[-1]
    transcript = re.sub(r'_i\d+_len.*','',transcript_full)
    transcript = re.sub(r'\n','',transcript) # remove teh new line char
    trans_dict[transcript] = taxid




# Now loop through rge RSEM results and replace gene IDS by taxids:
RSEM_df =  pd.read_csv(RSEM_res, index_col=0, sep='\t')
for index, row in RSEM_df.iterrows():
    gene_id = row['gene_id']
    if gene_id in trans_dict.keys():
        RSEM_df.at[index,'gene_id'] = trans_dict[gene_id]
    else:
        RSEM_df.at[index,'gene_id'] = 'unk_taxid'


pd.DataFrame.to_csv(RSEM_df, output_fp, sep='\t', header=True)
    

