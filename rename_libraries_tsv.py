#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rename columns of a OTU table (tsv file) based on a csv file with old and new names.

@ V.R.Marcelino
Created on 11 Sep 2019
"""

import pandas as pd

file2rename = "CCM_Species_table_Hunters.csv"
file_w_names = "FrCode2VRMs.csv"

# parse old names, adding the '..res', and add to a dictionary.
names_df = pd.read_csv(file_w_names, sep=',', index_col="FranCode")
names_dic = {}

for index, row in names_df.iterrows():
    names_dic[index] = row[0]
    

# now read the file o rename and replace teh names:
df = pd.read_csv("16S_OTU_table.txt", sep='\t')

new_df = df.rename(columns=names_dic)

pd.DataFrame.to_csv(new_df, "renamed_table.csv", index=False)

