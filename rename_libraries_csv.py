#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rename columns of a csv file based on a csv file with old and new names.

@ V.R.Marcelino
Created on 1 Aug 2019
"""

import pandas as pd

file2rename = "0_EukFamily_table_Hunters_no_contam_ed.csv"
file_w_names = "SRA2SampleID.csv"

# parse old names, adding the '..res', and add to a dictionary.
names_df = pd.read_csv(file_w_names, sep=',', index_col="Library_Name")
names_dic = {}

for index, row in names_df.iterrows():
    index = index + ".res"
    names_dic[index] = row[0]
    

# now read the file o rename and replace teh names:
df = pd.read_csv(file2rename, sep=',')

new_df = df.rename(columns=names_dic)

pd.DataFrame.to_csv(new_df, "1_EukFamily_table_Hunters_no_contam_ed_renmaed.csv", index=False)

