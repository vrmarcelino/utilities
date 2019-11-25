#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Produce links file for Circos fig
@ V.R.Marcelino
Created on Wed Jan 10 10:42:32 2018
"""

import pandas as pd

# read file with gene names and positions
bands = pd.read_csv("bands.csv")

# Files containing the co-occurring genes
genes_b1 = pd.read_csv("Bird1_genes.csv")
genes_b2 = pd.read_csv("Bird2_genes.csv")
genes_b3 = pd.read_csv("Bird3_genes.csv")
genes_b4 = pd.read_csv("Bird4_genes.csv")
genes_b5 = pd.read_csv("Bird5_genes.csv")
genes_b6 = pd.read_csv("Bird6_genes.csv")
genes_b7 = pd.read_csv("Bird7_genes.csv")
genes_b8 = pd.read_csv("Bird8_genes.csv")
genes_b9 = pd.read_csv("Bird9_genes.csv")
genes_b10 = pd.read_csv("Bird10_genes.csv")
genes_b11 = pd.read_csv("Bird11_genes.csv")


# function to add gene positions in a links pd.df table:
    
def add_to_table (genes_lib, colour):
    global links
    nrows = len(genes_lib.index)
  
    for i in range(nrows): 
        for j in range(i+1, nrows):
            gene1 = genes_lib.iloc[i,1]
            gene2 = genes_lib.iloc[j,1]

            # find info in the bands file
            target_line_g1 = bands.loc[bands['Gene'] == gene1]
            antib_g1 = target_line_g1.iloc[0,0]
            pos1_g1 = target_line_g1.iloc[0,2]
            pos2_g1 = target_line_g1.iloc[0,3]
        
            target_line_g2 = bands.loc[bands['Gene'] == gene2]
            antib_g2 = target_line_g2.iloc[0,0]
            pos1_g2 = target_line_g2.iloc[0,2]
            pos2_g2 = target_line_g2.iloc[0,3]
            c="color="
            new_line = [[antib_g1,pos1_g1,pos2_g1,antib_g2,pos1_g2,pos2_g2,c,colour]]
        
            links = links.append(new_line, ignore_index=True)


links = pd.DataFrame()

# Birds - colored by WTP (red) or not (grey):
add_to_table (genes_b1, 'rdylbu-11-div-1')
add_to_table (genes_b2, 'rdylbu-11-div-2')
add_to_table (genes_b3, 'greys-7-seq-2')
add_to_table (genes_b4, 'greys-7-seq-3')
add_to_table (genes_b5, 'rdylbu-11-div-3')
add_to_table (genes_b6, 'rdylbu-11-div-4')
add_to_table (genes_b7, 'greys-8-seq-2')
add_to_table (genes_b8, 'greys-7-seq-4')
add_to_table (genes_b9, 'greys-7-seq-5')
add_to_table (genes_b10, 'greys-7-seq-6')


# note that bird11 has to be added manually (only two genes)!
# Bla 11 12 Tet 4 5 color=greys-7-seq-7

#Save
links.to_csv(r'links.txt', header=None, index=None, sep=' ', mode='a')





