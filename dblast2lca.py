#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dblast2lca - calculates the "LCA" and produces a taxonomy report.

This series of scripts calculates the lowest classification per read, produces a summary of the results, 
and sorts the contigs into taxonomic bins

@ V.R.Marcelino
Created on Tue Oct 10 16:35:28 2017
"""

from ete3 import NCBITaxa
ncbi = NCBITaxa()
import pandas as pd
import time
start_time = time.time()
from argparse import ArgumentParser


# input
parser = ArgumentParser()
parser.add_argument('-i', '--input_diamond_blast_result', help='The path to the diamond blast result', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
blast_in = args.input_diamond_blast_result
output_fp = args.output


# parameters:
blast_headers = ["contigID","sseqid","stitle","staxid","stuff1","stuff2","stuff3","stuff4"]


# Function that returns the lineage and the consensus taxa 
# (simialr to Lowest Common Ancestor) for a set of taxid
def get_lca (taxid_list):
    results = []
    if len(taxid_list) > 1:
            tree=ncbi.get_topology(taxid_list)
            try:
                common_ancestor=tree.get_common_ancestor(taxid_list)
                lineage_ids = ncbi.get_lineage(common_ancestor.name)
            except ValueError:
                lineage_ids = ("taxid %s bug" %taxid_list)
                pass

    else:
        try:
            lineage_ids = ncbi.get_lineage(taxid_list[0])
        except ValueError:
            lineage_ids = "taxid bug"
            pass

    if lineage_ids == "taxid bug": # deal with ValueError of "node names not found"
        results = lineage_ids
    else:
        try:
            lin = ncbi.translate_to_names(lineage_ids)
            lowest_class = lin[-1]
        except ValueError:
            lin = lineage_ids
            lowest_class = "no idea"
            pass
        results = [lowest_class,lin]
        
    return results


# Read diamond blast output:
blast_tb = pd.read_table(blast_in, header = None, names = blast_headers )


# Get a list of contig IDs:
all_ctgs = blast_tb.contigID.unique()


# loop though unique ids, blast output, and return the first report
report_dic = {}
for ctg in all_ctgs:
    
    match_lines = blast_tb.loc[blast_tb['contigID'] == ctg]
    
    # get the taxonomic classification
    taxids = match_lines.staxid.astype(str).values.tolist()
    taxclass = get_lca(taxids)
    print (taxclass[0])
    
    # get the UniProt ID and seq for the first record (best scored)
    sseqid = match_lines.sseqid.astype(str).values.tolist()[0]
    UniProtID = sseqid.split("|")[-1]
    stitle = match_lines.stitle.astype(str).values.tolist()[0]
    staxid = match_lines.staxid.astype(str).values.tolist()[0]
    
    report_dic[ctg] = [sseqid,UniProtID,stitle,staxid,taxclass[0],taxclass[1]]


# Save dict to file
tax_report = str(output_fp) # name of this output file
sf = open(tax_report, "w")
sf.write("contigID"+'\t'+"full_sseqid"+'\t'+"UniProtID"+'\t'+"stitle"+'\t'+"staxid"+'\t'+"LCA"+'\t'+"taxrank1"+'\t'+"taxrank2"+'\t'+"taxrank3"+'\t'+"taxrank4"+'\t'+"taxrank5"+'\t'+"taxrank6"+'\t'+"taxrank7"+'\t'+"taxrank8"+'\t'+"taxrank9"+'\t'+"taxrank10"+'\n')
for k, v in report_dic.items():
    
    clean_values = str(v)
    clean_values = clean_values.replace("['","")
    clean_values = clean_values.replace("']]","")
    clean_values = clean_values.replace(";","_")
    clean_values = clean_values.replace("', '","\t")
    clean_values = clean_values.replace("', root","")
    
    sf.write(str(k) + '\t' + clean_values + '\n')
sf.close()


# Run time
time_min = round(((time.time() - start_time) * 0.0166667),2)
print(" Done in %s minutes" % (time_min))

