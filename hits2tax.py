#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uses the output of ISMA+A to filter the contigs for the taxon of interest.
Outputs a list of contigs IDs and their respective taxon (according to blast)
Works with Python 3.

USAGE: hits2tax.py -i trinity_snt.tax_uniqueAlignments.txt -c Bird.01.Trinity.Trinity.fasta -t Fungi -o fungi_out

Created on Thu Sep 15 10:42:59 2017
Last modified on Sep 25
@ V.R.Marcelino

"""
from ete3 import NCBITaxa
ncbi = NCBITaxa()
from Bio import SeqIO
from argparse import ArgumentParser
# if not installed, need to be done via PBS
#ncbi.update_taxonomy_database() # update the database regularly


# input files
parser = ArgumentParser()
parser.add_argument('-i', '--input_ismaa', help='The path to the isma-a result', required=True)
parser.add_argument('-c', '--input_contigs', help='Trinty contigs in fasta format', required=True)
parser.add_argument('-t', '--taxon_to_keep', help='The organism of interest', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
input_ismaa_res = args.input_ismaa
input_contigs = args.input_contigs
taxa2keep = args.taxon_to_keep
output_fp = args.output

# comment if not debugging
#input_ismaa_res = "trinity_snt.tax_uniqueAlignments.txt"
#input_contigs = "Bird.01.Trinity.Trinity.fasta"
#taxa2keep = "Fungi"
#output_fp = "fungi_out"


# Function that returns a taxonomy for each taxID
def get_lineage_from_taxonomy_db (taxid):
    lineage_ids = ncbi.get_lineage(taxid)
    lin = ncbi.translate_to_names(lineage_ids)
    print (taxid)
    return lin


# loop through the ismaa-a output and return the cotig IDs of the taxon of interest
QueryKeep_dic = {}
for line in open(input_ismaa_res):
    splits = line.strip().split('\t')
    taxonID = splits[2]
    
    # Exclude the "root" (-1)
    if taxonID != "-1":
        
        # run the function to get taxa data
        lineage = get_lineage_from_taxonomy_db(taxonID)
        
        # check if it is the target organism
        if taxa2keep in lineage:
            
            # only include seq if the contig ID is not in the list yet
            if splits[0] not in QueryKeep_dic.keys():
                QueryKeep_dic[ splits[0] ] = lineage


#save dict to file
contigIDs = str(output_fp + "_contigIDs.txt") # name of this output file
sf = open(contigIDs, "w")
for k, v in QueryKeep_dic.items():
    sf.write(str(k) + '\t' + str(v) + '\n')
sf.close()


### Get the contigs
filtered_contigs = []
print ("getting the contigs... ")
for seq_record in SeqIO.parse(input_contigs, "fasta"):
    if seq_record.id in QueryKeep_dic.keys():
        filtered_contigs.append (seq_record)


# Save contigs to file
filtered_contigs_out = str(output_fp + "_filtered_contigs.fas") # name of this output file
count_seqs = SeqIO.write(filtered_contigs, filtered_contigs_out, "fasta")
print ("")
print ("Saved wanted contig IDs in the %s file." %(contigIDs))
print ("Saved %i reads in the %s file." %(count_seqs, filtered_contigs_out))
print ("")

