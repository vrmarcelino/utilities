#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to get accession IDs from the RefSeq CDSs.
rename fasta files that only contain an UID - put a species name.
It does not check for repetitive species like the version 1 of this script.
Used for the complete RefSeq fungi (not only RPB1)
USAGE: rename_RefSeqs_v2.py -i <ref.fna> -e vrnarcelino@gmail.com -o <ref_renamed.fna>
@ V.R.Marcelino
Created on Wed Nov 15 14:06:43 2017
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Entrez
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--input_seqs', help='The path to the reference fasta file', required=True)
parser.add_argument('-e', '--email', help='your email, required for Entrez', required=True)
parser.add_argument('-o', '--output', help='output file path', required=True)

args = parser.parse_args()
input_seqs = args.input_seqs
Entrez.email = args.email
output_fp = args.output


input_seqs = "all_fung_refseqs.fna"
Entrez.email = "vrmarcelino@gmail.com"
output_fp = "accession_IDs.fas"


## loop through sequences and get accession numbers
list_accessions = []
for seq_record in SeqIO.parse(input_seqs, "fasta"):
    get_gi = re.split("_|\|", seq_record.id)
    accession = get_gi[1] + "_" + get_gi[2]
    list_accessions.append(accession)
    

unique_accessions = list(set(list_accessions))


# save list
with open("unique_accession_IDs.txt", 'w') as file_handler:
    for item in unique_accessions:
        file_handler.write("{}\n".format(item))


