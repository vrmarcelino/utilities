#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to rename fasta files that only contain an UID - put a species name.
It does not check for repetitive species like the version 1 of this script.
Used for the complete RefSeq fungio (not only RPB1)
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


#input_seqs = "GCF_000001985.1_JCVI-PMFA1-2.0_cds_from_genomic.fna"
#Entrez.email = "vrmarcelino@gmail.com"
#output_fp = "all_fung_refseqs_renamed.fas"

# Function that returns a species name from a accession ID:
def taxonomer (UID):
    handle = Entrez.efetch("nucleotide",rettype="gb", id=UID)
    results = SeqIO.read(handle, 'genbank')
    handle.close()
    desc_split = results.description.split(" ")
    species = desc_split[0] + "_" + desc_split[1]
    return species


# storage for the renamed seqs
new_fasta = []

# loop through sequences and rename them
previous_gi = ""
previous_seq_name = ""
for seq_record in SeqIO.parse(input_seqs, "fasta"):
    get_gi = re.split("_|\|", seq_record.id)
    accession = get_gi[1] + "_" + get_gi[2]

    #only access NCBI if it is a different accession ID
    if accession != previous_gi:
        species_name = taxonomer(accession)
        previous_seq_name = species_name
        previous_gi = accession

    else:
        species_name = previous_seq_name
    
    new_seq_name = species_name + "_" + seq_record.id
    seq = str(seq_record.seq)
    new_record = SeqRecord(Seq(seq), id=new_seq_name, description='')
    new_fasta.append(new_record)
    print (new_seq_name)
        
# Save renamed seqs
count = SeqIO.write(new_fasta,output_fp, "fasta")

print ("")
print ("Done. Saved %i sequences in the output file." % (count))
print ("")

