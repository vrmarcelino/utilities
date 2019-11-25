#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename simulated metagenomics datasets

@ V.R.Marcelino
Created on Tue Mar 27 15:33:54 2018
Modified on 02 - Jul - 2018

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fasta_fp_R1 = "sim_metatrans_R1_long_names.fas"
fasta_fp_R2 = "sim_metatrans_R2_long_names.fas"


def renamer(records, PE):
    for seq_record in records:
        
        old_id = seq_record.id
        new_id = old_id.replace( "/" , "_") + PE

        seq = str(seq_record.seq)

        new_rec = SeqRecord(Seq(seq), id=new_id, description="")
        yield new_rec
        

# Do it and save R1
fasta_parser_R1 = SeqIO.parse(fasta_fp_R1, "fasta")
renamed_fasta_R1 = "sim_metatrans_R1_renam.fas"
count = SeqIO.write(renamer(fasta_parser_R1,"/1"),renamed_fasta_R1, "fasta")

print ("")
print ("Saved %i sequences in the output file." % (count))
print ("")


# Do it and save R1
fasta_parser_R2 = SeqIO.parse(fasta_fp_R2, "fasta")
renamed_fasta_R2 = "sim_metatrans_R2_renam.fas"
count = SeqIO.write(renamer(fasta_parser_R2,"/2"),renamed_fasta_R2, "fasta")

print ("")
print ("Saved %i sequences in the output file." % (count))
print ("")


