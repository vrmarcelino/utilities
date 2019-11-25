#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rename simulated metagenomics datasets

@ V.R.Marcelino
Created on Tue Mar 27 15:33:54 2018
Modified on 22 - Jun - 2018

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


fastq_fp = "sim_metagenome_interleaved_raw.fastq"

fastq_parser = SeqIO.parse(fastq_fp, "fastq")


def renamer(records):
    for seq_record in records:
        old_id = seq_record.id
        ncbi_full = seq_record.description.split()[1]
        ncbi_ID = ncbi_full.split("=")[1]
    
        new_ID = ncbi_ID + "_" + old_id
        seq = str(seq_record.seq)
        description = seq_record.description
        quals = seq_record.letter_annotations

        new_rec = SeqRecord(Seq(seq), id=new_ID, description=description, letter_annotations=quals)
        yield new_rec
        

# Save
renamed_interleaved = "sim_metagenome_interleaved.fastq"
count = SeqIO.write(renamer(fastq_parser),renamed_interleaved, "fastq")

print ("")
print ("Saved %i sequences in the output file." % (count))
print ("")



####### split the fastq - for some reason, this is not working!! to be fixed. using deinterleave_fastq.sh instead
#paired_file = "sim_metagenome_interleaved.fastq"

#R1_file = "sim_metagen_R1.fastq"
#R2_file = "sim_metagen_R2.fastq"

#def splitter_R1(records):
#    for seq_record in records:
#        if "/1" in seq_record.id:
#            yield seq_record
            
 
#def splitter_R2(records):
#    for seq_record in records:
#        if "/2" in seq_record.id:
#            yield seq_record
           
#fastq_parser2 = SeqIO.parse(paired_file, "fastq")
#SeqIO.write(splitter_R1(fastq_parser2),R1_file, "fastq")

#fastq_parser2 = SeqIO.parse(paired_file, "fastq")
#SeqIO.write(splitter_R2(fastq_parser2),R2_file, "fastq")


