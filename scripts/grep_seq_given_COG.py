#!/bin/python3
from Bio import SeqIO
from parse_COG_calls import extract_specific_COG
import sys
import os

COG = sys.argv[1] # e.g. COG4766

if "COG_calling_method" not in os.environ.keys():
    print("the environment is not set up yet, do these 2 commands:")
    print("    source ./scripts/declare_file_location.sh")
    print(f"    source ./scripts/declare_file_location.sh {COG}")
    exit(1)

COG_calling_method = os.environ['COG_calling_method']
out_faa = os.environ['gene_seq_file']
all_seq_fasta = os.environ['all_seq_fasta']

record_dict = SeqIO.to_dict(SeqIO.parse(all_seq_fasta, "fasta"))
print(f"I finished parsing: {all_seq_fasta}!")

def give_COG_write_seq(record_dict, COG):
	with open(out_faa, 'w') as out_f:
		for ORF in extract_specific_COG(COG_of_interest = COG, call_from = COG_calling_method):
			out_f.write(f">{ORF}\n")
			out_f.write(str(record_dict[ORF].seq) + "\n")

# if you give me a lot of COGs, I can write this many many times
give_COG_write_seq(record_dict, COG)
