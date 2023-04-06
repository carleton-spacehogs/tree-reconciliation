#!/bin/python3
from Bio import SeqIO
from parse_COG_calls import extract_specific_COG
import sys

COG_of_interest = sys.argv[1] # e.g. COG4766
COG_calling_method = sys.argv[2] # diamond or deepNOG
out_faa = sys.argv[3]
all_seq_fasta = sys.argv[4] # "SingleLine_EnrichedGenomes.faa"

record_dict = SeqIO.to_dict(SeqIO.parse(all_seq_fasta, "fasta"))
print(f"I finished parsing: {all_seq_fasta}!")

def give_COG_write_seq(record_dict, COG_of_interest):
	with open(out_faa, 'w') as out_f:
		for ORF in extract_specific_COG(COG_of_interest = COG_of_interest, call_from = COG_calling_method):
			out_f.write(f">{ORF}\n")
			out_f.write(str(record_dict[ORF].seq) + "\n")

# if you give me a lot of COGs, I can write this many many times
give_COG_write_seq(record_dict, COG_of_interest)
