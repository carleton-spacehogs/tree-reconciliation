#!/bin/python3
import sys
from Bio import SeqIO

'''Given an alignment file, only keep sequences that have 
	above X percent of align sequence (instead of gaps)'''
def trim(in_file, out_file, percent = 80):
	if percent > 100 or percent < 0:
		print("The percentage must be between 0 to 100 (inclusive)")
	with open(out_file, 'w') as out_f:
		for seq_obj in SeqIO.parse(in_file, "fasta"):
			sequence = seq_obj.seq
			if sequence.count('-')/len(sequence)*100 <= 20:
				out_f.write(f">{seq_obj.name}\n")
				out_f.write(str(sequence) + "\n")
	print(f"Finish trimming, outputted to {out_file}")

def main():
	trim(in_file = sys.argv[1], out_file = sys.argv[2], percent = float(sys.argv[3]))

if __name__ == "__main__":
	main()
