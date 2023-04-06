#!/bin/python3
import sys
import os
import csv

deepnog_COG_call_file = "SingleLine_deepnog.csv"
confidence_level = 0.7

diamond_blast_output = "diamond-out.tsv"
COG_ref = "/researchdrive/zhongj2/cog-20.cog.csv"
requirements = "minimum 30% identity, maximum e-value of 10-5, minimum 70% subject and query alignment"

'''
return: a list of ORF names of a specific COG
'''
def extract_specific_COG(COG_of_interest, call_from = "diamond"):
	f = []
	if call_from == "diamond":
		f = open_diamond_COG_calls()
	if call_from == "deepNOG":
		f = open_deepNOG_COG_calls()
	ORFs_of_interested_COG = []
	for columns in f:
		ORF_name, COG, accuracy = columns # accuracy can be eval from diamond
		if COG == COG_of_interest:
			ORFs_of_interested_COG.append(ORF_name)
	return ORFs_of_interested_COG

'''
count the number of ORFs in each COG
'''
def count_COGs_calls(call_from = "diamond"):
	out_count_f = call_from + "_COG_count.txt"
	
	COG_calls = []
	COG_count_dict = {}
	if call_from == "deepNOG":
		COG_calls = open_deepNOG_COG_calls()
	elif call_from == "diamond":
		COG_calls = open_diamond_COG_calls()
	
	for r in COG_calls:
		COG = r[1]
		if COG in COG_count_dict:
			COG_count_dict[COG] += 1
		else:
			COG_count_dict[COG] = 1

	# sorted_COG_count_dict is a list of tuples
	sorted_COG_count_dict = sorted(COG_count_dict.items(), key=lambda x:x[1])
	print(f"go to {out_count_f} to see the ORF counts")
	with open(out_count_f, "w") as f:
		for tuple in sorted_COG_count_dict:
			f.write(f"{tuple[0]},{tuple[1]}\n")


# deepNOG specific stuff begins
def open_deepNOG_COG_calls():
	f = open(deepnog_COG_call_file)
	out = []
	num_total_row = 0
	num_confident_row = 0
	for line in f.readlines()[1:]:
		num_total_row += 1
		columns = line.split(',')
		if float(columns[2]) >= confidence_level:
			out.append(columns)
			num_confident_row += 1
	print("finished reading prediction")
	print(f"out of {num_total_row} ORFs, {num_confident_row} has a confidence level above {confidence_level}.")
	return out
# deepNOG specific stuff ends


# diamond specific stuff begins

# The diamond database is in /researchdrive/zhongj2/cog-20.dmnd
'''
return a dict with:

key: gene id (e.g. NP_682421.1 -> NP_682421_1)
value: COG number (e.g. COG5718)
'''
def get_id_COG_dict():
	id_COG_dict = {}
	f = open(COG_ref) # COG_ref: see def in the top
	for line in f.readlines():
		r = line.split(',')
		gene_id, COG = r[2], r[6]
		gene_id = gene_id.replace(".", "_")
		id_COG_dict[gene_id] = COG
	return id_COG_dict

'''
filter diamond alignment:
input: diamond-out.tsv -> align our representive genome with the COG2020 database
return: a 2D python array -> meet our requirements and with unique match with the lowest e-value

And will write this 2D array into diamond-out-unique.tsv 
'''
def filter_diamond_output(BLAST_output):
	outfile = "tmp/diamond-out-unique.csv"

	if os.path.isfile(outfile):
		print(f"we already did the cleaning and writes to {outfile}. Just read the file and be done")
		return list(csv.reader(open(outfile)))

	# actually doing the cleaning
	BLAST_default_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
	expected_num_column = len(BLAST_default_columns)
	num_total_row = 0

	f = open(BLAST_output)
	ORF_seen = {} # key: ORF name; value: a BLAST match of this ORF
	for line in f.readlines():
		num_total_row += 1
		if num_total_row % 330000 == 0:
			print(f"counting the {num_total_row} line")
		r = line.split() # r for row
		if len(r) <= 1:
			pass
		if len(r) != expected_num_column:
			print("The number of columns in diamond blast output is not what I expected")
			exit(1)
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = r
		query_percent_alignment = (int(qend) - int(qstart))/int(length)
		target_percent_alignment = (int(send) - int(sstart))/int(length)
		if float(pident) > 30 and query_percent_alignment > 0.3 and target_percent_alignment > 0.3:
			if qseqid not in ORF_seen:
				ORF_seen[qseqid] = r
			else:
				cur_eval = ORF_seen[qseqid][-2] # look for the match with the lowest e-value
				if float(evalue) < float(cur_eval):
					ORF_seen[qseqid] = r
	
	best_BLAST_hits = list(ORF_seen.values())
	with open(outfile, "w") as f:
		writer = csv.writer(f)
		writer.writerows(best_BLAST_hits)
	return best_BLAST_hits

'''
filter diamond alignment:
diamond-out.tsv -> align our representive genome with the COG2020 database
'''
def open_diamond_COG_calls():
	outfile = "tmp/diamond_ORFid_COG_evalue.csv"

	if os.path.isfile(outfile):
		print(f"we figured out {outfile}. Just read the file and be done")
		return list(csv.reader(open(outfile)))

	# doing the heavy lifting
	gene_id_dict = get_id_COG_dict()
	print("finish opening the gene_id:COG dictionary")

	unique_diamond_output = filter_diamond_output(diamond_blast_output)

	out = []
	for r in unique_diamond_output:
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore = r
		cog_db_gene_id = sseqid.replace(".","_")
		COG = gene_id_dict[cog_db_gene_id]
		out.append([qseqid, COG, evalue])

	with open(outfile, "w") as f:
		writer = csv.writer(f)
		writer.writerows(out)
	return out

def main():
	ORF_calling_method = sys.argv[1]
	count_COGs_calls(ORF_calling_method)

if __name__ == "__main__":
	main()
