#!/bin/python3
import sys
import os
import csv
import datetime
from Bio import SeqIO

'''
Jimmy (Juntao) Zhong, Spring 2023, under Professor Rika Anderson at Carleton College

This script take a COG number, and record the summary statistics of this alignment.
It can only be ran after the reconciliation is finished.

example (if record of COG4100 already exist the script does nothing):
./summarize_reconciliation.py COG4100 

If you want to overwrite the old record of COG4100, do:
./summarize_reconciliation.py COG4100 rewrite
'''

summary_file = "COG_reconciliation_summary.csv" # with the following columns:
# COG_name,num_ORFs_begin,num_ORFs_in_tree,raw_alignment_length,trimmed_alignment_length,average_ORF_length,num_events

COG = sys.argv[1]
rewrite = True if len(sys.argv) == 3 and sys.argv[2] == "rewrite" else False

original_align = f"gene_alignments/{COG}.afa"
trimmed_align = f"gene_alignments/trimv2_{COG}.afa"
event_dates_f = f"ecceTERA_analysis/{COG}_symmetric.events_event_dates.txt"
iqtree_log = f"iqtree_gene_trees/{COG}.log"

def record_exists(COG, summary):
    for l in summary:
        if l[0] == COG:
            return True

def check_required_files():
    for file in [original_align, trimmed_align, event_dates_f, iqtree_log]:
        if not os.path.isfile:
            print(f"I need this file: {file}, but it is not there")
            exit(1)

def get_alignment_length(alignment_dict):
    for val in alignment_dict.values():
        return len(val.seq)

def get_average_gene_len(alignment_dict):
    total_bp = 0
    for val in alignment_dict.values():
        aa_sequence = str(val.seq).replace("-", "")
        aa_sequence = aa_sequence.replace("*", "")
        total_bp += len(aa_sequence)
    return total_bp/len(alignment_dict)

def count_events(event_dates_f):
    num_events = -1 # excluding the first line column names
    with open(event_dates_f) as events:
        for event in events:
            if "?" not in event:
                num_events += 1
    return num_events

def get_earliest_event(event_dates_f):
    # inspired by getGeneBirthDate.py [symmetric.events_event_dates.txt file] from Garrett Chappell, Summer 2021
    earliest_date_line = []
    curr_earliest_date = 0
    with open(event_dates_f) as events:
        for event in events:
            if "?" not in event and "event" not in event: # "event" for the first line; "?" for null data
                col = event.split('\t')
                midpoint_date = float(col[5].strip())
                if midpoint_date > curr_earliest_date:
                    curr_earliest_date = midpoint_date
                    earliest_date_line = col
    event_type, left_node, right_node, left_date, right_date, mid_date = earliest_date_line
    return [event_type, left_date, mid_date.strip(), right_date]

def get_num_reconciliation(iqtree_log):
    msg = "Total number of iterations: "
    with open(iqtree_log) as log:
        for l in log:
            if msg in l:
                return int(l.replace(msg, ""))

if __name__ == "__main__":
    summary = list(csv.reader(open(summary_file)))
    if (not rewrite) and record_exists(COG, summary):
        print("record exists, not doing anything. Use the 'rewrite' to overwrite old record")
        exit(0)

    check_required_files()

    original_align_dict = SeqIO.to_dict(SeqIO.parse(original_align, "fasta"))
    trimmed_align_dict = SeqIO.to_dict(SeqIO.parse(trimmed_align, "fasta"))

    num_ORFs_begin = len(original_align_dict)
    num_ORFs_in_tree = len(trimmed_align_dict)
    raw_alignment_length = get_alignment_length(original_align_dict)
    trimmed_alignment_length = get_alignment_length(trimmed_align_dict)
    average_ORF_length = get_average_gene_len(original_align_dict)
    num_events = count_events(event_dates_f)
    num_reconciliation = get_num_reconciliation(iqtree_log)

    row = [COG, num_ORFs_begin, num_ORFs_in_tree, average_ORF_length, raw_alignment_length, trimmed_alignment_length, 
           num_events, num_reconciliation, datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")]
    row += get_earliest_event(event_dates_f)
    row = [str(x) for x in row]
    print(row)

    if rewrite and record_exists(COG, summary):
        for line_num in range(len(summary)):
            if summary[line_num][0] == COG:
                summary[line_num] = row
        with open(summary_file, 'w') as s:
            wr = csv.writer(s)
            wr.writerows(summary)
    else:
        with open(summary_file, 'a') as s:
            write_row = ",".join(row) + "\n"
            s.write(write_row)
