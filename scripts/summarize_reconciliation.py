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

COG = sys.argv[1]
error_status = sys.argv[2] # 0: sucessful; 1: ecceTERA memory fail; 2: 20 < ORFs after trimming

rewrite = True if len(sys.argv) == 3 and sys.argv[2] == "rewrite" else False

if "pre_trim" not in os.environ.keys():
    print("the environment is not set up yet, do these 2 commands:")
    print("    source ./scripts/declare_file_location.sh")
    print(f"    source ./scripts/declare_file_location.sh {COG}")
    exit(1)

original_align = os.environ["pre_trim"]
trimmed_align = os.environ["trimv2"]
event_dates_f = os.environ["sym_event_date_f"]
iqtree_log = os.environ["iqtree_log"]
summary_file = os.environ["COG_summary"] # with the following columns:
# COG_name,num_ORFs_begin,num_ORFs_in_tree,raw_alignment_length,trimmed_alignment_length,average_ORF_length,num_events

def record_exists(COG, summary):
    for l in summary:
        if l[0] == COG:
            return True

def check_required_files():
    def print_and_exit(file):
        print(f"I need this file: {file}, but it is not there")
        exit(1)
    if not os.path.isfile(original_align):
        print_and_exit(original_align)
    if error_status == 0:
        if not os.path.isfile(event_dates_f):
            print_and_exit(event_dates_f)
    elif error_status == 1:
        for file in [trimmed_align, iqtree_log]:
            if not os.path.isfile(file):
                print_and_exit(file)
    print("I think have all files I need")

def get_alignment_length(alignment_dict):
    for val in alignment_dict.values():
        return len(val.seq)

def get_average_gene_len(alignment_dict):
    total_bp = 0
    for val in alignment_dict.values():
        aa_sequence = str(val.seq).replace("-", "")
        aa_sequence = aa_sequence.replace("*", "")
        total_bp += len(aa_sequence)
    if len(alignment_dict) == 0:
        return 0
    else:
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
                midpoint_date = 0
                if col[5].strip() != "NA":
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

def basic_summary(original_align):
    original_align_dict = SeqIO.to_dict(SeqIO.parse(original_align, "fasta"))
    num_ORFs_begin = len(original_align_dict)
    raw_alignment_length = get_alignment_length(original_align_dict)
    average_ORF_length = get_average_gene_len(original_align_dict)
    return [COG, num_ORFs_begin, average_ORF_length, raw_alignment_length]

def enough_ORFs_summary(trimmed_align, iqtree_log):
    trimmed_align_dict = SeqIO.to_dict(SeqIO.parse(trimmed_align, "fasta"))
    num_ORFs_in_tree = len(trimmed_align_dict)
    trimmed_alignment_length = get_alignment_length(trimmed_align_dict)
    num_reconciliation = get_num_reconciliation(iqtree_log)
    return [num_ORFs_in_tree, trimmed_alignment_length, num_reconciliation]

def ecceTERA_summary(event_dates_f):
    return [count_events(event_dates_f)] + get_earliest_event(event_dates_f)

if __name__ == "__main__":
    summary = list(csv.reader(open(summary_file)))
    if (not rewrite) and record_exists(COG, summary):
        print("record exists, not doing anything. Use the 'rewrite' to overwrite old record")
        exit(0)

    check_required_files()

    row = []
    if error_status == "0":
        row = basic_summary(original_align) + enough_ORFs_summary(trimmed_align, iqtree_log) + ecceTERA_summary(event_dates_f) + ["everything_successful"]
    elif error_status == "1":
        row = basic_summary(original_align) + enough_ORFs_summary(trimmed_align, iqtree_log) + ["NA"] * 5 + ["ecceTERA_failed"]
    elif error_status == "2":
        row = basic_summary(original_align) + ["NA"] * 8 + ["less_than_20_ORFs"]
    elif error_status == "3":
        row = basic_summary(original_align) + ["NA"] * 8 + ["less_than_100_alignment_len"]
    else:
        print("I only accept 3 error status: 0, 1, 2, 3")
        print("0: sucessful; 1: ecceTERA memory fail; 2: ORFs < 20 after trimming; 3: alignment length < 100 after trimming")
        exit(1)

    row.append(datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
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
