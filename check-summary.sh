#!/bin/bash

source scripts/declare_file_location.sh
source scripts/utils.sh

echo $COG_summary
validate_required_folders

echo "Checkpoint 1: fixing the missing summaries in $COG_summary"

running_COGs=$(ps aux | egrep -i "iqtree|ecceTERA" | grep -oP 'COG\d+' | uniq)

for afa in $(ls gene_alignments/COG*.afa); do
	COG=$(echo $afa | grep -oP 'COG\d+')
	source scripts/declare_file_location.sh $COG
	echo $$running_COGs | grep $COG > /dev/null
	if [ $? -ne 0 ]; then # it is not running
		res=$(grep $COG $COG_summary)
		if [ $? -ne 0 ]; then
			echo "we don't have the record for $COG, but I think we ran it before"
			scripts/declare_file_location.sh $COG
			if test -f $R_plot; then
				echo "Recording for $COG as an ecceTERA success"
				scripts/summarize_reconciliation.py $COG 0
			elif test -f $iqtree_ufboot; then
				echo "Recording for $COG as an ecceTERA fail"
				scripts/summarize_reconciliation.py $COG 1
			else
				validate_alignment $trimv2 $num_seq_min $alignment_len_min 1
			fi
		else
			echo found: $res
		fi
	else
		echo $COG is running, skip it
	fi
done

echo "Checkpoint 1: done"
echo "Checkpoint 2: whether there are missing graphs with ecceTERA output"

for graph in $(ls R-plots/histogram); do
	COG=$(echo $graph | awk -F"-" '{print $3}')
	grep $COG $COG_summary > /dev/null
	if [ $? -ne 0 ]; then
		echo $COG has the histogram:
		ls -lh R-plots/histogram/$graph
		echo but $COG not in $COG_summary, so I re-do the analysis
		# ./reconcile.sh --gen_graph_only $COG
	fi
done

echo "Checkpoint 2: done"
