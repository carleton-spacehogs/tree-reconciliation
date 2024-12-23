#!/bin/bash

clock_model=$1

if [ -z $clock_model]; then
	echo You must select one of the following clock model: ugam1, cir1, ln3
	echo "For example:
	./check-summary.sh ugam1
exiting
"
	exit 1
fi


source scripts/declare_file_location.sh --clock_model $clock_model
source scripts/utils.sh

echo $COG_summary
validate_required_folders

echo "Checkpoint 1: fixing missing summaries in $COG_summary"

running_COGs=$(ps aux | egrep -i "iqtree|ecceTERA" | grep -oP 'COG\d+' | uniq)

for afa in $(ls gene_alignments/COG*.afa); do
	COG=$(echo $afa | grep -oP 'COG\d+')
	source scripts/declare_file_location.sh --clock_model $clock_model --gene_name $COG
	echo $running_COGs | grep $COG > /dev/null
	if [ $? -ne 0 ]; then # it is not running
		res=$(grep $COG $COG_summary)
		if [ $? -ne 0 ]; then
			echo "No record for $COG, but $afa exists"
			if test -f $ecceTERA_sym; then
				echo "Recording for $COG as an ecceTERA success"
				scripts/summarize_reconciliation.py $COG 0
			elif test -f $iqtree_ufboot; then
				# echo "Recording for $COG as an ecceTERA fail"
				# scripts/summarize_reconciliation.py $COG 1
				./reconcile.sh --gene_tree $iqtree_ufboot --clock_model $clock_model
			else
				validate_res=$(validate_alignment $trimv2 $num_seq_min $alignment_len_min 1)
				echo $validate_res
				echo $validate_res | grep "validated" > /dev/null
				if [ $? -eq 0 ]; then
					echo "
$trimv2 ---> can proceed to tree building, but
$iqtree_ufboot ---> is missing (incomplete run)
consider doing:
./reconcile.sh --alignment $trimv2 --clock_model $1
"
				fi
			fi
		# else
		# 	echo "found : $res"
		fi
	else
		echo $COG is running, skip it
	fi
done

echo "Checkpoint 1: done"
echo "Checkpoint 2: whether there are missing graphs with ecceTERA output"

success_COGs=$(grep success $COG_summary | awk -F "," '{print $1}')

for COG in $success_COGs; do
	source scripts/declare_file_location.sh --clock_model $clock_model --gene_name $COG
	if [ ! -f $R_plot ]; then
		echo I think we are missing $R_plot, remaking it
		./reconcile.sh --gen_graph_only $COG --clock_model $clock_model
	fi
done

echo "Checkpoint 2: done"



