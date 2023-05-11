#!/bin/bash

# COG4833,20
# COG3364,20
# COG5733,21
# COG3331,23
# COG3388,23
# COG3305,24
# COG4337,25
# COG3312,25
# COG3433,26
# COG3397,26
# COG3356,26
# COG3365,26
# COG3133,26
# COG3317,29
# COG4338,30
# COG5338,30
# COG3310,31
# COG3341,32
# COG3337,32
# COG0833,32
# COG3933,32
# COG4333,32
# COG5702,33
# COG5746,33
# COG5559,33
# COG3188,33
# COG5771,33
# COG4240,33
# COG4521,33
# COG4154,33
# COG1590,33
# COG4916,33
# COG4246,33
# COG3277,33
# COG1665,33
# COG4766,33
# COG4933,33
# COG4397,33
# COG3357,33
# COG5925,33
# COG3069,33
# COG5029,33
# COG1103,33
# COG5397,33
# COG3037,33
# COG3297,33
gene1=COG4549
gene2=COG4310

clock_model="ln3"

run_and_test() {
	command="$1"
	required_file="$2"
	$command

	if [ $? -ne 0 ]; then
		echo the execution of:
		echo $command
		echo returns non-zero status. Something maybe wrong, exiting
		exit 1
	fi

	if test -f $required_file; then
		echo the execution of:
		echo $command
		echo looks correct
		echo
	else
		echo I did not find $required_file, which was expected from
		echo $command
		echo exiting
		exit 1
	fi
}

source ./scripts/declare_file_location.sh --clock_model $clock_model --gene_name $gene1

t1="./reconcile.sh --stop_before_tree_building --COG $gene1 --clock_model $clock_model"
run_and_test "$t1" $pre_trim

t2="./reconcile.sh --stop_before_tree_building --gene_sequences tmp/$gene1.faa --clock_model $clock_model"
run_and_test "$t2" $pre_trim

t3="./reconcile.sh --stop_before_reconciliation --alignment $pre_trim --clock_model $clock_model"
run_and_test "$t3" $iqtree_ufboot

t4="./reconcile.sh --gene_tree $iqtree_ufboot --clock_model $clock_model"
run_and_test "$t4" $R_plot

# # start over and do it again as a whole pipeline
# # rm $trimv2 $iqtree_ufboot $store_gene_tree_filename $ecceTERA_sym $chronogram_internal_nodes_f $sym_event_f $sym_event_date_f $R_plot

# Run the whole pipeline with another gene
source ./scripts/declare_file_location.sh --clock_model $clock_model --gene_name $gene2
t5="./reconcile.sh --COG $gene2 --clock_model $clock_model"
run_and_test "$t5" $R_plot

