#!/bin/bash

# COG5403,30
gene=COG5403

run_and_test() {
	$ccommand=$1
	$required_file=$2
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

source ./scripts/declare_file_location.sh
source ./scripts/declare_file_location.sh $gene

t1="./reconcile.sh --stop_before_tree_building --COG $gene"
run_and_test $t1 $alignment

t2="./reconcile.sh --stop_before_tree_building --gene_sequences tmp/$gene.faa"
run_and_test $t2 $alignment

t3="./reconcile.sh --stop_before_reconciliation --alignment $alignment"
run_and_test $t3 $iqtree_ufboot

t4="./reconcile.sh --gene_tree $iqtree_ufboot"
run_and_test $t3 $R_plot

# # start over and do it again as a whole pipeline
rm $trimv2 $iqtree_ufboot $store_gene_tree_filename $ecceTERA_sym $chronogram_internal_nodes_f $sym_event_f $sym_event_date_f $R_plot

t5="./reconcile.sh --COG $gene"
run_and_test $t5 $R_plot

