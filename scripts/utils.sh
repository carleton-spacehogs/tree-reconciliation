#!/bin/bash

validate_required_folders() {
	for folder in $require_folders; do
		if [ ! -d $folder ]; then
			echo $folder does not exist, I am creating it.
			mkdir $folder
		fi
	done

	if [ ! -f $COG_summary ]; then
		echo $COG_summary does not exist, creating a new file
		echo -en "COG,num_ORFs_begin,average_ORF_length,raw_alignment_length,num_ORFs_in_tree,trimmed_alignment_length,num_iterations,num_events,earliest_event_type,left_date,mid_date,right_date,exit_status,time_stamp\n" > $COG_summary
	fi
}

less_20_warning() {
	num_seq_min=$1
	num_seq_exist=$2
	alignment_for_tree=$3
	echo There are only $num_seq_exist ORFs exist in the alignment file $alignment_for_tree
	echo Tree reconciliation requires at least $num_seq_min sequences for analysis
	echo skipping
	./scripts/summarize_reconciliation.py $gene_name 2
}

validate_alignment() {
	alignment_for_tree=$1
	num_seq_min=$2
	alignment_len_min=$3
	do_not_exit=$4
	if test -f $alignment_for_tree; then
		num_seq_exist=$(grep -c ">" $alignment_for_tree)
		if (( $num_seq_min > $num_seq_exist )); then
			less_20_warning $num_seq_min $num_seq_exist $alignment_for_tree
			if [ $do_not_exit -ne 1 ]; then
				exit 0
			fi
		fi

		first_alignment_length=$( head -n 20 $alignment_for_tree | tail -n +2 | grep --max-count 1  -b '>' | cut -f1 -d:)
		if (( $alignment_len_min > $first_alignment_length - 1)); then
			echo The length of the alignment --- $alignment_for_tree --- is only $first_alignment_length
			echo We need the alignment length of at least $alignment_len_min amino acids
			echo skipping
			./scripts/summarize_reconciliation.py $gene_name 3
			if [ $do_not_exit -ne 1 ]; then
				exit 0
			fi
		fi
	else
		less_20_warning $num_seq_min $num_seq_exist $alignment_for_tree
		if [ $do_not_exit -ne 1 ]; then
			exit 0
		fi
	fi

	echo Alignment file --- $alignment_for_tree --- is validated
}

run_ecceTERA()
{
	spcies_tree=$1
	gene_tree=$2

	if test -f $ecceTERA_sym; then
		echo "$ecceTERA_sym exists, use the old symmetric reconciliation"
	else
		ecceTERA species.file=$spcies_tree gene.file=$gene_tree verbose=true print.reconciliations=1 recPhyloXML.reconciliation=true amalgamate=true
		ecceTERA_status=$?
		echo $ecceTERA_status,$gene_tree >> eceeTERA-exit-code.tmp
		if [ $ecceTERA_status -ne 0 ]; then
			echo ecceTERA failed exit status $ecceTERA_status
			echo skipping
			./scirpts/summarize_reconciliation.py $gene_name 1
			exit $ecceTERA_status
		else 
			mv $old_e_sym $ecceTERA_sym
			mv $old_e_asym $ecceTERA_asym
			mv $old_e_ran $ecceTERA_ran
		fi
	fi
}


inhouse_scripts_processing()
{
	python3 scripts/getInternalNodeDates.py $chronogram $ecceTERA_sym
	mv $old_internal_node_f $chronogram_internal_nodes_f

	python3 scripts/py3-scripts/recPhyloXMLEventSummary.py -i $ecceTERA_sym -o $sym_event_f --include.transfer.departure
	sed -r -i "s|\s+|\t|g" $sym_event_f

	python3 scripts/recphyloxmlinterpreterspecV2.py $sym_event_f $chronogram_internal_nodes_f $ecceTERA_sym

	python3 scripts/summarize_reconciliation.py $gene_name 0
	if [ $? -ne 0 ]; then
		echo scripts/summarize_reconciliation.py failed. stop here
		exit 1
	fi
}



