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
	./scripts/summarize_reconciliation.py $gene_name 2
}

should_I_exit() {
	do_not_exit=$1
	if [ $do_not_exit -ne 1 ]; then
		echo skipping
		conda deactivate
		exit 1
	else
		echo 1
	fi
}

validate_alignment() {
	alignment_for_tree=$1
	num_seq_min=$2
	alignment_len_min=$3
	do_not_exit=$4
	if test -f $alignment_for_tree; then
		num_seq_exist=$(grep -c ">" $alignment_for_tree)
		first_alignment_length=$( head -n 20 $alignment_for_tree | tail -n +2 | grep --max-count 1  -b '>' | cut -f1 -d:)
		if (( $num_seq_min > $num_seq_exist )); then
			less_20_warning $num_seq_min $num_seq_exist $alignment_for_tree
			should_I_exit $do_not_exit
		elif (( $alignment_len_min > $first_alignment_length - 1)); then
			echo The length of the alignment --- $alignment_for_tree --- is only $first_alignment_length
			echo We need the alignment length of at least $alignment_len_min amino acids
			./scripts/summarize_reconciliation.py $gene_name 3
			should_I_exit $do_not_exit
		else
			echo $alignment_for_tree is validated
		fi
	else
		less_20_warning $num_seq_min $num_seq_exist $alignment_for_tree
		should_I_exit $do_not_exit
	fi
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
		if test -f $old_e_sym; then
			mv $old_e_sym $ecceTERA_sym
			mv $old_e_asym $ecceTERA_asym
			mv $old_e_ran $ecceTERA_ran
		else
			echo ecceTERA failed? ecceTERA exit status is : $ecceTERA_status
			echo $ecceTERA_status,$gene_tree >> eceeTERA-exit-code.tmp
			echo skipping
			./scripts/summarize_reconciliation.py $gene_name 1
			exit $ecceTERA_status
		fi
	fi
}


analysis()
{
	if [ ! -f "$sym_event_f" ]; then
		python3 scripts/getInternalNodeDates.py $chronogram $ecceTERA_sym
		mv $old_internal_node_f $chronogram_internal_nodes_f

		python3 scripts/py3-scripts/recPhyloXMLEventSummary.py -i $ecceTERA_sym -o $sym_event_f --include.transfer.departure
		sed -r -i "s|\s+|\t|g" $sym_event_f
	else
		echo I already have $sym_event_f, not doing it again.
	fi

	python3 scripts/recphyloxmlinterpreterspecV2.py $sym_event_f $chronogram_internal_nodes_f $ecceTERA_sym

	python3 scripts/summarize_reconciliation.py $gene_name 0 rewrite
	if [ $? -ne 0 ]; then
		echo scripts/summarize_reconciliation.py failed. stop here
		exit 1
	fi

	Rscript --vanilla scripts/plot-gene-events-histogram.R
	Rscript --vanilla scripts/plot-gene-events-histogram-topBottom.R
	Rscript --vanilla scripts/plot-gene-timeline.R
}



