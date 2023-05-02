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
		echo -en "COG,num_ORFs_begin,num_ORFs_in_tree,average_ORF_length,raw_alignment_length,trimmed_alignment_length,num_events,num_iterations,time_stamp,earliest_event_type,left_date,mid_date,right_date\n" > $COG_summary
	fi
}

validate_alignment() {
	alignment_for_tree=$1
	num_seq_min=$2
	alignment_len_min=$3
	num_seq_exist=$(grep -c ">" $alignment_for_tree)
	if (( $num_seq_min > $num_seq_exist )); then
		echo There are only $num_seq_exist ORFs exist in the alignment file $alignment_for_tree
		echo Tree reconciliation requires at least $num_seq_min sequences for analysis
		echo skipping
		exit 0
	fi
	
	first_alignment_length=$( head -n 20 $alignment_for_tree | tail -n +2 | grep --max-count 1  -b '>' | cut -f1 -d:)
	if (( $alignment_len_min > $first_alignment_length - 1)); then
		echo The length of the alignment --- $alignment_for_tree --- is only $first_alignment_length
		echo We need the alignment length of at least $alignment_len_min amino acids
		echo skipping
		exit 0
	fi

	echo Alignment file --- $alignment_for_tree --- is validated
	# echo Alignment file --- $alignment_for_tree --- is validated >> tmp/tmp2.txt
	echo Jimmy building tree...
}

generate_gene_tree()
{
	method=$1 # raxml or iqtree
	trimmed_alignment=$2
	gene_name=$3
	num_core=$4
	gene_tree_f=$5

	echo using \"$method\" to make Max likelihood tree for \"$gene_name\"
	if [[ $method == "raxml" ]]; then
		echo raxml_gene_trees/RAxML_bipartitions.$gene_name.tree > $gene_tree_f
		# took it from baross: 
		# I have to write "-w $(pwd)/gene_trees -n ${gene_name}.tree" instead of "gene_trees/${gene_name}.tree"
		../raxmlHPC -f a -m PROTGAMMAAUTO -p 33 -x 33 -s $trimmed_alignment -w ../raxml_gene_trees -n $gene_name.tree -T $num_core -N 100 # num_bootstrap set to 100
	elif [[ $method == "iqtree" ]]; then
		which iqtree
		if [ $? -ne 0 ]; then
			echo iqtree is not usable in this environment.
			echo do \"which iqtree\" and see the problem
			exit 1
		fi
		echo $iqtree_ufboot > $gene_tree_f
		iqtree -s $trimmed_alignment --nmax 3000 -m MFP -bb 1000 -wbtl -ntmax 10 -nt AUTO -T $num_core --prefix $outfile_prefix -madd C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60 -mrate E,I,G,I+G,R
	else
		echo "ERROR!! \$method can only be \"iqtree\" or \"raxml\""
		exit 1
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
		if test -f $old_e_sym; then
			mv $old_e_sym $ecceTERA_sym
			mv $old_e_asym $ecceTERA_asym
			mv $old_e_ran $ecceTERA_ran
		else
			# if for some reason ecceTERA fails, stop here
			echo I do not have ---- $old_f ----, I will stop here
			exit 52
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

	python3 scripts/summarize_reconciliation.py $gene_name # output to COG_reconciliation_summary.csv
}



