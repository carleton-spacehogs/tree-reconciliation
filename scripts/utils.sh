#!/bin/bash

validate_required_folders() {
	require_folders="ecceTERA_analysis ecceTERA_output gene_alignments iqtree_gene_trees tmp R-plots R-plots/histogram R-plots/timeline"
	for folder in $require_folders; do
		if [ ! -d $folder ]; then
			echo $folder does not exist, I am creating it.
			mkdir $folder
		fi
	done
	
	COG_summary="COG_reconciliation_summary.csv"
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
		outfile_prefix=iqtree_gene_trees/$gene_name
		echo ${outfile_prefix}.ufboot > $gene_tree_f
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
	gene_name=$3

	ecceTERA_sym="ecceTERA_output/${gene_name}_symmetric_reconciliation.recPhyloXML"
	if test -f $ecceTERA_sym; then
		echo "$ecceTERA_sym exists, use the old symmetric reconciliation"
	else
		ecceTERA species.file=$spcies_tree gene.file=$gene_tree verbose=true print.reconciliations=1 recPhyloXML.reconciliation=true amalgamate=true
		for type in asymmetric random symmetric; do
			old_f=reconciliationsFile_canonical_${type}.recPhyloXML
			if test -f $old_f; then
				mv $old_f ecceTERA_output/${gene_name}_${type}_reconciliation.recPhyloXML
			else
				# if for some reason ecceTERA fails, stop here
				echo I do not have ---- $old_f ----, I will stop here
				exit 52
			fi
		done
	fi
}


inhouse_scripts_processing()
{
	chronogram=$1
	gene_name=$2

	recPhyloXML_file=ecceTERA_output/${gene_name}_symmetric_reconciliation.recPhyloXML
	chronogram_internal_nodes_f=ecceTERA_analysis/${chronogram%%_sample.chronogram}_${gene_name}_internal_nodes.txt # %% remove suffix
	sym_event_f=ecceTERA_analysis/${gene_name}_symmetric.events
	sym_event_date_f=ecceTERA_analysis/${gene_name}_symmetric.events_event_dates.txt
	birthdate_f=allGeneBirthDate.txt

	python3 scripts/getInternalNodeDates.py $chronogram $recPhyloXML_file
	# produce file: ugam1_ChenParamsEarth_sample.chronogram_internal_nodes.txt
	mv ${chronogram}_internal_nodes.txt $chronogram_internal_nodes_f

	python3 scripts/py3-scripts/recPhyloXMLEventSummary.py -i $recPhyloXML_file -o $sym_event_f --include.transfer.departure
	sed -r -i "s|\s+|\t|g" $sym_event_f

	python3 scripts/recphyloxmlinterpreterspecV2.py $sym_event_f $chronogram_internal_nodes_f $recPhyloXML_file

	python3 scripts/summarize_reconciliation.py $gene_name # output to COG_reconciliation_summary.csv
}



