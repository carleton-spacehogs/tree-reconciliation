#!/bin/bash

gene_name=$1

if [ -z $gene_name ]; then
	echo declaring the hard-coded variables.
	# used in scirpts/utils.sh -> validate_required_folders
	num_core=10
	num_seq_min=20
	alignment_len_min=100

	require_folders="ecceTERA_analysis ecceTERA_output gene_alignments iqtree_gene_trees tmp R-plots R-plots/histogram R-plots/timeline"
	chronogram=ugam1_ChenParamsEarth_sample.chronogram
	old_internal_node_f=${chronogram}_internal_nodes.txt

	all_seq_fasta=SingleLine_EnrichedGenomes.faa
	COG_summary="COG_reconciliation_summary_v2.csv"
	gene_tree_method=iqtree
	COG_calling_method=diamond

	conda_env_base="/Accounts/zhongj2/miniconda3"
	conda_sh=${conda_env_base}/etc/profile.d/conda.sh
	conda_R_env=${conda_env_base}/envs/anvio-dev
else
	echo declaring default gene_name specific filenames and variables
	gene_seq_file=tmp/$gene_name.faa

	pre_trim=gene_alignments/${gene_name}.afa
	trimv1=gene_alignments/trimv1_$gene_name.afa
	trimv2=gene_alignments/trimv2_$gene_name.afa
	# store_gene_tree_filename="tmp/${gene_name}_${gene_tree_method}_gene_tree_filename.txt"

	outfile_prefix=iqtree_gene_trees/$gene_name
	iqtree_ufboot=${outfile_prefix}.ufboot
	iqtree_log=${outfile_prefix}.log

	old_e_sym="reconciliationsFile_canonical_symmetric.recPhyloXML"
	old_e_asym="reconciliationsFile_canonical_asymmetric.recPhyloXML"
	old_e_ran="reconciliationsFile_canonical_random.recPhyloXML"

	ecceTERA_sym="ecceTERA_output/${gene_name}_symmetric_reconciliation.recPhyloXML"
	ecceTERA_asym="ecceTERA_output/${gene_name}_asymmetric_reconciliation.recPhyloXML"
	ecceTERA_ran="ecceTERA_output/${gene_name}_random_reconciliation.recPhyloXML"

	chronogram_internal_nodes_f=ecceTERA_analysis/${chronogram%%_sample.chronogram}_${gene_name}_internal_nodes.txt # %% remove suffix
	sym_event_f=ecceTERA_analysis/${gene_name}_symmetric.events
	sym_event_date_f=ecceTERA_analysis/${gene_name}_symmetric.events_event_dates.txt

	# for validation purpose
	R_plot="R-plots/histogram/${COG_calling_method}-${gene_tree_method}-${gene_name}-eventsHistogram.png"
fi

export gene_name gene_seq_file all_seq_fasta pre_trim trimv2 sym_event_f sym_event_date_f iqtree_log COG_summary gene_tree_method COG_calling_method # for the python/R scripts
