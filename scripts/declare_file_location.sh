#!/bin/bash

options=$(getopt -o g:c: --long gene_name:,clock_model: -- "$@")
eval set -- "$options"

while true; do
	case "$1" in
		-g|--gene_name) # COG2801
			gene_name="$2"
			shift 2
			;;
		-c|--clock_model) # "ugam1, cir1, or ln3"
			clock_model_input="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Invalid option: $1" >&2
			exit 1
			;;
	esac
done

# echo declaring the hard-coded variables.
num_core=10
num_seq_min=15
alignment_len_min=100

# minimum 30% identity, maximum e-value of 10-5, minimum 70% subject and query alignment
min_seq_identity=30
min_percent_alignment=70

gene_tree_method=iqtree
COG_calling_method=diamond
clock_model="ugam1" # default

if [ "$clock_model_input" = "ugam1" ] || [ "$clock_model_input" = "cir1" ] || [ "$clock_model_input" = "ln3" ]; then
	# echo selecting the clock $clock_model_input
	clock_model=$clock_model_input
else
	echo "\"$clock_model_input\" is not a valid clock model. Options are ugam1, cir1, or ln3"
	echo "defaulting to $clock_model"
fi

chronogram="${clock_model}_ChenParamsEarth_sample.chronogram"
old_internal_node_f=${chronogram}_internal_nodes.txt
e_output="${clock_model}_ecceTERA_output"
e_analysis="${clock_model}_ecceTERA_analysis"
COG_summary="COG_reconciliation_summary_${clock_model}.csv"

conda_env=/workspace/data/Space_Hogs_shared_workspace/env_tree_reconciliation

large_file_dir=${conda_env}/jimmy_files
all_seq_fasta=${large_file_dir}/SingleLine_EnrichedGenomes.faa
diamond_COG_match=${large_file_dir}/diamond-out.tsv
deepNOG_COG_match=${large_file_dir}/SingleLine_deepnog.csv
COG_ref=${large_file_dir}/cog-20.cog.csv

for f in $chronogram $all_seq_fasta $diamond_COG_match $deepNOG_COG_match; do
	if [ ! -f $f ]; then
		echo I cannot find the file: $f, 
		echo which is required for running this program, exiting
		exit 1
	fi
done

require_folders="$e_analysis $e_output gene_alignments iqtree_gene_trees tmp R-plots R-plots/histogram R-plots/timeline R-plots/topBottom_histogram"

if [ ! -z $gene_name ]; then
	# echo declaring default filenames for : $gene_name
	gene_seq_file=tmp/$gene_name.faa

	pre_trim=gene_alignments/${gene_name}.afa
	trimv1=gene_alignments/trimv1_$gene_name.afa
	trimv2=gene_alignments/trimv2_$gene_name.afa

	outfile_prefix=iqtree_gene_trees/$gene_name
	iqtree_ufboot=${outfile_prefix}.ufboot
	iqtree_log=${outfile_prefix}.log

	old_e_sym="reconciliationsFile_canonical_symmetric.recPhyloXML"
	old_e_asym="reconciliationsFile_canonical_asymmetric.recPhyloXML"
	old_e_ran="reconciliationsFile_canonical_random.recPhyloXML"

	ecceTERA_sym="${e_output}/${gene_name}_symmetric_reconciliation.recPhyloXML"
	ecceTERA_asym="${e_output}/${gene_name}_asymmetric_reconciliation.recPhyloXML"
	ecceTERA_ran="${e_output}/${gene_name}_random_reconciliation.recPhyloXML"

	chronogram_internal_nodes_f=${e_analysis}/${chronogram%%_sample.chronogram}_${gene_name}_internal_nodes.txt # %% remove suffix
	sym_event_f=${e_analysis}/${gene_name}_symmetric.events
	sym_event_date_f=${e_analysis}/${gene_name}_symmetric.events_event_dates.txt

	# for validation purpose
	R_plot="R-plots/histogram/${gene_name}-${clock_model}-eventsHistogram.png"
fi

export min_seq_identity min_percent_alignment gene_name conda_env gene_seq_file all_seq_fasta diamond_COG_match deepNOG_COG_match COG_ref pre_trim trimv2 sym_event_f sym_event_date_f iqtree_log COG_summary clock_model chronogram e_output e_analysis COG_calling_method gene_tree_method ecceTERA_sym # for the python/R scripts
