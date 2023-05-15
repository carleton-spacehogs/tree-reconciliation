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
num_seq_min=20
alignment_len_min=100

gene_tree_method=iqtree
COG_calling_method=diamond

if [ "$clock_model_input" = "ugam1" ] || [ "$clock_model_input" = "cir1" ] || [ "$clock_model_input" = "ln3" ]; then
	echo selecting the clock $clock_model_input
	clock_model=$clock_model_input
else
  echo "\"$clock_model_input\" is not a valid clock model. Options are ugam1, cir1, or ln3"
  clock_model="ugam1" # default
  echo "defaulting to $clock_model"
fi

chronogram="${clock_model}_ChenParamsEarth_sample.chronogram"
e_output="${clock_model}_ecceTERA_output"
e_analysis="${clock_model}_ecceTERA_analysis"

require_folders="$e_analysis $e_output gene_alignments iqtree_gene_trees tmp R-plots R-plots/histogram R-plots/timeline"
old_internal_node_f=${chronogram}_internal_nodes.txt

all_seq_fasta=SingleLine_EnrichedGenomes.faa
COG_summary="COG_reconciliation_summary_${clock_model}.csv"

conda_env_base="/Accounts/zhongj2/miniconda3"
conda_sh=${conda_env_base}/etc/profile.d/conda.sh
conda_R_env=${conda_env_base}/envs/anvio-dev

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

export gene_name gene_seq_file all_seq_fasta pre_trim trimv2 sym_event_f sym_event_date_f iqtree_log COG_summary clock_model chronogram e_output e_analysis COG_calling_method gene_tree_method ecceTERA_sym # for the python/R scripts
