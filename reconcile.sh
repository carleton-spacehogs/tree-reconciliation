#!/bin/bash

source ./scripts/utils.sh

Help()
{
	# Display Help
	echo
	echo "Jimmy's effort to automate tree conciliation with ecceTERA
working under Professor Rika Anderson in Carleton College, MN. Started March 2023

choose only one of the entry point to start this program

--COG [COG_number]
	e.g.: $(pwd) --COG COG0035

--gene_sequence [amino acids sequences in fasta format]
	e.g.: $(pwd) --gene rubisco_gene_sequences.faa

--alignment [mulitple_sequence_alignment.afa]
	e.g.: $(pwd) --alignment alignment.afa

--gene_trees [a collection of newick gene trees; output of iqtree or RAXML]
	e.g.: $(pwd) --gene_trees COG0035_trees.ufboot

COG:,gene_tree:,species_tree:,alignment:,gene_sequence


Other options:
--stop_before_reconciliation
	does everything before the ecceTERA step (memory intense)

--species_tree
	provide path to your own species tree.
	Default: ugam1_ChenParamsEarth_sample.chronogram

-h --help : to display this message
"
}


chronogram=ugam1_ChenParamsEarth_sample.chronogram
all_seq_fasta=SingleLine_EnrichedGenomes.faa
num_core=10
gene_tree_method=iqtree
COG_calling_method=diamond


options=$(getopt -o c:gt:s:a:gs:h --long COG:,gene_tree:,species_tree:,alignment:,gene_sequence:,stop_before_reconciliation,help -- "$@")
eval set -- "$options"

while true; do
	case "$1" in
		-h|--help)
			Help
			exit
			;;
		-c|--COG) # COG2801
			COG="$2"
			shift 2
			;;
		-gt|--gene_tree)
			gene_tree="$2"
			shift 2
			;;
		-s|--spcies_tree)
			chronogram="$2"
			shift 2
			;;
		-a|--alignment)
			alignment="$2"
			shift 2
			;;
		-gs|--gene_sequence)
			gene_sequence="$2"
			shift 2
			;;
		--stop_before_reconciliation)
			stop_before_reconciliation=true
			shift
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


parse_gene_name() { echo $1 | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}'; }


reconcile_and_analysis() {
	gene_tree=$1
	echo "
	gene tree: $gene_tree
	species tree: $chronogram
	now reconciling them
	"

	gene_name=$(parse_gene_name $gene_tree)
	# e.g. gene_tree=iqtree_gene_trees/COG0774.ufboot
	# ---> gene_name=COG0774

	run_ecceTERA $chronogram $gene_tree $gene_name
	inhouse_scripts_processing $chronogram $gene_name
	# inhouse_scripts_processing must executed in Jimmy's conda base environment
	# XML has some version changes, only the version in my base environment worked...

	# currently, R is broken in my base environment
	source activate anvio-dev
	Rscript --vanilla plot-gene-events-histogram.R $gene_name $gene_tree_method $COG_calling_method
	Rscript --vanilla plot-gene-timeline.R $gene_name $gene_tree_method $COG_calling_method
	source deactivate
}


makeTree_and_reconcile() {
	alignment=$1
	gene_name=$(parse_gene_name $alignment)

	echo trimming the original alignment $alignment
	tmp_alignment=tmp/$gene_name.alignment
	# trim out positions with mostly gaps
	~/miniconda3/bin/trimal -in $alignment -out $tmp_alignment -gt 0.15 # -automated1
	
	# we excluded genes with more than 20% gaps in the trimmed alignment
	trimmed_alignment=gene_alignments/trimmed_$COG.afa
	./keep_seq_geq_Xpercent.py $tmp_alignment $trimmed_alignment 80
	rm $tmp_alignment

	# we store the filename of the outputted gene tree here:
	store_gene_tree_filename="tmp/${gene_name}_${gene_tree_method}_gene_tree_filename.txt"
	generate_gene_tree $gene_tree_method $trimmed_alignment $gene_name $num_core $store_gene_tree_filename

	if [ "$stop_before_reconciliation" = true ]; then
		echo I see --stop_before_reconciliation, so I stop
		exit $?
	fi

	gene_tree=$(cat $store_gene_tree_filename)
	reconcile_and_analysis $gene_tree
}


align_makeTree_and_reconcile() {
	seq_file=$1 # gene_sequence
	gene_name=$(parse_gene_name ${seq_file})
	alignment=gene_alignments/${gene_name}.afa
	echo aligning this gene: $gene_name, and alignment output: $alignment
	# muscle is localed in this folder, executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
	./muscle5.1 -align $seq_file -output $alignment

	makeTree_and_reconcile $alignment
}


extractSeq_align_makeTree_and_reconcile() {
	COG=$1
	# sanity check
	echo $COG | grep COG
	if [[ $? != 0 ]]; then
		echo The COG number $COG must be in the format of "COGXXXX"
		exit 1
	fi

	gene_seq_file=tmp/$COG.faa
	python3 grep_seq_given_COG.py $COG $COG_calling_method $gene_seq_file $all_seq_fasta

	align_makeTree_and_reconcile $gene_seq_file
}


# if given both gene and species tree, just do the reconciliation
if [ ! -z "$gene_tree" ]; then
	reconcile_and_analysis $gene_tree
	exit $?
elif [ ! -z "$alignment" ]; then
	makeTree_and_reconcile $alignment
	exit $?
elif [ ! -z "$gene_sequence" ]; then
	align_makeTree_and_reconcile $gene_sequence
	exit $?
elif [ ! -z "$COG" ]; then
	extractSeq_align_makeTree_and_reconcile $COG
	exit $?
else
	echo You have to gimme an option to work with. Need help???
	echo
	Help
fi

