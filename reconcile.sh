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
	e.g.: $(pwd)/reconcile.sh --COG COG0035

--gene_sequence [amino acids sequences in fasta format]
	e.g.: $(pwd)/reconcile.sh --gene rubisco_gene_sequences.faa

--alignment [mulitple_sequence_alignment.afa]
	e.g.: $(pwd)/reconcile.sh --alignment alignment.afa

--gene_tree [a collection of newick gene trees; output of iqtree or RAXML]
	e.g.: $(pwd)/reconcile.sh --gene_tree COG0035_trees.ufboot

--gen_graph_only [gene_name]

Other options:
--stop_before_reconciliation
	does everything before the ecceTERA step (memory intense)

--stop_before_tree_building

--species_tree
	provide path to your own species tree.
	Default: ugam1_ChenParamsEarth_sample.chronogram

-h --help : to display this message
"
}

options=$(getopt -o c:gg:gt:s:a:gs:h --long COG:,gen_graph_only:,gene_tree:,species_tree:,alignment:,gene_sequences:,stop_before_reconciliation,stop_before_tree_building,overwrite,help -- "$@")
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
		-gg|--gen_graph_only) # followed by gene_name
			gen_graph_only="$2"
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
		-gs|--gene_sequences)
			gene_sequences="$2"
			shift 2
			;;
		--stop_before_reconciliation)
			stop_before_reconciliation=true
			shift
			;;
		--stop_before_tree_building)
			stop_before_tree_building=true
			shift
			;;
		--overwrite)
			overwrite=true
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

parseGeneName_and_declareFilenames() {
	# e.g. iqtree_gene_trees/COG0774.ufboot -> COG0774
	gene_name=$(echo $1 | awk -F '/' '{print $NF}' | awk -F '.' '{print $1}')
	source ./scripts/declare_file_location.sh $gene_name
	echo $gene_name
}


analysis() {
	gene_name=$1
	source $conda_sh
	conda activate $conda_env_base
	inhouse_scripts_processing $chronogram $gene_name
	# inhouse_scripts_processing must executed in Jimmy's conda base environment
	# XML has some version changes, only the version in my base environment worked...

	# currently, R is broken in my base environment
	conda activate $conda_R_env
	Rscript --vanilla scripts/plot-gene-events-histogram.R
	Rscript --vanilla scripts/plot-gene-timeline.R
	conda deactivate
}


reconcile_and_analysis() {
	gene_tree=$1
	echo "
	gene tree: $gene_tree
	species tree: $chronogram
	now reconciling them
	"
	run_ecceTERA $chronogram $gene_tree

	if [ $? -ne 0 ]; then
		echo ecceTERA failed... Stop here
		exit 1
	fi

	analysis $gene_name
}


makeTree_and_reconcile() {
	alignment=$1
	# gene_name=$(parse_gene_name $alignment)

	echo trimming the original alignment $alignment
	# trim out positions with mostly gaps
	~/miniconda3/bin/trimal -in $alignment -out $trimv1 -gt 0.15 # -automated1
	
	# we excluded genes with more than 20% gaps in the trimmed alignment
	python3 scripts/keep_seq_geq_Xpercent.py $trimv1 $trimv2 80
	rm $trimv1

	validate_alignment $trimv2 $num_seq_min $alignment_len_min
	
	if [ "$stop_before_tree_building" = true ]; then
		echo I see --stop_before_tree_building, so I stop
		exit $?
	fi
	
	# we store the filename of the outputted gene tree here:
	# store_gene_tree_filename="tmp/${gene_name}_${gene_tree_method}_gene_tree_filename.txt"
	generate_gene_tree $gene_tree_method $trimv2 $gene_name $num_core $store_gene_tree_filename

	if [ "$stop_before_reconciliation" = true ]; then
		echo I see --stop_before_reconciliation, so I stop
		exit $?
	fi

	gene_tree=$(cat $store_gene_tree_filename)
	reconcile_and_analysis $gene_tree
}


align_makeTree_and_reconcile() {
	seq_file=$1 # gene_sequences
	echo aligning this gene: $gene_name, and alignment output: $alignment
	# muscle is localed in this folder, executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
	./muscle5.1 -align $seq_file -output $alignment

	if [ $? -ne 0 ]; then
		echo alignment of $seq_file not successful, stop here
		exit 1
	fi

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

	echo $all_seq_fasta
	python3 scripts/grep_seq_given_COG.py $COG $COG_calling_method $gene_seq_file $all_seq_fasta

	if [ $? -ne 0 ]; then
		echo grabbing the sequence for $COG not successful, stop here
		rm $gene_seq_file
		exit 1
	fi

	align_makeTree_and_reconcile $gene_seq_file
}

am_I_done() {
	if [ "$overwrite" = true ]; then
		echo I am over writing stuff
	elif [ -f "$R_plot" ]; then
		echo "$R_plot exists. I think I am done"
		exit 0
	fi
}

source scripts/declare_file_location.sh # declare hard-coded variables
# gene specific variables declared in parseGeneName_and_declareFilenames
validate_required_folders

# if given both gene and species tree, just do the reconciliation
if [ ! -z "$gen_graph_only" ]; then
	source scripts/declare_file_location.sh $gen_graph_only
	analysis $gen_graph_only
	exit $?
elif [ ! -z "$gene_tree" ]; then
	parseGeneName_and_declareFilenames $gene_tree
	am_I_done
	reconcile_and_analysis $gene_tree
	exit $?
elif [ ! -z "$alignment" ]; then
	parseGeneName_and_declareFilenames $alignment
	am_I_done
	makeTree_and_reconcile $alignment
	exit $?
elif [ ! -z "$gene_sequences" ]; then
	parseGeneName_and_declareFilenames $gene_sequences
	am_I_done
	align_makeTree_and_reconcile $gene_sequences
	exit $?
elif [ ! -z "$COG" ]; then
	source ./scripts/declare_file_location.sh $COG
	am_I_done
	extractSeq_align_makeTree_and_reconcile $COG
	exit $?
else
	echo You have to gimme an option to work with. Need help???
	echo
	Help
fi

