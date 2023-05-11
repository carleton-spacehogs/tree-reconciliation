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

--clock_model
	select either ugam1, cir1, or ln3. If nothing is selected, ugam1 is used
	the chronogram file is set to:
		ugam1_ChenParamsEarth_sample.chronogram, or
		cir1_ChenParamsEarth_sample.chronogram, or
		ln3_ChenParamsEarth_sample.chronogram,
	accordingly

-h --help : to display this message
"
}

options=$(getopt -o c:gg:gt:s:a:gs:h --long COG:,gen_graph_only:,gene_tree:,clock_model:,alignment:,gene_sequences:,stop_before_reconciliation,stop_before_tree_building,overwrite,help -- "$@")
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
		-s|--clock_model)
			clock_model="$2"
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
	pre_trim=$1

	echo trimming the original alignment $pre_trim
	# trim out positions with mostly gaps
	~/miniconda3/bin/trimal -in $pre_trim -out $trimv1 -gt 0.15 # -automated1
	
	# we excluded genes with more than 20% gaps in the trimmed alignment
	python3 scripts/keep_seq_geq_Xpercent.py $trimv1 $trimv2 80
	rm $trimv1

	validate_alignment $trimv2 $num_seq_min $alignment_len_min

	if [ "$stop_before_tree_building" = true ]; then
		echo I see --stop_before_tree_building, so I stop
		exit 0
	fi

	which iqtree
	if [ $? -ne 0 ]; then
		echo can not use iqtree. Do \"which iqtree\" and see the problem.
		exit 1
	fi

	iqtree -s $trimv2 --nmax 3000 -m MFP -bb 1000 -wbtl -ntmax 10 -nt AUTO -T $num_core --prefix $outfile_prefix -madd "C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60" -mrate $rate_models "E,I,G,I+G,R" # List of 1. additional models 2. rate model for heterogeneity among sites

	if [ "$stop_before_reconciliation" = true ]; then
		echo I see --stop_before_reconciliation, so I stop
		exit 0
	fi

	reconcile_and_analysis $iqtree_ufboot
}


align_makeTree_and_reconcile() {
	seq_file=$1 # gene_sequences
	echo aligning this gene: $gene_name, and alignment output: $pre_trim
	# muscle is localed in this folder, executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
	./muscle5.1 -align $seq_file -output $pre_trim

	if [ $? -ne 0 ]; then
		echo alignment of $seq_file not successful, stop here
		exit 1
	fi

	makeTree_and_reconcile $pre_trim
}


extractSeq_align_makeTree_and_reconcile() {
	COG=$1
	# sanity check
	echo $COG | grep COG
	if [[ $? != 0 ]]; then
		echo The COG number $COG must be in the format of "COGXXXX"
		exit 1
	fi

	if test -f $iqtree_ufboot; then
		echo I found the gene_tree: skipping to the reconciliation.
		reconcile_and_analysis $iqtree_ufboot
	else
		echo $all_seq_fasta
		python3 scripts/grep_seq_given_COG.py $COG
		if [ $? -ne 0 ]; then
			echo grabbing the sequence for $COG not successful, stop here
			rm $gene_seq_file
			exit 1
		fi
		align_makeTree_and_reconcile $gene_seq_file
	fi
}

am_I_done() {
	if [ "$overwrite" = true ]; then
		echo I am over writing stuff
	elif [ -f "$R_plot" ]; then
		echo "$R_plot exists. I think I am done"
		python3 scripts/summarize_reconciliation.py $gene_name 0 rewrite
		exit 0
	fi
}

source scripts/declare_file_location.sh --clock_model $clock_model
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
	source ./scripts/declare_file_location.sh --clock_model $clock_model --gene_name $COG
	am_I_done
	extractSeq_align_makeTree_and_reconcile $COG
	exit $?
else
	echo You have to gimme an option to work with. Need help???
	echo
	Help
fi

