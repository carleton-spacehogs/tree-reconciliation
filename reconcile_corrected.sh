#!/bin/bash

source ./scripts/utils.sh

Help()
{
	# Display Help
	echo
	echo "Jimmy's effort to automate tree conciliation with ecceTERA
working under Professor Rika Anderson in Carleton College, MN. Started March 2023

choose only one of the entry points to start this program

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

exit_msg_errorCode "You have to gimme an option to work with. Need help???" 0
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

exit_msg_errorCode () {
	conda deactivate
	echo $1
	echo exit $2
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
		exit_msg_errorCode "ecceTERA failed... Stop here" 1
	fi

	analysis $gene_name # analysis is in scripts/utils.sh
}


makeTree_and_reconcile() {
	pre_trim=$1

	echo trimming the original alignment $pre_trim
	# trim out positions with mostly gaps
	trimal -in $pre_trim -out $trimv1 -gt 0.15 # -automated1
	
	# we excluded genes with more than 20% gaps in the trimmed alignment
	python3 scripts/keep_seq_geq_Xpercent.py $trimv1 $trimv2 80
	rm $trimv1

	validate_alignment $trimv2 $num_seq_min $alignment_len_min 0 # 0: exit when there is a problem

	if [ "$stop_before_tree_building" = true ]; then
		exit_msg_errorCode "I see --stop_before_tree_building, so I stop" 0
	fi

	which iqtree
	if [ $? -ne 0 ]; then
		exit_msg_errorCode "can not use iqtree. Do \"which iqtree\" and see the problem." 1
	fi

	iqtree -s $trimv2 --nmax 3000 -m MFP -bb 1000 -wbtl -ntmax 20 -nt AUTO -T $num_core --prefix $outfile_prefix -madd "C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,UL2,UL3,EX_EHO,LG4M,LG4X,CF4,LG+C10,LG+C20,LG+C30,LG+C40,LG+C50,LG+C60" -mrate $rate_models "E,I,G,I+G,R" # List of 1. additional models 2. rate model for heterogeneity among sites

	if [ "$stop_before_reconciliation" = true ]; then
		exit_msg_errorCode "I see --stop_before_reconciliation, so I stop" 0
	fi

	reconcile_and_analysis $iqtree_ufboot
}


align_makeTree_and_reconcile() {
	# if [ ! -f $trimv2 ]; then
	# 	echo aligning this gene: $gene_name, and alignment output: $pre_trim
	# 	# muscle executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
	# 	./muscle5.1 -align $gene_seq_file -output $pre_trim
	# fi

	if [ ! -f $trimv2 ]; then
		echo aligning this gene: $gene_name, and alignment output: $pre_trim
		# muscle executable downloaded from https://github.com/rcedgar/muscle/releases/tag/5.1.0
		# Count the number of lines in the gene sequence file
		num_genes=$(wc -l < "$gene_seq_file")
		# Check if the number of genes is greater than 1000
		if [ "$num_genes" -gt 1000 ]; then
		echo doingsuper5
			./muscle5.1 -super5 $gene_seq_file -output $pre_trim -threads 20
		else
			./muscle5.1 -align $gene_seq_file -output $pre_trim -threads 20
		fi
	fi
	
	if [ $? -ne 0 ]; then
		exit_msg_errorCode "alignment of $gene_seq_file not successful, stop here" 1
	fi

	makeTree_and_reconcile $pre_trim
}


extractSeq_align_makeTree_and_reconcile() {
	COG=$1
	# sanity check
	echo $COG | grep COG
	
	if [[ $? != 0 ]]; then
		
		exit_msg_errorCode "The COG number $COG must be in the format of \"COGXXXX\"" 1
	fi

	if test -f $iqtree_ufboot; then
		#if [ "$stop_before_reconciliation" = false ]; then
		echo I found the gene_tree: skipping to the reconciliation.
		reconcile_and_analysis $iqtree_ufboot
		#fi
	else
		echo $all_seq_fasta
		python3 scripts/grep_seq_given_COG.py $COG
		if [ $? -ne 0 ]; then
			exit_msg_errorCode "grabbing the sequence for $COG not successful, stop here" 1
		fi
		align_makeTree_and_reconcile
	fi
}

am_I_done() {
	if [ "$overwrite" = true ]; then
		echo I am over writing stuff
	elif [ -f "$R_plot" ]; then
		python3 scripts/summarize_reconciliation.py $gene_name 0 rewrite
		exit_msg_errorCode "$R_plot exists. I think I am done" 0
	fi
	activate_conda_env
}

source scripts/declare_file_location.sh --clock_model $clock_model
validate_required_folders

# if given both gene and species tree, just do the reconciliation
if [ ! -z "$gen_graph_only" ]; then
	source scripts/declare_file_location.sh --gene_name $gen_graph_only --clock_model $clock_model
	activate_conda_env
	analysis
	exit_msg_errorCode "Done" $?
elif [ ! -z "$gene_tree" ]; then
	parseGeneName_and_declareFilenames $gene_tree
	am_I_done
	reconcile_and_analysis $gene_tree
	exit_msg_errorCode "Done" $?
elif [ ! -z "$alignment" ]; then
	parseGeneName_and_declareFilenames $alignment
	am_I_done
	makeTree_and_reconcile $alignment
	exit_msg_errorCode "Done" $?
elif [ ! -z "$gene_sequences" ]; then
	parseGeneName_and_declareFilenames $gene_sequences
	am_I_done
	align_makeTree_and_reconcile
	exit_msg_errorCode "Done" $?
elif [ ! -z "$COG" ]; then
	source ./scripts/declare_file_location.sh --clock_model $clock_model --gene_name $COG
	am_I_done
	extractSeq_align_makeTree_and_reconcile $COG
	exit_msg_errorCode "Done-COG" $?
else
	Help
fi
