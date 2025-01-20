#!/bin/bash

options=$(getopt -o if:c:cm: --long help,flags:,COG_list:,clock_model: -- "$@")
eval set -- "$options"

while true; do
	case "$1" in
		-h|--help)
			echo Help
			exit
			;;
		-f|--flags)
			flags="$2";
			shift 2
			;;
		-c|--COG_list)
			COG_list="$2";
			shift 2
			;;
		-cm|--clock_model)
			CM="$2"
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

dispatch_example="./reconcile_TN.sh --clock_model ugam1 --flags \" --COGs\" --COG_list \"COG2046 COG2059\""

is_memory_ok() {
	free_memory_lower_bound=$1
	percent_free_memory=$(free | grep Mem | awk '{print $4/$2 * 100.0}')
	# ahh, string comparison in shell script is a pain...
	if (( $(bc <<<"$free_memory_lower_bound > $percent_free_memory") )); then
		echo false
	else
		echo true
	fi
}

job_running() {
	PID=$1
	ps $PID | grep -E "^T" # T for stop jobs
	if [ $? -eq 0 ]; then
		echo false
	else
		echo true
	fi
}

acknowledge_execute() {
	flags=$1
	clockModel=$2
	source ./scripts/declare_file_location.sh --clock_model $clockModel --gene_name $COG
	if [[ "$flags" == *"--gene_tree"* ]]; then
		input=$iqtree_ufboot
	elif [[ "$flags" == *"--alignment"* ]]; then
		input=$alignment
	elif [[ "$flags" == *"--gene_sequence"* ]]; then
		input=$gene_seq_file
	elif [ -z "$flags" ]; then
		echo "you must give me --flags, for example:"
		echo $dispatch_example
	else # for "--COG" and "--gen_graph_only"
		input=$COG
	fi

	command="./reconcile_corrected.sh --clock_model $clockModel $flags $input"
	echo I am doing this command:
	echo " "
    echo "  _____  ____   _____           "
    echo " / ____ / __ \ / ____|          "
    echo "| |    | |  | | |  __           "
    echo "| |    | |  | | | |_ |          "
    echo "| |____| |__| | |__| |          "
    echo " \_____ \____/ \_____|          "
    echo "                                "
		echo $command
	eval $command # >> run_${COG}_job.log
}

if [ -z $(echo $COG_list | sed 's/ //g') ]; then
	echo gimme a list of COGs, for example:
	echo $dispatch_example
	exit 1
fi

for COG in $COG_list; do
	echo "curCOG is $COG"
	if [ "$CM" = "all" ] && [[ "$flags" != *"--stop_before"* ]] ; then
		echo You told me to run all clock models, so I am doing it:
		acknowledge_execute $flags ugam1
		acknowledge_execute $flags cir1 & # should just be doing ecceTERA
		acknowledge_execute $flags ln3 & # so it we can do it concurrently
	elif [ "$CM" = "ugam1 cir1" ] && [[ "$flags" != *"--stop_before"* ]] ; then
		echo You told me to run ugam1 cir1, so I am doing it:
		acknowledge_execute $flags ugam1
		acknowledge_execute $flags cir1 & # should just be doing ecceTERA
	else
		acknowledge_execute $flags $CM
	fi
done

echo "done"
