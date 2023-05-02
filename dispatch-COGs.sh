#!/bin/bash

free_memory_start_point=50
free_memory_stop_point=25

options=$(getopt -o if:c: --long help,flags:,COG_list: -- "$@")
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

if [ -z $(echo $COG_list | sed 's/ //g') ]; then
	echo gimme a list of COGs, for example:
	echo ./reconcile.sh --COGs \"COG2046 COG2059\"
	exit 1
fi

for COG in $COG_list; do
	source ./scripts/declare_file_location.sh $COG
	if [[ "$flags" == *"--gene_tree"* ]]; then
		input=$iqtree_ufboot
	elif [[ "$flags" == *"--alignment"* ]]; then
		input=$alignment
	elif [[ "$flags" == *"--gene_sequence"* ]]; then
		input=$gene_seq_file
	else # "--COG" and "--gen_graph_only"
		input=$COG
	fi

	command="./reconcile.sh $flags $input"
	echo I am doing this command:
		echo $command
	eval $command >> run_${COG}_job.log
done

# COG_list=$(grep -E ",40$" diamond_COG_count.txt | awk -F ',' '{print $1}' | head -n 5)
# COG_list="COG5677 COG5626 COG4291 COG1421 COG0269"

# 20: COG5734 COG5674 COG5768 COG5778 COG5790
# 40: COG5677 COG5626 COG4291 COG1421 COG0269
# 60: COG3265 COG4691 COG4619 COG3010 COG1439
# 80: COG5716 COG1448 COG1631 COG1889 COG3155
# 160: COG4587 COG1125 COG1464 COG0692 COG2850
# 320: COG2022 COG2453 COG5659 COG1966

