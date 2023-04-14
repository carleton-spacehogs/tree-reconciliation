#!/bin/bash

free_memory_start_point=50
free_memory_stop_point=25

ignore_memory=false

options=$(getopt -o ih --long ignore_memory,help -- "$@")
eval set -- "$options"

while true; do
	case "$1" in
		-h|--help)
			echo Help
			exit
			;;
		-i|--ignore_memory)
			ignore_memory=true;
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


max_jobs=3
# COG_list="COG2146 COG1785 COG0783 COG0221 COG1668 COG2181 COG0239 COG2132 COG0477 COG1151" # tree1
# COG_list="COG0248 COG0003 COG0025 COG0053 COG0155 COG0168 COG0288 COG0306 COG0310 COG0370" # tree2
# COG_list="COG0376 COG0428 COG0475 COG0529 COG0530 COG0555 COG0573 COG0581 COG0598 COG0605" # tree3
# COG_list="COG0659 COG0704 COG0725 COG0753 COG0798 COG0803 COG0855 COG0861 COG1006 COG1055" # tree4
# COG_list="COG1108 COG1121 COG1178 COG1218 COG1226 COG1230 COG1283 COG1320 COG1324 COG1393" # inorganic1
COG_list="COG1496 COG1528 COG1553 COG1633 COG1814 COG1840 COG1863 COG1910 COG1914 COG1918" # inorganic2
# COG_list="COG1971 COG2046 COG2059 COG2072 COG2111 COG2128 COG2193 COG2212 COG2221 COG2223" # \\ waiting for inorganic3
# COG_list="COG2239 COG2346 COG2382 COG2710 COG2895 COG2897 COG2906 COG2998 COG3221 COG3420" # inorganic4
# COG_list="COG3540 COG3639 COG3696 COG3746 COG4149 COG4651" # deming work1

num_done=0
num_COGs=$(echo $COG_list | awk -F 'COG' '{print NF - 1}')

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

if [ $ignore_memory = true ]; then
	for COG in $COG_list; do
		./reconcile.sh --COG $COG --stop_before_reconciliation
	done
else
	for COG in $COG_list; do
		./reconcile.sh --COG $COG &
		# ./reconcile.sh --gene_tree iqtree_gene_trees/${COG}.ufboot &
		COG_PID=$!
		job_str="$COG_PID for $COG"
		is_done=false

		echo Dispatching job $job_str
		while [[ $is_done = false ]]; do
			# echo I am waiting, check in 100 seconds >> tmp/tmp.txt
			sleep 3
			ps $COG_PID &>/dev/null
			if [ $? = 1 ]; then is_done=true; fi

			if [ $is_done = false ]; then
				if [ $(is_memory_ok $free_memory_start_point) = true ]; then
					if [ $(job_running $COG_PID) = false ]; then 
						echo $job_str was paused, now continue >> tmp/tmp2.txt
						kill -CONT $COG_PID
					else
						echo $job_str memory ok and job running >> tmp/tmp2.txt
					fi
				elif [ $(is_memory_ok $free_memory_stop_point) = false ]; then # memory not ok
					if [ $(job_running $COG_PID) = true ]; then
						echo $job_str is running, now paused due to memory >> tmp/tmp2.txt
						kill -STOP $COG_PID
					else
						echo $job_str memory not ok and job not running >> tmp/tmp2.txt
					fi
				else
					echo $job_str memory between start and stop point, not doing anything >> tmp/tmp2.txt
				fi
			else
				echo $job_str is done
			fi
		done
	done
fi

# COG_list=$(grep -E ",40$" diamond_COG_count.txt | awk -F ',' '{print $1}' | head -n 5)
# COG_list="COG5677 COG5626 COG4291 COG1421 COG0269"

# 20: COG5734 COG5674 COG5768 COG5778 COG5790
# 40: COG5677 COG5626 COG4291 COG1421 COG0269
# 60: COG3265 COG4691 COG4619 COG3010 COG1439
# 80: COG5716 COG1448 COG1631 COG1889 COG3155
# 160: COG4587 COG1125 COG1464 COG0692 COG2850
# 320: COG2022 COG2453 COG5659 COG1966

