#!/bin/bash

free_memory_start_point=50
free_memory_stop_point=75
max_jobs=3
COG_list=$(grep -E ",320$" diamond_COG_count.txt | awk -F ',' '{print $1}' | head -n 5)

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

for COG in $COG_list; do
	./reconcile.sh --COG $COG &
	COG_PID=$!
	job_str="$COG_PID for $COG"
	is_done=false

	echo Dispatching job $job_str

	while [ $is_done = false ]; do
		echo I am waiting, check in 100 seconds >> tmp/tmp.txt
		sleep 100
		ps $COG_PID &>/dev/null
		if [ $? = 1 ]; then is_done=true; fi

		if [ $is_done = false ]; then
			if [ $(is_memory_ok $free_memory_start_point) = true ]; then
				if [ $(job_running $COG_PID) = false ]; then 
					echo $job_str was paused, now continue
					kill -CONT $COG_PID
				else
					echo memory ok and job running
				fi
			elif [ $(is_memory_ok $free_memory_stop_point) = false ]; then # memory not ok
				if [ $(job_running $COG_PID) = true ]; then
					echo $job_str is running, now paused due to memory
					kill -STOP $COG_PID
				else
					echo memory not ok and job not running
				fi
			else
				echo memory between start and stop point, not doing anything
			fi
		else
			echo $job_str is done
		fi
	done

done



# grep -n ,20 diamond_COG_count.txt
# 640:COG5734,20 -> done
# 641:COG5674,20 -> done
# 642:COG5768,20 -> done
# 643:COG5778,20 -> done
# 644:COG5790,20 -> done

# COG_list=$(grep -E ",40$" diamond_COG_count.txt | awk -F ',' '{print $1}' | head -n 5)
# COG_list="COG5677 COG5626 COG4291 COG1421 COG0269"

# 40: COG5677 COG5626 COG4291 COG1421 COG0269
# 60: COG3265 COG4691 COG4619 COG3010 COG1439
# 80: COG5716 COG1448 COG1631 COG1889 COG3155
# 160: COG4587 COG1125 COG1464 COG0692 COG2850
# 320: COG2022 COG2453 COG5659 COG1966

