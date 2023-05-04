#!/bin/bash

source scripts/declare_file_location.sh

COG_graphs=$(ls R-plots/histogram)

for graph in $COG_graphs; do
	COG=$(echo $graph | awk -F"-" '{print $3}')
	grep $COG $COG_summary > /dev/null
	if [ $? -ne 0 ]; then
		echo $COG has the histogram:
		ls -lh R-plots/histogram/$graph
		echo but $COG not in $COG_summary, so I re-do the analysis
		./reconcile.sh --gen_graph_only $COG
	fi
done