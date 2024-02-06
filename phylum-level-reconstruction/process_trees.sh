#!/bin/bash

# Loop through all tree files
for tree_file in *rooted.treefile; do
    # Define the output CSV file name
    output_csv="${tree_file%.treefile}.csv"

    # Call the Python script
    python extract_genomes.py "$tree_file" "$output_csv"
done

