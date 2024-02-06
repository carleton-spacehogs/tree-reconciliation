#!/bin/bash

predictions_csv="Predictions_tuned.csv"

# Loop through all genome CSV files
for genome_csv in *.csv; do
    if [ "$genome_csv" != "$predictions_csv" ]; then
        # Define the output CSV file name
        output_csv="${genome_csv%.csv}_merged.csv"

        # Call the Python script for merging
        python merge_csv.py "$genome_csv" "$predictions_csv" "$output_csv"
    fi
done
