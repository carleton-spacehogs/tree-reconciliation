#!/bin/bash

# Loop through all .treefile files in the current directory
for tree_file in *.treefile; do
    # Define the output file name (replace .treefile extension with .nexus)
    output_file="${tree_file%.treefile}.nexus"

    # Call the Python script for conversion
    python convert_to_nexus.py "$tree_file" "$output_file"
done
