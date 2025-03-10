#!/bin/bash

input_fasta="$1"  # Input FASTA file
output_file="$2"  # Output TSV file

# Process the FASTA file
grep -v "^>" "$input_fasta" | tr -d '\n' | awk '{
    for (i=1; i<=length($0); i++) 
        print i "\t" substr($0, i, 1)
}' > "$output_file"

echo "Conversion completed: Output saved to $output_file"