#!/bin/bash
##set -euo pipefail

# Color definitions
    green='\033[32m'
    red='\033[31m'
    cyan='\033[36m'
    purple='\033[35m'
    orange='\033[33m'
    nocolor='\033[m'

# Function to check file existence
    check_file() {
        if [[ ! -f "$1" ]]; then
            echo "Error: File $1 not found" >&2
            exit 1
            else echo "File $1 found"
        fi
    }

# Check required files
    check_file "pairwise_analysis.list.csv"
    check_file "selected_main-lineage_split.list"
    check_file "selected_sub-lineage_split.list"
    check_file "run_sample_ids.txt"

# Get the main lineages
    while IFS= read -r filt_lineage; do
        grep -E "${filt_lineage}" pairwise_analysis.list.csv | while IFS=';' read -r sampleID main_lineage sub_lineage; do
            echo "${filt_lineage},${sampleID}" >> tmp.lineage_samples_tuple.main.csv
        done
    done < selected_main-lineage_split.list

# Get the sub-lineages
    while IFS= read -r filt_lineage; do
        grep -E "${filt_lineage}" pairwise_analysis.list.csv | while IFS=';' read -r sampleID main_lineage sub_lineage; do
            echo "${filt_lineage},${sampleID}" >> tmp.lineage_samples_tuple.sub.csv
        done
    done < selected_sub-lineage_split.list

# Merge the files
    cut -d ',' -f2 tmp.lineage_samples_tuple.sub.csv > tmp.lineage_samples_tuple.sub.sampleID
    grep -v -f tmp.lineage_samples_tuple.sub.sampleID tmp.lineage_samples_tuple.main.csv > tmp.lineage_samples_tuple.main-final.csv
    cat tmp.lineage_samples_tuple.main-final.csv tmp.lineage_samples_tuple.sub.csv > tmp.lineage_samples_tuple.1.csv

# Create tuple for analysis
    grep -f run_sample_ids.txt tmp.lineage_samples_tuple.1.csv | cut -d ',' -f1 | sort | uniq > tmp.lineages.to.keep

# Remove lineages with less than 3 genomes
    while IFS=';' read -r filt_lineage; do
        count=$(grep -c "${filt_lineage}" tmp.lineage_samples_tuple.1.csv)
        if [[ ${count} -ge 3 ]]; then
            echo -e "${cyan}${filt_lineage}${green}: > 3 genomes, can be clustered (count: ${cyan}${count}${green})${nocolor}"
        elif [[ ${count} -lt 3 ]]; then
            echo "${filt_lineage}" >> tmp.lineages.to.remove.from.keep.list
            echo -e "${purple}${filt_lineage}${red}: < 3 genomes, cannot be clustered (count: ${purple}${count}${red})${nocolor}"
        fi
    done < tmp.lineages.to.keep

# In case any of the filetr files are empty
    touch tmp.lineages.to.remove.from.keep.list # incase this is empty
    touch tmp.lineages.to.keep.1
    touch tmp.lineages.to.keep

# Process the generated tuples
    grep -v -f tmp.lineages.to.remove.from.keep.list tmp.lineages.to.keep > tmp.lineages.to.keep.1
    grep -f tmp.lineages.to.keep.1 tmp.lineage_samples_tuple.1.csv > final.lineage_samples_tuple.csv
    grep -v -f tmp.lineages.to.keep tmp.lineage_samples_tuple.1.csv > final.skipped-lineages_tuple.csv