process COMPILE_SEQUENCING_STATS {

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/", mode: 'copy'

    input:
        val(runID)
        path(tbdb_results)
        path(who_results)
        path(mtbseq_compiled_strains)
        path(mtbseq_compiled_map_stats)
        path(lineage_fractions)
        val(sample_ids)

    output:
        path("${runID}.sequencing_summary.csv")
        path("sequencing_summary.csv"),         emit: analysis_summary
        path("who_resistance_summary.csv"),     emit: who_resistance
        path("tbdb_resistance_summary.csv"),    emit: tbdb_resistance
        path("lineage_samples_tuple.csv"),           emit: lineage_sample_tuple
        path("skipped-lineages.csv"),                emit: skipped_lineages

    script:

    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Create the lienage fraction strings
        Rscript ${params.r_script_dir}/tbprofiler_lineage_fractions.R \\
                        --tbprofiler    tbdb-tbprofiler.txt \\
                        --lineages      lineages.fractions.txt

    # Convert the list of sample IDs to a format suitable for grep
        echo '${sample_ids.join("\n")}' > run_sample_ids.txt

    # Generate summary statistics and create the sampleID,lineage df for
    ## creating into a channel 
        Rscript ${params.r_script_dir}/compile-sequencing-statistics.R \\
                    --minimum_coverage ${params.mtbseq_min_cov} \\
                    --runID ${runID} \\
                    --dictionary_path ${params.r_script_dir} \\
                    1>>.command.out \\
                    2>>.command.err || true # NOTE This is a hack to overcome the exit status 1
                    

    # Seperate out the genomes from this run into their own results file
        grep -f run_sample_ids.txt ${runID}.sequencing_summary.csv > tmp.${runID}.sequencing_summary.csv
        mv tmp.${runID}.sequencing_summary.csv ${runID}.sequencing_summary.csv

    # extract the lineages from the params.config file
        echo '${params.lineage_pairwise_main.join('\n')}' > selected_lineage_split.main
        
        echo '${params.lineage_pairwise_sub.join('\n')}' > selected_lineage_split.sub

        while read -r main_lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv

            grep -E "\${main_lineage}" pairwise_analysis.list.csv | while IFS=';' read -r sampleID main_lineage sub_lineage; do
                # Check if sub_lineage starts with main_lineage
                
                if [[ "\${sub_lineage}" == "\${main_lineage}"* ]]; then

                    # Append the result to the output file
                    echo "\${main_lineage},\${sampleID}" >> lineage_samples_tuple.csv   

                fi
            done
        done < selected_lineage_split.list

    # Remove that lineage if there are less than 3 genomes (minimum needewd for MTBSeq pairwise analysis)
        awk -F',' '
            {
                count[\$1]++      # Count occurrences of each lineage in column 1
                lines[NR] = \$0   # Store the entire line
                lineage[NR] = \$1 # Store the lineage (column 1)
            }
            END {
                for (i = 1; i <= NR; i++) {
                    if (count[lineage[i]] >= 3) {
                        print lines[i] # Print lines for lineages with >= 3 entries
                    }
                }
            }' "lineage_samples_tuple.csv" > "lineage_samples_tuple.csv.tmp"

        cp lineage_samples_tuple.csv lineage_samples_tuple.unfiltered.csv
        cp lineage_samples_tuple.csv.tmp lineage_samples_tuple.csv

        touch skipped-lineages.csv

    """
}
