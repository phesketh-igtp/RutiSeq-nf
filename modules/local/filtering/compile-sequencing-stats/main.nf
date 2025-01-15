process COMPILE_SEQUENCING_STATS {

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/tbprofiler/pairwise/", mode: 'copy'

    input:
        val(runID)
        path tbdb_results
        path who_results
        path mtbseq_compiled_strains
        path mtbseq_compiled_map_stats
        path lineage_fractions

    output:
        path("${runID}.sequencing_summary.csv")
        path("sequencing_summary.csv"),                 emit: analysis_summary
        path("who_resistance_summary.csv"),             emit: who_resistance
        path("tbdb_resistance_summary.csv"),            emit: tbdb_resistance
        path("lineage_samples_tuple.csv"),              emit: lineage_sample_tuple

    script:
    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Create the lienage fraction strings
    Rscript ${params.r_script_dir}/tbprofiler_lineage_fractions.R \\
                    --tbprofiler    tbdb-tbprofiler.txt \\
                    --lineages      lineages.fractions.txt


    # Generate summary statistics and create the sampleID,lineage df for
    ## creating into a channel 
    Rscript ${params.r_script_dir}/compile-sequencing-statistics.R \\
                --mtbseq_statistics     "Mapping_and_Variant_Statistics".tab \\
                --mtbseq_classification "Strain_Classification".tab \\
                --tbprofiler_tbdb       "tbdb-tbprofiler.txt" \\
                --tbprofiler_who        "who-tbprofiler.txt" \\
                --lineage_fractions     "tbprofiler.lineages.fractions.txt" \\
                --minimum_coverage      ${params.mtbseq_min_cov} \\
                --runID                 ${runID} \\
                --dictionary_path       ${params.r_script_dir} \\
                ${additional_args}
                
    # extract the lineages from the params.config file
    echo '${params.lineage_pairwise.join('\n')}' > selected_lineage_split.list

    while read -r lineage; do
        # Use grep to find matching lines from pairwise_analysis.list.csv
        grep "\${lineage}" pairwise_analysis.list.csv | while IFS=';' read -r sampleID sub_lineage; do
            # Check if sub_lineage contains the lineage
            if [[ "\${sub_lineage}" == *"\${lineage}"* ]]; then
                # Append the result to the output file
                echo "\${lineage},\${sampleID}" >> lineage_samples_tuple.csv
            fi
        done
    done < selected_lineage_split.list

    # Remove that lineage is there are less than 3 genomes (minimum needewd for MTBSeq pairwise analysis)
        
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

    mv lineage_samples_tuple.csv lineage_samples_tuple.unfiltered.csv

    mv lineage_samples_tuple.csv.tmp lineage_samples_tuple.csv

    cat lineage_samples_tuple.csv

    """
}