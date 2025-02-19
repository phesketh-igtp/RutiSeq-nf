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
        path("archive/${runID}.sequencing_summary.csv")
        path("main/sequencing_summary.csv"),         emit: analysis_summary
        path("main/who_resistance_summary.csv"),     emit: who_resistance
        path("main/tbdb_resistance_summary.csv"),    emit: tbdb_resistance
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
        ##echo '${params.lineage_pairwise_exceptions.join('\n')}' > selected_lineage_exceptions.list

        while read -r lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv
            grep -E "\${lineage}" pairwise_analysis.list.csv | while IFS=';' read -r sampleID sub_lineage; do
                # Check if sub_lineage contains the lineage
                if [[ "\${sub_lineage}" == *"\${lineage}"* ]]; then
                    # Append the result to the output file
                    echo "\${lineage},\${sampleID}" >> lineage_samples_tuple.csv
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

        mv lineage_samples_tuple.csv lineage_samples_tuple.unfiltered.csv
        mv lineage_samples_tuple.csv.tmp lineage_samples_tuple.csv

    # Filter the lineages to contain ONLY the lineages from the newest run
        mv lineage_samples_tuple.csv all-lineage_samples_tuple.csv

        # Collect the lineage IDs as define from the tuple.csv
            grep -f run_sample_ids.txt all-lineage_samples_tuple.csv \\
                | cut -d ',' -f1 \\
                | sort \\
                | uniq \\
                | sed "s/\$/,/g" \\
                > run_sample_ids_taxonomy.txt

        # Create tuple for pairwise analysis
            grep -f run_sample_ids_taxonomy.txt \\
                all-lineage_samples_tuple.csv \\
                > lineage_samples_tuple.csv
        # Create tuple for skipped analysis
            grep -v -f run_sample_ids_taxonomy.txt \\
                all-lineage_samples_tuple.csv \\
                > skipped-lineages.csv

    # Move the outputs into folders
        mkdir -p main archive
        mv sequencing_summary.csv main/
        mv tbdb_resistance_summary.csv main/
        mv who_resistance_summary.csv main/
        mv ${runID}.sequencing_summary.csv archive/

    """
}