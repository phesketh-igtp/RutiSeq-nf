process COMPILE_SEQUENCING_STATS3 {

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/", mode: 'copy'

    input:
        val(runID)
        path(tbdb_results)
        path(who_results)
        path(mtbseq_compiled_strains)
        path(mtbseq_compiled_map_stats)
        path(lineage_fractions)
        val(sampleID_list)

    output:
        path("${runID}.sequencing_summary.csv")

        path("sequencing_summary.csv"),             emit: analysis_summary
        path("who_resistance_summary.csv"),         emit: who_resistance
        path("tbdb_resistance_summary.csv"),        emit: tbdb_resistance

        path("final.lineage_samples_tuple.csv"),    emit: lineage_sample_tuple
        path("skipped-lineages.csv"),               emit: skipped_lineages

    script:

    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Create the lienage fraction strings
        Rscript ${params.r_script_dir}/tbprofiler_lineage_fractions.R \\
                        --tbprofiler    ${tbdb_results} \\
                        --lineages      lineages.fractions.txt

    # Convert the list of sample IDs to a format suitable for grep
        echo '${sampleID_list.join("\n")}' > run_sample_ids.txt

    # Generate summary statistics and create the sampleID,lineage df for
    ## creating into a channel 
        Rscript ${params.r_script_dir}/compile-sequencing-statistics.R \\
                    --minimum_coverage ${params.mtbseq_min_cov} \\
                    --runID ${runID} \\
                    --dictionary_path ${params.r_script_dir}

    #··································································································#

    # Seperate out the genomes from this run into their own results file
        grep -f run_sample_ids.txt ${runID}.sequencing_summary.csv > tmp.${runID}.sequencing_summary.csv
        mv tmp.${runID}.sequencing_summary.csv ${runID}.sequencing_summary.csv

    #··································································································#

    # extract the lineages from the params.config file

        echo '${params.lineage_pairwise_sub.join('\n')}' > selected_sub-lineage_split.list
        echo '${params.lineage_pairwise_main.join('\n')}' > selected_main-lineage_split.list

        # Get the main lineages 
        while read -r filt_lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv
            grep -E "\${filt_lineage}" pairwise_analysis.list.csv | \\
                                                while IFS=';' read -r sampleID main_lineage sub_lineage; do

                    
                if [[ "\${filt_lineage}" == "\${main_lineage}"* ]]; then

                    # Append the result to the output file
                    echo "\${filt_lineage},\${sampleID}" >> lineage_samples_tuple.main.csv   

                fi
            done
        done < selected_main-lineage_split.list

        # get the sub-lineages
        while read -r filt_lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv
            grep -E "\${filt_lineage}" pairwise_analysis.list.csv | \\
                                                while IFS=';' read -r sampleID main_lineage sub_lineage; do

                # Append the result to the output file
                echo "\${filt_lineage},\${sampleID}" >> lineage_samples_tuple.sub.csv   

            done
        done < selected_sub-lineage_split.list

    #··································································································#

        cut -d ',' -f2 lineage_samples_tuple.sub.csv > tmp.lineage_samples_tuple.sub.sampleID
        grep -v -f tmp.lineage_samples_tuple.sub.sampleID lineage_samples_tuple.main.csv > lineage_samples_tuple.main-final.csv

        cat lineage_samples_tuple.main-final.csv lineage_samples_tuple.sub.csv > final.lineage_samples_tuple.csv
    """
}    