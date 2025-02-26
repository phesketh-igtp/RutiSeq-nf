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
        val(sampleID_list)

    output:
        path("${runID}.sequencing_summary.csv")

        path("sequencing_summary.csv"),      emit: analysis_summary
        path("who_resistance_summary.csv"),  emit: who_resistance
        path("tbdb_resistance_summary.csv"), emit: tbdb_resistance

        path("pairwise_analysis.list.csv"),  emit: pairwise_analysis_list

    script:

    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Create the lienage fraction strings
        Rscript ${params.r_script_dir}/tbprofiler_lineage_fractions.R \\
                        --tbprofiler    ${tbdb_results} \\
                        --lineages      lineages.fractions.txt

    # Generate summary statistics and create the sampleID,lineage df for
    ## creating into a channel 
        Rscript ${params.r_script_dir}/compile-sequencing-statistics.R \\
                    --minimum_coverage ${params.mtbseq_min_cov} \\
                    --runID ${runID} \\
                    --dictionary_path ${params.r_script_dir}

    # Convert the list of sample IDs to a format suitable for grep
        echo '${sampleID_list.join("\n")}' > run_sample_ids.txt

    # Seperate out the genomes from this run into their own results file
        grep -f run_sample_ids.txt ${runID}.sequencing_summary.csv > tmp.${runID}.sequencing_summary.csv
        mv tmp.${runID}.sequencing_summary.csv ${runID}.sequencing_summary.csv

    # Create the file to go to the tuple seperation
        awk -F "\t" '{ if ( \$14 > ${params.mtbseq_min_cov} ) print \$4 }' Mapping_and_Variant_Statistics.tab | sort | uniq > min.qual.genomes
        sed 's/\t/,/g' tbdb-tbprofiler.txt | cut -d ',' -f1,2,3 > tmp.pairwise_analysis.list.csv
        grep -f min.qual.genomes tmp.pairwise_analysis.list.csv > pairwise_analysis.list.csv

    """
}    