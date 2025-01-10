process COMPILE_SEQUENCING_STATS {

    tag "${runID}"

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
        path("pairwise_analysis_lineage_split.list"),   emit: pairwise_list
        path("lineage_samples_paths.csv"),              emit: lineage_sample_path

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
                --minimum_coverage      ${params.mtbseq.min_cov} \\
                --runID                 ${runID} \\
                --dictionary_path       ${params.r_script_dir} \\
                ${additional_args} \\
                

    # Create sample list of all the MTB lineages to be analyzed
    if [[ -f lineage_samples_paths.csv ]]; then rm lineage_samples_paths.csv; fi
    
    for lineage in ${params.lineage_pairwise.join(' ')}; do
        awk -F',' '\$2 ~ /'\$lineage'/ {print \$1}' pairwise_analysis_lineage_split.list |\\
                sed 's@_@\t@g' \\
                > \$lineage.samples.txt
        echo "\$lineage,\$lineage.samples.txt" >> lineage_samples_paths.csv
    done
    """
}