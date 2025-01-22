process GENERATE_SUMMARY_REPORT {

    publishDir "${params.outdir}/bbdd/results/${runID}", mode: 'copy'

    input:
        val runID
        path pairwise_clusters_processed
        path analysis_summary
        path who_resistance
        path tbdb_resistance

    output:
        path("${runID}_RutiSeq-results.xlsx")

    script:

        """
        Rscript ${params.r_script_dir}/generate_summary_report.R \\
                --summary ${analysis_summary} \\
                --who_res ${who_resistance} \\
                --tbdb_res ${tbdb_resistance} \\
                --clusters ${pairwise_clusters_processed} \\
                --output ${runID}_RutiSeq-results.xlsx
        """
}