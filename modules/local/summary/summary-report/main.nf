process GENERATE_SUMMARY_REPORT {

    publishDir "${params.outdir}/bbdd/results/", mode: 'copy'

    tag "${runID}"

    input:
        val runID
        path pairwise_clusters
        path pairwise_matrix
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
                --clusters ${pairwise_clusters} \\
                --matrices ${pairwise_matrix} \\
                --output ${runID}_RutiSeq-results.xlsx
        """
}