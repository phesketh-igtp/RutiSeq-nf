process GENERATE_SUMMARY_REPORT {

    conda params.r_stats_env

    publishDir "${params.outdir}/results/excel-out", mode: 'copy'

    input:
        val runID
        path pairwise_clusters_processed
        path sequencing_summary
        path who_resistance
        path tbdb_resistance

    output:
        path("${runID}_RutiSeq-results.xlsx")

    script:

        """
        Rscript ${params.r_script_dir}/generate_summary_report.R \\
                --summary ${sequencing_summary} \\
                --who_res ${who_resistance} \\
                --tbdb_res ${tbdb_resistance} \\
                --clusters ${pairwise_clusters_processed} \\
                --rlibrary ${params.r_script_dir} \\
                --output ${runID}_RutiSeq-results.xlsx
        """
}