process GENERATE_SUMMARY_REPORT {

    tag "${runID}"

    input:
        val runID
        path pairwise_clusters,
        path pairwise_matrix,
        path analysis_summary,
        path who_resistance,
        path tbdb_resistance

    output:
        path("${YYMMDD}_${runID}_results.xlsx")

    shell

        """
        Rscript \\
                --summary ${analysis_summary} \\
                --who_res ${who_resistance} \\
                --tbdb_res ${tbdb_resistance} \\
                --clusters ${pairwise_clusters} \\
                --matrices ${pairwise_matrix} \\
                --output ${YYMMDD}_${runID}_results.xlsx

        """
}