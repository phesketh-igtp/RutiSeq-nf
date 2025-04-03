process COMPILE_CN_READS_SUMMARY {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/negative-controls/", mode: 'copy'

    input:
        val(runID)
        path(all_mtbseq_class)
        path(all_mtbseq_stats)
        path(tbprofile_compiled)
        path(k2_combined)
        path(stats_combined)

    output:
        path "${runID}_combined_qc_results.csv"
        path "combined_qc_results.csv", emit: combined_cn_results

    script:
        """



        """
}