process COMPILE_CN_READS_SUMMARY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-03
    @version: 0.1
    @description: 
        This module runs an R script that combines all the results into the 
        format that can then be compiled into the final XLSX sheet produced
        to summarise the entire genome collection in then SUMMARY_WF()
*/

    tag "${runID}"

    conda params.r_stats_env

    publishDir "${params.outdir}/negative-controls/", mode: 'copy'

    input:
        val(runID)
        path(tbprofile_compiled)

        path(all_mtbseq_class)
        path(all_mtbseq_stats)
        
        path(k2_combined)
        path(stats_combined)

    output:
        path("${runID}_combined_qc_results.csv")
        path("combined_qc_results.csv"), emit: combined_cn_results

    script:
        """
        #Rscript ${params.r_script_dir}/negative-control-compile.R

        exit 1
        """

}