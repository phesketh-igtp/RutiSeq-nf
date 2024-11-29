process COMBINE_QC_RESULTS {

    publishDir "${params.outdir}/combined_qc", mode: 'copy'

    input:
        path qc_files
        val runID

    output:
    path "${runID}_${timestamp}_combined_qc_results.csv", emit 

    script:
    timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())

    """

    echo -e "RunID\tSampleID\torig.R1_reads\torig.R1_aveQ\torig.R2_reads\torig.R1_aveQ\torig.MTB_perc\tfilt.R1_reads\tfilt.R1_aveQ\tfilt.R2_reads\tfilt.R2_aveQ" > ${runID}_${timestamp}_combined_qc_results.csv
    awk -v OFS=',' -v runid="${runID}" '{print runid,\$0}' ${qc_files} >> ${runID}_${timestamp}_combined_qc_results.csv

    """
}