process COMBINE_QC_RESULTS {

    tag "${runID}"

    publishDir "${params.outdir}/combined_qc", mode: 'copy'

    input:
    path qc_files
    val runID

    output:
    path "${runID}_combined_qc_results.csv", emit: combined_qc

    script:
    """
        set -e
        set -x

        echo -e "RunID,SampleID,orig.R1_reads,orig.R1_aveQ,orig.R2_reads,orig.R1_aveQ,orig.MTB_perc,filt.R1_reads,filt.R1_aveQ,filt.R2_reads,filt.R2_aveQ" > ${runID}_combined_qc_results.csv

        echo "Combining QC results"
        awk -v OFS=',' -v runid="${runID}" '{print runid,\$0}' *.qc.out >> ${runID}_combined_qc_results.csv

        sort ${runID}_combined_qc_results.csv | uniq > tmp
        mv tmp ${runID}_combined_qc_results.csv

    """
}