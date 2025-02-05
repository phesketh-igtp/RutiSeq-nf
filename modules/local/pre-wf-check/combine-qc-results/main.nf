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

        echo "Current directory: \$(pwd)"
        echo "Listing input files:"
        #ls -l ${qc_files}

        echo "Creating combined QC results file"
        echo -e "RunID\tSampleID\torig.R1_reads\torig.R1_aveQ\torig.R2_reads\torig.R1_aveQ\torig.MTB_perc\tfilt.R1_reads\tfilt.R1_aveQ\tfilt.R2_reads\tfilt.R2_aveQ" > ${runID}_combined_qc_results.csv

        echo "Combining QC results"
        awk -v OFS=',' -v runid="${runID}" '{print runid,\$0}' *.qc.out >> ${runID}_combined_qc_results.csv

        echo "Contents of combined QC results file:"
        cat ${runID}_combined_qc_results.csv

        echo "Final directory contents:"
        ls -l
    """
}