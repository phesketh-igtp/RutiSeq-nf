process COMBINE_QC_RESULTS {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/negative-controls/", mode: 'copy'

    input:
        val(runID)
        path(qc_files)
        path(kaiju_files)


    output:
    path "${runID}_combined_qc_results.csv", emit: combined_qc

    script:
        """
        set -e
        set -x

        #Creating combined QC results file
            echo -e "RunID\tSampleID\torig.R1_reads\torig.R1_aveQ\torig.R2_reads\torig.R1_aveQ\torig.MTB_perc\tfilt.R1_reads\tfilt.R1_aveQ\tfilt.R2_reads\tfilt.R2_aveQ" > ${runID}_combined_nc_qc_results.negative-control.csv
            awk -v OFS='\t,' -v runid="${runID}" '{print runid\$0}' ${qc_files} >> ${runID}_combined_nc_qc_results.negative-control.csv

        # Creating Kaiju results
            echo -e "RunID\tSampleID\tpercent	reads\ttaxon_id\ttaxon_name" > ${runID}_combined_kaiju_results.negative-control.tsv
            awk -v OFS='\t,' -v runid="${runID}" '{print runid\$0}' ${kaiju_files} >> ${runID}_combined_kaiju_results.negative-control.tsv

        """
}