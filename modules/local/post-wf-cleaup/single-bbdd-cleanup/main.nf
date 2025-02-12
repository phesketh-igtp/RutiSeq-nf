process POST_SINGLE_BBDD_CLEANUP {

    /*
        This module clean up all the intermediate files that are carried over through the process
            of emits/publications - this has been an awkward workaround to ensuring that one 
            genome that do not have the necessary outputs for the PAIRWISE_WF() are analysed.
        Unfortunately at each module step then the tuple containing all the paths to the previous
            are published in the output directory (mode: 'copy') - this final module just ensures that 
            all intermediate files are removed from the publish directory - even though they would have
            been remove during the process.
    */

    tag "${sampleID}"

    array 100

    input:
        tuple val(sampleID)

    script:

        """
        # Check and remove files only if they exist
            rm -f ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R1.fastq.gz
            rm -f ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R2.fastq.gz

            rm -f  ${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R1.fastq.gz       
            rm -f  ${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R2.fastq.gz

            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R1.fastq.gz       
            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R2.fastq.gz
            rm -f  ${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}/tbdb-${sampleID}.results.txt

            rm -f  ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${sampleID}_mtbc_R1.fastq.gz
            rm -f  ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${sampleID}_mtbc_R2.fastq.gz
            rm -f  ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/tbdb-${sampleID}.results.txt
            rm -f  ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/who-${sampleID}.results.txt

            rm -f ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}*R1.fastq.gz
            rm -f ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}*R2.fastq.gz
            rm -f ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/tbdb-${sampleID}.results.txt
            rm -f ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/who-${sampleID}.results.txt
            rm -f ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/*.tab
        """

}