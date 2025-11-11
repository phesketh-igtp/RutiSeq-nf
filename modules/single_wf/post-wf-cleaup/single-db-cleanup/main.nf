process POST_SINGLE_DB_CLEANUP {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: DEPRECIATED - CURRENTLY NOT IN USE!
        This module clean up all the intermediate files that are carried over through the process
        of emits/publications - this has been an awkward workaround to ensuring that one 
        genome that do not have the necessary outputs for the PAIRWISE_WF() are analysed.
        Unfortunately at each module step then the tuple containing all the paths to the previous
        are published in the output directory (mode: 'copy') - this final module just ensures that 
        all intermediate files are removed from the publish directory - even though they would have
        been remove during the process.
*/

    tag "${sampleID}"

    array 50

    input:
        tuple val(sampleID)

    script:

        """
        # Check and remove files only if they exist
            rm -f ${params.outDir}/db/read-qc/mtbc_reads/${sampleID}*.fastq.gz

            rm -f ${params.outDir}/db/tbprofiler/${sampleID}*.fastq.gz

            rm -f ${params.outDir}/db/tbprofiler/who-only/${sampleID}*.fastq.gz       
            rm -f ${params.outDir}/db/tbprofiler/who-only/${sampleID}/tbdb-${sampleID}.results.txt

            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/${sampleID}*.fastq.gz
            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/tbdb-${sampleID}.results.txt
            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/who-${sampleID}.results.txt

            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}*.fastq.gz
            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/SNP-Profiles/tbdb-${sampleID}*results.txt
            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/SNP-Profiles/who-${sampleID}*results.txt
            rm -f ${params.outDir}/db/mtbseq/samples/${sampleID}/SNP-Profiles/*.tab
        """

}