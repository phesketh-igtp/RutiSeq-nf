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

    array 100

    input:
        val(sampleID)

    script:

        """
        # Check and remove files only if they exist
            rm -f ${params.outDir}/db/sample/${sampleID}/tbprofiler/*.fastq.gz
            rm -f ${params.outDir}/db/sample/${sampleID}/tbprofiler/bam

            rm -f ${params.outDir}/db/sample/${sampleID}/mtbseq/${sampleID}*.fastq.gz
            rm -f ${params.outDir}/db/sample/${sampleID}/mtbseq/tbdb-${sampleID}.results.txt
            rm -f ${params.outDir}/db/sample/${sampleID}/mtbseq/who-${sampleID}.results.txt

            rm -f ${params.outDir}/db/sample/${sampleID}/snippy/${sampleID}*.fastq.gz
            rm -f ${params.outDir}/db/sample/${sampleID}/snippy/tbdb-*.txt
            rm -f ${params.outDir}/db/sample/${sampleID}/snippy/who-*.txt
        """

}