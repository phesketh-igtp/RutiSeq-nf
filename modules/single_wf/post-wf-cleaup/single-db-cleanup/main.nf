process POST_SINGLE_DB_CLEANUP {

    tag "${sampleID}"

    array 100

    input:
        val(sampleID)

    script:

        """
        # Check and remove files only if they exist
        rm -f ${params.outDir}/db/samples/${sampleID}/*_R*fastq.gz
        ## TB-Profiler cleanup
        rm -f ${params.outDir}/db/samples/${sampleID}/tbprofiler/*_R*fastq.gz
        rm -rf ${params.outDir}/db/samples/${sampleID}/tbprofiler/bam/
        rm -rf ${params.outDir}/db/samples/${sampleID}/tbprofiler/vcf/

        ## MTBSeq cleanup
        rm -f ${params.outDir}/db/samples/${sampleID}/mtbseq/*_R*fastq.gz
        rm -f ${params.outDir}/db/samples/${sampleID}/mtbseq/tbdb-${sampleID}.results.txt
        rm -f ${params.outDir}/db/samples/${sampleID}/mtbseq/who-${sampleID}.results.txt

        ## Snippy cleanup
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/*_R*fastq.gz
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}*gatk_position*
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.Mapping_and_Variant_Statistics.tab
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/${sampleID}.Strain_Classification.tab
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/tbdb-*
        rm -f ${params.outDir}/db/samples/${sampleID}/snippy/who-*
        """

}

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