process MTBSEQ_STATS_COMPILE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process compiles the statistics and classification files from the MTBseq analysis.
        It concatenates the files from all samples and removes the header lines.
        The output files are saved in the specified output directory.
        The output files are:
            - Mapping_and_Variant_Statistics.tab
            - Strain_Classification.tab
*/

    publishDir "${params.outDir}/db/comparison/mtbseq/", 
        mode: 'copy',
        overwrite: true

    input:
        val(sampleID_list)
        //path(stats_files)
        //path(mtbseq_class_files)

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats

    script:
        """
        # Concatenate all files, excluding the header
        cat ${params.outDir}/db/samples/*/mtbseq/Statistics/*.tab \\
            | sed '/^Date/d' \\
            | sed "s/'//g" > Mapping_and_Variant_Statistics.tab
        
        # Concatenate all files, excluding the header
        cat ${params.outDir}/db/samples/*/mtbseq/Classification/*.tab \\
            | sed '/^Date/d' \\
            | sed "s/'//g" > Strain_Classification.tab
        """

}