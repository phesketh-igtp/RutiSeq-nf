process MTBSEQ_STATS_COMPILE {

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'move'

    input:
        path stats_files
        path mtbseq_class_files

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats

    script:
        """

        # Concatenate all files, excluding the header
        cat ${params.outdir}/bbdd/mtbseq/samples/*/Statistics/*.tab \\
                            | sed '/^Date/d' \\
                            | sed "s/'//g" > Mapping_and_Variant_Statistics.tab
        
        # Concatenate all files, excluding the header
        cat ${params.outdir}/bbdd/mtbseq/samples/*/Classification/*.tab \\
                            | sed '/^Date/d' \\
                            | sed "s/'//g" > Strain_Classification.tab

        """

}