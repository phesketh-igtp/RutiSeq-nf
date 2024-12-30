process MTBSEQ_STATS_COMPILE {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        path stats_files
        path mtbseq_class_files
        val runID

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats

    script:
    """
    # Ensure the header is only included once
    head -n 1 ${stats_files[0]} >> Strain_Classification.tab
    
    # Concatenate all files, excluding the header
    for file in ${mtbseq_class_files[0]}
    do
        tail -n +2 \$file >> Mapping_and_Variant_Statistics.tab
    done
    """

}