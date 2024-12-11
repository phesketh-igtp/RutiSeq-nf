process MTBSEQ_STATS_COMPILE {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/mtbseq/", mode: 'copy'

    input:
        path stats_files
        path mtbseq_class_files

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats
        path("sample_coverage.tsv"),                    emit: sample_coverage

    script:
    """
    # Ensure the header is only included once
    head -n 1 ${stats_files[0]} > concatenated_mtbseq_stats.tsv
    
    # Concatenate all files, excluding the header
    for file in ${stats_files}
    do
        tail -n +2 \$file >> concatenated_mtbseq_stats.tsv
    done
    """

}