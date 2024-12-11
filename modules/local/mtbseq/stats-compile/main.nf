process MTBSEQ_STATS_COMPILE {

    tag "${runID}"

    publishDir "${params.outdir}/bbdd/mtbseq/", mode: 'copy'

    input:
        val runID
        path strain_classifications
        path mapping_statistics

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats
        path("sample_coverage.tsv"),                    emit: sample_coverage

    script:
    """

    # Process Strain Classification files
    head -n 1 ${strain_classifications[0]} > Strain_Classification.tab
    for file in ${strain_classifications}; do
        tail -n +2 \$file >> Strain_Classification.tab
    done

    # Process Mapping and Variant Statistics files
    head -n 1 ${mapping_statistics[0]} > Mapping_and_Variant_Statistics.tab
    for file in ${mapping_statistics}; do
        tail -n +2 \$file >> Mapping_and_Variant_Statistics.tab
    done

    # Create tuple_mapping_stats.tsv
    cut -f4,19 Mapping_and_Variant_Statistics.tab | sed '1d' > sample_coverage.tsv

    """

}