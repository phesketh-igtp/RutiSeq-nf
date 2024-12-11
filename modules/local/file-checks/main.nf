process FILE_CHECK {
    tag "$sampleID"

    input:
    tuple val(sampleID), path(forward), path(reverse)

    output:
        path("pairwise_samples.txt"),   emit: pairwise_input,   optional: true
        path("single_samples.txt"),     emit: single_input,     optional: true

    script:

    """
    mtbseq_class=${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab)
    mtbseq_stats=${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab
    mtbseq_pos=${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab
    mtbseq_vars=${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_cf4_cr4_fr75_ph4_outmode001.tab
    tbdb_out=${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt
    who_out=${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt

    if [ -f \${mtbseq_class} ] && [ -f \${mtbseq_stats} ] && [ -f \${mtbseq_pos} ] && [ -f \${mtbseq_vars} ] && [ -f \${tbdb_out} ] && [ -f \${who_out} ]; then
        
        echo "${sampleID},\$(realpath \$mtbseq_class),\$(realpath \$mtbseq_stats),\$(realpath \$mtbseq_pos),\$(realpath \$mtbseq_vars),$(realpath \$tbdb_out),\$(realpath \$who_out)" > pairwise_samples.txt

    else

        echo "${sampleID},${forward},${reverse}" > single_samples.txt

    fi
    """
}