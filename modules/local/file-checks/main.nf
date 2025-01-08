process FILE_CHECK {
    tag "$sampleID"

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        path("samples.tuple.csv"), emit: sample_paths

    script:
        def forward_path = forward.toRealPath()
        def reverse_path = reverse.toRealPath()
        def mtbseq_class = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab"
        def mtbseq_stats = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab"
        def mtbseq_pos   = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab"
        def mtbseq_vars  = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_cf4_cr4_fr75_ph4_outmode011.tab"
        def tbdb_out     = "${params.outdir}/bbdd/tbprofiler/results/tbdb-${sampleID}.results.txt"
        def who_out      = "${params.outdir}/bbdd/tbprofiler/who-only/results/who-${sampleID}.results.txt"
        def mtbseq_vcf   = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}.gatk.vcf.gz"

        """
        if [ -f ${mtbseq_class} ] && [ -f ${mtbseq_stats} ] && [ -f ${mtbseq_pos} ] && [ -f ${mtbseq_vars} ] && [ -f ${tbdb_out} ] && [ -f ${who_out} ] && [ -f ${mtbseq_vcf} ]; then       
            echo "${sampleID},,,${mtbseq_class},${mtbseq_stats},${mtbseq_pos},${mtbseq_vars},${tbdb_out},${who_out},${mtbseq_vcf}" > samples.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2    
        else
            echo "${sampleID},${forward_path},${reverse_path},,,,,,," > samples.tuple.csv
            echo "DEBUG: Added to samples.txt: ${sampleID}" >&2
        fi

        # Ensure the process doesn't fail if one of the files doesn't exist
        touch samples.tuple.csv
        """
}