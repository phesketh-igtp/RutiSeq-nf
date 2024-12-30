process FILE_CHECK {
    tag "$sampleID"

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        path("pairwise_samples.txt"),   emit: pairwise_input,   optional: true
        path("single_samples.txt"),     emit: single_input,     optional: true

    script:

        def forward_path = forward.toRealPath()
        def reverse_path = reverse.toRealPath()
        def mtbseq_class = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab"
        def mtbseq_stats = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab"
        def mtbseq_pos   = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab"
        def mtbseq_vars  = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_cf4_cr4_fr75_ph4_outmode011.tab"
        def tbdb_out     = "${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt"
        def who_out      = "${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt"
        def mtbseq_vcf   = "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/${sampleID}.gatk.vcf.gz"

        """
        if [ -f ${mtbseq_class} ] && [ -f ${mtbseq_stats} ] && [ -f ${mtbseq_pos} ] && [ -f ${mtbseq_vars} ] && [ -f ${tbdb_out} ] && [ -f ${who_out} ] && [ -f ${mtbseq_vcf} ]; then
            
            echo "${sampleID},${mtbseq_class},${mtbseq_stats},${mtbseq_pos},${mtbseq_vars},${tbdb_out},${who_out},${mtbseq_vcf}" > pairwise_samples.txt
            
            echo "DEBUG: Added to pairwise_samples.txt: ${sampleID}" >&2
        
        else
            
            echo "${sampleID},${forward_path},${reverse_path}" > single_samples.txt
            
            echo "DEBUG: Added to single_samples.txt: ${sampleID}" >&2
        
        fi

        echo "DEBUG: Content of pairwise_samples.txt:" >&2
        cat pairwise_samples.txt >&2 || echo "pairwise_samples.txt does not exist" >&2
        echo "DEBUG: Content of single_samples.txt:" >&2
        cat single_samples.txt >&2 || echo "single_samples.txt does not exist" >&2

        # Ensure the process doesn't fail if one of the files doesn't exist
            touch pairwise_samples.txt single_samples.txt
        """
}