process FILE_CHECK {
    tag "$sampleID"

    input:
        tuple val(sampleID), path(forward), path(reverse)

    output:
        tuple val(sampleID), 
            path("${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab"), 
            path("${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab"), 
            path("${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab"),
            path("${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_*.tab"),
            path("${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt"),
            path("${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt"),      emit: pairwise_input, optional: true
        tuple val(sampleID), path(forward), path(reverse),                                          emit: single_input, optional: true

    script:
        """
        if [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab" ] && \
        [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab" ] && \
        [ -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab" ] && \
        compgen -G "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_*.tab" > /dev/null && \
        [ -f "${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt" ] && \
        [ -f "${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt" ]; then
            
            echo "${sampleID}" >> pairwise_samples.txt
            # Create symbolic links to the existing files
            ln -s "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Classification/Strain_Classification.tab" .
            ln -s "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Statistics/Mapping_and_Variant_Statistics.tab" .
            ln -s "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Position_Tables/${sampleID}.gatk_position_table.tab" .
            ln -s "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Called/${sampleID}.gatk_position_variants_"*.tab .
            ln -s "${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt" .
            ln -s "${params.outdir}/bbdd/tbprofiler/who-only/results/${sampleID}.results.txt" .
            
            # Create a flag file for pairwise samples
            touch pairwise_${sampleID}.flag
        else
            echo "${sampleID}" >> single_samples.txt
            # Create a flag file for single samples
            touch single_${sampleID}.flag
        fi
        """
}