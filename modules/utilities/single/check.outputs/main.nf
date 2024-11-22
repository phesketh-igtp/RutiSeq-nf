// modules/check_outputs.nf

process CHECK_EXISTING_OUTPUTS {
    tag "Checking outputs for ${sampleID}"
    
    input:
    tuple val(sampleID), 
    path(${params.outdir}/bbdd/mtbseq/samples/${sampleID})
    path(${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.txt)
    path(${params.outdir}/bbdd/tbprofiler/results/${sampleID}.results.json)

    output:

    script:
    """

    """
}