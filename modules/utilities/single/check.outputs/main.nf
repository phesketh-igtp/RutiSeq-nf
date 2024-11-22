// modules/check_outputs.nf

process CHECK_OUTPUTS {
    tag "Checking outputs for ${sampleID}"
    
    input:
    tuple val(sampleID), 
    path(${params.outdir}/bbdd/mtbseq/samples/${sampleID})
    path()

    output:
    tuple val(sampleID), path(expected_outputs), emit: to_process
    tuple val(sampleID), path("${sampleID}_existing_outputs.txt"), emit: existing

    script:
    """
    touch ${sampleID}_existing_outputs.txt
    for output in ${expected_outputs}
    do
        if [[] -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${output}" && \ 
               -f "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${output}" ]]
        then
            echo "${output}" >> ${sampleID}_existing_outputs.txt
            rm ${output}
        fi
    done
    """
}