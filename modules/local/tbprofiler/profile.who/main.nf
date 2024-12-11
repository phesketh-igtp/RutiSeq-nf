process TBPROFILER_PROFILE_WHO {
    
    tag "$sampleID"
    
    conda 'bioconda::tb-profiler==6.5.0'

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/who-only", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)

    output:
        path "results/${sampleID}.results.txt",     emit: who_out
        path "results/${sampleID}.results.json"

    script:
    
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """

    # Check the DB is updated

    # Run the WHO-only analysis using the VCFs
    tb-profiler profile \\
        -1 ${mtbc_forward} \\
        -2 ${mtbc_reverse} \\
        -p ${sampleID} \\
        --txt --dir . \\
        --db ${params.tbprofiler_who} \\
        --threads ${task.cpus} \\
        ${additional_args}

    rm -rf bam/ vcf/
    """

    stub:
    """
    mkdir -p bam vcf results
    touch results/${sampleID}.results.txt
    touch results/${sampleID}.results.json
    """
}