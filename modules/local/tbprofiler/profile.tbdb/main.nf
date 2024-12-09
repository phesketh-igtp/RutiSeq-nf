process TBPROFILER_PROFILE_TBDB {
    
    tag "$sampleID"

    //conda "../conda/tb-profiler.yml"
    conda '/imppc/labs/emlab/phesketh/miniforge3/envs/tb-profiler'

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)
        path(tbprofiler_db)

    output:
        tuple val(sampleID), path("bam/${sampleID}.bam"), emit: tbprof_tbdb_bam
        tuple val(sampleID), path("vcf/${sampleID}.targets.vcf.gz"), emit: tbprof_tbdb_vcf
        tuple val(sampleID), path("results/${sampleID}.results.txt"), emit: tbprof_tbdb_res
        tuple val(sampleID), path("results/${sampleID}.results.json"), emit: tbprof_tbdb_json

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    tb-profiler profile \\
        -1 ${mtbc_forward} \\
        -2 ${mtbc_reverse} \\
        -p ${sampleID} \\
        --txt --dir . \\
        --db ${tbprofiler_db}/tbdb \\
        --threads ${task.cpus} \\
        ${additional_args}

    """

    stub:
    """
    mkdir -p bam vcf results
    touch bam/${sampleID}.bam
    touch vcf/${sampleID}.targets.vcf.gz
    touch results/${sampleID}.results.txt
    touch results/${sampleID}.results.json
    """
}