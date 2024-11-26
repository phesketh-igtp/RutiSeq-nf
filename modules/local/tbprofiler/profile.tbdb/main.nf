process TBPROFILER_PROFILE_TBDB {
    
    tag "$sampleID"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler" : "./modules/local/tbprofiler/tb-profiler.yml" }

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'link'

    input:
        tuple val(sampleID), path(forward), path(reverse)
        path tbprofiler_db

    output:
        tuple val(sampleID), path("bam/${sampleID}.bam"),              emit: tbprof_tbdb_bam
        tuple val(sampleID), path("vcf/${sampleID}.targets.vcf.gz"),    emit: tbprof_tbdb_vcf
        tuple val(sampleID), path("results/${sampleID}.results.txt"),  emit: tbprof_tbdb_res
        tuple val(sampleID), path("results/${sampleID}.results.json"), emit: tbprof_tbdb_json

    script:
    def args = task.ext.args ?: ""
    """
    tb-profiler profile \\
        -1 ${forward.toRealPath()} \\
        -2 ${reverse.toRealPath()} \\
        -p ${sampleID} \\
        --txt --dir . \\
        --db ${tbprofiler_db}/tbdb \\
        --threads ${task.cpus} \\
        ${args}

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