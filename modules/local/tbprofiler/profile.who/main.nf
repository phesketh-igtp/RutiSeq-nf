process TBPROFILER_PROFILE_WHO {
    
    tag "$sampleID"
    
    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler" : "./modules/local/tbprofiler/tb-profiler.yml" }
   
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'move'

    input:
    tuple val(sampleID), path(vcf)

    output:
    path "results/${sampleID}.results.txt",  emit: tbprof_tbdb_res
    path "results/${sampleID}.results.json", emit: tbprof_tbdb_json

    script:
    def args = task.ext.args ?: ""
    """
    tb-profiler profile \\
        --vcf ${vcf} \\
        -p ${sampleID} \\
        --txt \\
        --db who \\
        --dir . \\
        ${args} \\
        --threads ${task.cpus} \\
        --ram ${task.memory.toGiga()}

    rm -rf bam/ vcf/
    """

    stub:
    """
    mkdir -p bam vcf results
    touch results/${sampleID}.results.txt
    touch results/${sampleID}.results.json
    """
}