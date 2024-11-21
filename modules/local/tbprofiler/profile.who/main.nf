process TBPROFILER_PROFILE_WHO {
    tag "$sampleID"
    label 'process_medium'

    conda "bioconda::tb-profiler=6.3.0"

    container 'oras://community.wave.seqera.io/library/tb-profiler:6.3.0--4f362e6be5d39a05'

    publishDir "${params.outdir}/bbdd/tbprofiler/who-only", mode: 'copy'

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
    """

    stub:
    """
    mkdir -p bam vcf results
    touch results/${sampleID}.results.txt
    touch results/${sampleID}.results.json
    """
}