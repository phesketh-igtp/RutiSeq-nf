process TBPROFILER_PROFILE_TBDB {
    tag "$sampleID"
    label 'process_medium'

    conda "bioconda::tb-profiler=6.3.0"
    container 'community.wave.seqera.io/library/tb-profiler:6.3.0--35b5c369eb6e0d52'

    input:
    tuple val(old_name), val(sampleID), path(forward), path(reverse)

    output:
    path "bam/${sampleID}.bam",              emit: tbprof_tbdb_bam
    path "vcf/${sampleID}.target.vcf.gz",    emit: tbprof_tbdb_vcf
    path "results/${sampleID}.results.txt",  emit: tbprof_tbdb_res
    path "results/${sampleID}.results.json", emit: tbprof_tbdb_json

    script:
    def args = task.ext.args ?: ""
    """
    tb-profiler profile \\
        -1 ${forward} \\
        -2 ${reverse} \\
        -p ${sampleID} \\
        --txt \\
        --dir . \\
        ${args} \\
        --threads ${task.cpus} \\
        --ram ${task.memory.toGiga()}
    """

    stub:
    """
    mkdir -p bam vcf results
    touch bam/${sampleID}.bam
    touch vcf/${sampleID}.target.vcf.gz
    touch results/${sampleID}.results.txt
    touch results/${sampleID}.results.json
    """
}