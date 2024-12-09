process TBPROFILER_PROFILE_WHO {
    
    tag "$sampleID"
    
    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/tb-profiler" : "./modules/local/tbprofiler/tb-profiler.yml" }
    
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'

    publishDir "${params.outdir}/bbdd/tbprofiler/who-only", mode: 'copy'

    input:
        tuple val(sampleID), path(vcf)
        path tbprofiler_db

    output:
        path "results/${sampleID}.results.txt",  emit: tbprof_who_txt
        path "results/${sampleID}.results.json", emit: tbprof_who_json

    script:
    
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """

    # Check the DB is updated

    # Run the WHO-only analysis using the VCFs
    tb-profiler profile \\
        --vcf ${vcf} \\
        -p ${sampleID} \\
        --txt \\
        --db ${tbprofiler_db}/who \\
        --txt --dir . \\
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