process TBPROFILER_PROFILE_TBDB {
    
    tag "$sampleID"

    conda params.tbprofiler_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf)

    output:
        path("bam/tbdb-${sampleID}.bam")
        path("vcf/tbdb-${sampleID}.targets.vcf.gz")
        path("results/tbdb-${sampleID}.results.json")
        path("results/tbdb-${sampleID}.results.txt")

        // tuple for updating the sample ch
        tuple val(sampleID), 
                path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars),  
                path("results/tbdb-${sampleID}.results.txt"), // generated in this module
                path(who_out), path(mtbseq_vcf),                            emit: updated_sample_ch2

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Run TB-Proiler using TBDB database
            tb-profiler profile \\
                    -1 ${mtbc_forward} \\
                    -2 ${mtbc_reverse} \\
                -p tbdb-${sampleID} \\
                --txt --dir . \\
                --db ${params.tbprofiler_tbdb} \\
                --threads ${task.cpus} \\
                ${additional_args}

        """
}