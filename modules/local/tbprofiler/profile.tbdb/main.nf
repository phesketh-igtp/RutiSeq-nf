process TBPROFILER_PROFILE_TBDB {
    
    tag "$sampleID"

    conda params.tbprofiler_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
            } else { 'quay.io/biocontainers/tb-profiler' }
    }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)

    output:
        tuple val(sampleID), path("bam/${sampleID}.bam")
        tuple val(sampleID), path("vcf/${sampleID}.targets.vcf.gz"),    emit: tbdb_vcf
        tuple val(sampleID), path("results/${sampleID}.results.txt"),   emit: tbdb_out
        tuple val(sampleID), path("results/${sampleID}.results.json")

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        echo "Debug: Checking input files"
        ls -l ${mtbc_forward} ${mtbc_reverse}

        ln -s ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_R1.fastq.gz ${sampleID}_R1.fastq.gz
        ln -s ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_R2.fastq.gz ${sampleID}_R2.fastq.gz 

        echo "Debug: Running tb-profiler"
        tb-profiler profile \\
            -1 ${sampleID}_R1.fastq.gz \\
            -2 ${sampleID}_R2.fastq.gz  \\
            -p ${sampleID} \\
            --txt --dir . \\
            --db ${params.tbprofiler_tbdb} \\
            --threads ${task.cpus} \\
            ${additional_args}

        echo "Debug: Listing contents of working directory after tb-profiler"
        ls -l
        """
}