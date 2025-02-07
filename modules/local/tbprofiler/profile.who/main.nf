process TBPROFILER_PROFILE_WHO {
    
    tag "$sampleID"

    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
    
    publishDir "${params.outdir}/bbdd/tbprofiler/who-only/", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf)

    output:
        path("results/who-${sampleID}.results.json")
        path("results/who-${sampleID}.results.txt")

        // tuple for updating the sample ch
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars),  
                path(tbdb_out), 
                path("results/who-${sampleID}.results.txt"), // generated in this module
                path(mtbseq_vcf),                            emit: updated_sample_ch3

    script:
    
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Run main function
        tb-profiler profile \\
                -1 ${mtbc_forward} \\
                -2 ${mtbc_reverse} \\
                -p who-${sampleID} \\
            --txt --dir . \\
            --db ${params.tbprofiler_who} \\
            --threads ${task.cpus} \\
            ${additional_args}


        # remove the published files from the previous module:
            for file in \\
                "${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R1.fastq.gz" \\
                "${params.outdir}/bbdd/tbprofiler/who-only/${sampleID}_mtbc_R2.fastq.gz" \\
                "${params.outdir}/bbdd/tbprofiler/who-only/tbdb-${sampleID}.results.txt";
            do
                if [ -f "\${file}" ] || [ -e "\${file}" ]; then rm "\${file}"; fi
            done
        """
}