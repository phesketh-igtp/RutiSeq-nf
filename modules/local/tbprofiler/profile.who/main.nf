process TBPROFILER_PROFILE_WHO {

        /*
        In this module TBProfiler is run for a single genome using the WHO database only - this is used to 
            get WHO acceptd resistance profiles. The output prefix of module file is defined as who-${sampleID} to 
            prevent overlapping file names later down the line.I am certain that I tried stageAs: as some point to circumvent this issue, 
            but it caused other problems. 
            TODO: revisit this issue.
        Like the previous modules, the input tuple for this module is the paths to all the files
            needed for a sample to proceed into the PAIRWISE_WF(), which litters the publish directory
            these excess files and are removed at the very end. This was originally done with a IF argument
            to only remove them if the files did not exist - this prevents nextflow complaining of they 
            were not generated - but this created some sybtax errors, so I changed it to touching the 
            files then removing them manually. 
            TODO: might be to revisit this at a later date.
    */
    
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
            rm -f  ${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R1.fastq.gz
            rm -f  ${params.outdir}/bbdd/tbprofiler/${sampleID}_mtbc_R2.fastq.gz
        """
}