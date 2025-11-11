process TBPROFILER_PROFILE_WHO {

/*
        @author: Poppy J Hesketh Best
        @date: 2025-04-08
        @version: 1.0.1
        @description: 
                This module performs TB-Profiler using TBDB results to get MT lineage and 
                resistance genes using the WHO database.
        @changelog:
            v1.0.1-2025-04-08: Fixed - correct tb-profiler db paths
*/
    
    tag "$sampleID"

    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
    
    publishDir "${params.outDir}/db/tbprofiler/who-only/", mode: 'copy'

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
        # Update the database
            tb-profiler update_tbdb --branch who

        # Run main function
            tb-profiler profile \\
                    -1 ${mtbc_forward} \\
                    -2 ${mtbc_reverse} \\
                    -p who-${sampleID} \\
                --txt --dir . \\
                --db tbdb/who ${additional_args}

        # remove the published files from the previous module:
            rm -f  ${params.outDir}/db/tbprofiler/${sampleID}_mtbc_R1.fastq.gz
            rm -f  ${params.outDir}/db/tbprofiler/${sampleID}_mtbc_R2.fastq.gz
            
        """
}

//${params.outDir}/db/tbprofiler/tbdb/who
// \\ --threads 2 --ram ${task.cpus} - this has just decided to stop working 