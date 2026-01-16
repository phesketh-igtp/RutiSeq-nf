process TBPROFILER_PROFILE_TBDB {

/*
@author: Poppy J Hesketh Best
@date: 2025-04-08
@version: 1.0
@description: 
    This module performs TB-Profiler using TBDB results to get MT lineage and 
    resistance genes using the TBDB database.
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
        
    publishDir "${params.outDir}/db/samples/${sampleID}/tbprofiler/", mode: 'move'

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
            path("results/who-${sampleID}.results.txt"), // generated in this module
            path(mtbseq_vcf),                            emit: updated_sample_ch2

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Update the database
            tb-profiler update_tbdb
            tb-profiler update_tbdb --branch who

        # Run TB-Proiler using TBDB database
                tb-profiler profile \\
                    -1 ${mtbc_forward} \\
                    -2 ${mtbc_reverse} \\
                    -p tbdb-${sampleID} \\
                    --txt --dir . --platform illumina \\
                    --db tbdb/tbdb ${additional_args}

        # Touch the output files to the correct directory
                touch bam/tbdb-${sampleID}.bam
                touch vcf/tbdb-${sampleID}.targets.vcf.gz
                touch results/tbdb-${sampleID}.results.json
                touch results/tbdb-${sampleID}.results.txt

        # Run the WHO database next:
            tb-profiler profile \\
                    -1 ${mtbc_forward} \\
                    -2 ${mtbc_reverse} \\
                    -p tbdb-${sampleID} \\
                    --txt --dir . --platform illumina \\
                    --db tbdb/who ${additional_args}

        # remove the published files from the previous module:
            rm -f ${params.outDir}/db/read-qc/mtbc_reads/${sampleID}_mtbc_R1.fastq.gz
            rm -f ${params.outDir}/db/read-qc/mtbc_reads/${sampleID}_mtbc_R2.fastq.gz
        """
}

// ${params.outDir}/db/tbprofiler/tbdb/
// --threads 2 --ram ${task.cpus} 