process TBPROFILER_PROFILE {

/*
@author: Poppy J Hesketh Best
@date: 2025-04-08
@version: 1.0
@description: 
    This module performs TB-Profiler using TBDB results to get MT lineage and 
    resistance genes using the TBDB database.
@changelog:
    v1.0.1-2025-04-08: Fixed - correct tb-profiler db paths
    v2.0.0-2025-11-13: Merged both TBDB and WHO profiling into a single module
                        Added support for single-end reads
*/

    tag "$sampleID"

    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
        
    publishDir "${params.outDir}/db/samples/${sampleID}/tbprofiler/", mode: 'copy'

    input:
        tuple val(sampleID), 
            path(fastq_1), path(fastq_2), path(mtbseq_class), 
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(snippy_vcf)

    output:
        path("results/tbdb-${sampleID}.results.json")
        path("results/tbdb-${sampleID}.results.txt")

        // tuple for updating the sample ch
        tuple val(sampleID), 
            path(fastq_1), path(fastq_2), path(mtbseq_class), 
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars),  
            path("results/tbdb-${sampleID}.results.txt"), // generated in this module
            path("results/who-${sampleID}.results.txt"), // generated in this module
            path(snippy_vcf),                            emit: updated_sample_ch2

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Update the database
            tb-profiler update_tbdb
            tb-profiler update_tbdb --branch who

        # Run TB-Proiler using TBDB database
        if [[ ! ${fastq_2}]]; then
            echo "Single-end reads detected, running TB-Profiler WHO database with single-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p tbdb-${sampleID} \\
                    --txt --dir . --platform nanopore \\
                    --db tbdb/tbdb ${additional_args}

        else
            echo "paired-end reads detected, running TB-Profiler WHO database with single-end mode."
                tb-profiler profile \\
                        -1 ${fastq_1} \\
                        -2 ${fastq_2} \\
                        -p tbdb-${sampleID} \\
                        --txt --dir . --platform illumina \\
                        --db tbdb/tbdb ${additional_args}
        fi

        # Touch the output files to the correct directory
                touch bam/tbdb-${sampleID}.bam
                touch vcf/tbdb-${sampleID}.targets.vcf.gz
                touch results/tbdb-${sampleID}.results.json
                touch results/tbdb-${sampleID}.results.txt

        # Run the WHO database next:
        if [[ ! ${fastq_2}]]; then
            echo "Single-end reads detected, running TB-Profiler WHO database with single-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p tbdb-${sampleID} \\
                    --txt --dir . --platform nanopore \\
                    --db tbdb/who ${additional_args}

        else
            echo "paired-end reads detected, running TB-Profiler WHO database with single-end mode."
                tb-profiler profile \\
                        -1 ${fastq_1} \\
                        -2 ${fastq_2} \\
                        -p tbdb-${sampleID} \\
                        --txt --dir . --platform illumina \\
                        --db tbdb/who ${additional_args}
        fi
        """
}

// ${params.outDir}/db/tbprofiler/tbdb/
// --threads 2 --ram ${task.cpus} 