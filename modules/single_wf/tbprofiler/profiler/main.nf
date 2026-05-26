process TBPROFILER_PROFILE {

    tag "$sampleID"

    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
        
    publishDir "${params.outDir}/db/samples/${sampleID}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2),
            val(type),
            path(mtbseq_class), 
            path(mtbseq_stats), 
            path(mtbseq_pos), 
            path(mtbseq_vars), 
            path(tbdb_out), 
            path(who_out), 
            path(snippy_vcf)
        path(tbprofiler_db)

    output:
        path("tbprofiler/*")

        // tuple for updating the sample ch
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2),
            val(type),
            path(mtbseq_class), 
            path(mtbseq_stats), 
            path(mtbseq_pos), 
            path(mtbseq_vars),  
            path("tbprofiler/tbdb-${sampleID}.results.txt"), // generated in this module
            path("tbprofiler/who-${sampleID}.results.txt"), // generated in this module
            path(snippy_vcf), emit: updated_sample_ch2

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Create output directories
        mkdir -p tbprofiler/

        # Check if TB-Profiler has already been run for this sample by looking for key output files. 
        ## handles situations when the workflow is re-run and prevent each single-wf steps 
        ## from re-running unnecessarily.

        if [[ -f "${params.outDir}/db/samples/${sampleID}/tbprofiler/tbdb-${sampleID}.results.txt" \\
                && -f "${params.outDir}/db/samples/${sampleID}/tbprofiler/who-${sampleID}.results.txt" \\
                ]]; then
            
            echo -e "TB-Profiler results already exist for sample ${sampleID}, skipping TB-Profiler profiling step..."
            echo -e ""
            echo -e "Copying over TBDB and WHO results"
            cp ${params.outDir}/db/samples/${sampleID}/tbprofiler/tbdb-${sampleID}.results.txt tbprofiler/
            cp ${params.outDir}/db/samples/${sampleID}/tbprofiler/who-${sampleID}.results.txt tbprofiler/

        else

            echo "Running TB-Profiler profiling for sample ${sampleID}."    

            # Run TB-Profiler using TBDB database
            if [[ ! -f "${fastq_2}" ]]; then
                echo "Single-end reads detected, 
                running TB-Profiler TBDB database with single-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p tbdb-${sampleID} \\
                    --txt --platform nanopore \\
                    --dir tbprofiler/ \\
                    --call_whole_genome \\
                    --db tbdb/tbdb ${additional_args}

                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p who-${sampleID} \\
                    --txt --platform nanopore \\
                    --db tbdb/who \\
                    --dir tbprofiler/ \\
                    --call_whole_genome \\
                    ${additional_args}

                cp tbprofiler/results/* tbprofiler/
                
            else
                echo "Paired-end reads detected,
                running TB-Profiler TBDB database with paired-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -2 ${fastq_2} \\
                    -p tbdb-${sampleID} \\
                    --txt --platform illumina \\
                    --db tbdb/tbdb \\
                    --dir tbprofiler/ \\
                    --call_whole_genome \\
                    ${additional_args}
                
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -2 ${fastq_2} \\
                    -p who-${sampleID} \\
                    --txt --platform illumina \\
                    --db tbdb/who \\
                    --dir tbprofiler/ \\
                    --call_whole_genome \\
                    ${additional_args}

                cp tbprofiler/results/* tbprofiler/

            fi

        fi

        echo "TB-Profiler complete"
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-01-19
@version: 2.1.0
@description: 
    This module performs TB-Profiler using TBDB results to get MT lineage and 
    resistance genes using the TBDB database.
@changelog:
    v1.0.1-2025-04-08: Fixed - correct tb-profiler db paths
    v2.0.0-2025-11-13: Merged both TBDB and WHO profiling into a single module
                        Added support for single-end reads
    v2.1.0-2026-01-19: Updated to check for existing results to avoid re-running
    v2.1.1-2026-03-04: Corrected some file movement issues cauring errors and collapsing the module.
*/