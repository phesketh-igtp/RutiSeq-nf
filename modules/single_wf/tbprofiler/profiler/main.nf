process TBPROFILER_PROFILE {

    tag "$sampleID"

    conda params.tbprofiler_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
            else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
            else return null
        }
        
    publishDir "${params.outDir}/db/samples/${sampleID}/tbprofiler/", 
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
        // Main outputs
        path("tbprofiler.csv")
        path("results/tbdb-${sampleID}.results.*")
        path("results/who-${sampleID}.results.*")

        // tuple for updating the sample ch
        tuple val(sampleID), 
            path(fastq_1), 
            path(fastq_2),
            val(type),
            path(mtbseq_class), 
            path(mtbseq_stats), 
            path(mtbseq_pos), 
            path(mtbseq_vars),  
            path("results/tbdb-${sampleID}.results.txt"), // generated in this module
            path("results/who-${sampleID}.results.txt"), // generated in this module
            path(snippy_vcf), emit: updated_sample_ch2

    script:
        def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

        """
        # Create output directories
        version=\$(tb-profiler --version | sed 's@version @v.@g')
        tbdb_v=\$(tb-profiler list_db | grep 'tbdb' | cut -f2)
        who_v=\$(tb-profiler list_db | grep 'who_v2+' | cut -f2)
        version_string=\$(echo "\${version} (TBDB:\${tbdb_v} WHO:\${who_v})")
        safe_version=\$(printf '%s\n' "\$version_string" | sed 's/[&|]/\\&/g')

        # Check if TB-Profiler has already been run for this sample by looking for key output files. 
        ## handles situations when the workflow is re-run and prevent each single-wf steps 
        ## from re-running unnecessarily.

        if [[ -f "${params.outDir}/db/samples/${sampleID}/tbprofiler/tbprofiler.csv"  \\
                &&-f "${params.outDir}/db/samples/${sampleID}/tbprofiler/results/tbdb-${sampleID}.results.txt" \\
                && -f "${params.outDir}/db/samples/${sampleID}/tbprofiler/results/who-${sampleID}.results.txt"
                ]]; then
            
            echo -e "TB-Profiler results already exist for sample ${sampleID}, skipping TB-Profiler profiling step..."
            echo -e ""
            echo -e "Copying over TBDB and WHO results"
            mkdir results/
            ln -s ${params.outDir}/db/samples/${sampleID}/tbprofiler/tbprofiler.csv .
            ln -s ${params.outDir}/db/samples/${sampleID}/tbprofiler/results/tbdb-${sampleID}.results.txt results/
            ln -s ${params.outDir}/db/samples/${sampleID}/tbprofiler/results/who-${sampleID}.results.txt results/

        else

            echo "Running TB-Profiler profiling for sample ${sampleID}."    

            # Run TB-Profiler using TBDB database
            if [[ ! -f "${fastq_2}" ]]; then
                echo "Single-end reads detected, 
                running TB-Profiler TBDB database with single-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p tbdb-${sampleID} \\
                    --no_trim \\
                    --txt --platform nanopore \\
                    --dir . --db "tbdb" \\
                    --threads ${params.cpus} \\
                    ${additional_args}

                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -p who-${sampleID} \\
                    --no_trim \\
                    --txt --platform nanopore \\
                    --dir . --db "who_v2+" \\
                    --threads ${params.cpus} \\
                    ${additional_args}

                # collate the results into a single file
                tb-profiler collate --format csv

                # Add version value to all rows
                awk -v v="\${safe_version}" 'NR==1 {print "versions," \$0; next} \\
                    {print v "," \$0}' tbprofiler.csv \\
                    > tmp && mv tmp tbprofiler.csv

            else
                echo "Paired-end reads detected,
                running TB-Profiler TBDB database with paired-end mode."
                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -2 ${fastq_2} \\
                    --no-trim \\
                    -p tbdb-${sampleID} \\
                    --txt --platform illumina \\
                    --dir . --db "tbdb" \\
                    --threads ${params.cpus} \\
                    ${additional_args}

                tb-profiler profile \\
                    -1 ${fastq_1} \\
                    -2 ${fastq_2} \\
                    -p who-${sampleID} \\
                    --no_trim \\
                    --txt --platform illumina \\
                    --dir . --db "who_v2+" \\
                    --threads ${params.cpus} \\
                    ${additional_args}

                # collate the results into a single file
                tb-profiler collate --format csv

                # Add version value to all rows
                awk -v v="\${safe_version}" 'NR==1 {print "versions," \$0; next} \\
                    {print v "," \$0}' tbprofiler.csv \\
                    > tmp && mv tmp tbprofiler.csv

            fi

        fi

        echo "TB-Profiler complete"
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-06-03
@version: 2.2.0
@description: 
    This module performs TB-Profiler using TBDB results to get MT lineage and 
    resistance genes using the TBDB database.
@changelog:
    v1.0.1-2025-04-08: Fixed - correct tb-profiler db paths
    v2.0.0-2025-11-13: Merged both TBDB and WHO profiling into a single module
                        Added support for single-end reads
    v2.1.0-2026-01-19: Updated to check for existing results to avoid re-running
    v2.1.1-2026-03-04: Corrected some file movement issues cauring errors and collapsing the module.
    v2.2.0-2026-06-03: Moved collate into this module to avoid using TB-Profiler in downstream 
                    (pairwise) summarisation, due to incompatibility between software 
                    and database versions.
    v2.2.1-2026-06-05: Simplified structure of the module to prevent nested files
*/