process CN_TBPROFILER_TBDB {

/*
        @author: Poppy J Hesketh Best
        @date: 2025-04-01
        @version: 0.1
        @description: 
                This module runs TB-Profiler on a single sample using the TBDB database. 
                It is designed to be used in the context of the negative control workflow. 
                It takes a tuple of sampleID, forward read file, and reverse read file as 
                input. The output is the TB-Profiler classification and statistics files.
                The TBDB database is a custom database for TB-Profiler that is used to 
                identify Mycobacterium tuberculosis complex (MTBC) strains and their 
                resistance profiles.
                In the negative control, --no-delly is used as this causes a lot of failure
                for the fasta files when they have very few reads
*/

        tag "${sampleID}"
        
        conda params.tbprofiler_env

        container { 
                if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
                else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
                else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
                else return null
        }
        
        publishDir "${params.outdir}/bbdd/negative-controls/tbprofiler/", mode: 'copy'

        input:
                tuple val(sampleID), 
                        path(forward),
                        path(reverse)
                path(tbprofiler_update_db)
                        
        output:
                path("results/tbdb-${sampleID}.results.txt"), emit: tbprofiler_results
                path("results/tbdb-${sampleID}.results.json"), optional: true
                path("bam/*"), optional: true
                path("vcf/*"), optional: true
                path("${sampleID}_tb_profiler_status.txt")

        script:
                def additional_args = task.ext.additional_args ?: ''

        """
        # Run TB-Profiler using TBDB database
        set +e #tells the shell not to exit immediately if a command fails
        tb-profiler profile \\
                -1 ${forward} \\
                -2 ${reverse} \\
                -p tbdb-${sampleID} \\
                --txt --dir . \\
                --db ${params.outdir}/db/tbprofiler/tbdb \\
                --threads ${task.cpus} \\
                --no_delly ${additional_args} \\
                > tb-profiler.out 2> tb-profiler.err
        set -e #restores default command fails checks

        # Check if the results file was created
        if [[ -f results/tbdb-${sampleID}.results.txt ]]; then
                echo "${sampleID},SUCCESS" > ${sampleID}_tb_profiler_status.txt
        else
                echo "${sampleID},FAILED" > ${sampleID}_tb_profiler_status.txt
                # Create empty files to satisfy output requirements
                echo "" > results/tbdb-${sampleID}.results.txt
                echo "" > results/tbdb-${sampleID}.results.json
        fi

        # Always exit with status 0 to prevent pipeline failure
        exit 0
        """

}