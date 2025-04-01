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

        errorStrategy 'ignore'

        input:
                tuple val(sampleID), 
                        path(forward),
                        path(reverse)
                path(tbprofiler_update_db)
                        
        output:
                path("bam/tbdb-${sampleID}.bam", optional: true)
                path("vcf/tbdb-${sampleID}.targets.vcf.gz", optional: true)
                path("results/tbdb-${sampleID}.results.json", optional: true)
                path("results/tbdb-${sampleID}.results.txt", optional: true)
                path("${sampleID}_tb_profiler.log", optional: true)
                path("${sampleID}_tb_profiler_status.txt")

        script:
                def additional_args = task.ext.additional_args ?: ''

        """
        # Run TB-Proiler using TBDB database
        (
                tb-profiler profile \\
                -1 ${forward} \\
                -2 ${reverse} \\
                -p tbdb-${sampleID} \\
                --txt --dir . \\
                --db ${params.outdir}/db/tbprofiler/tbdb \\
                --threads ${task.cpus} ${additional_args}
        ) > >(tee ${sampleID}_tb_profiler.log) 2>&1

        exit_status=\$?

        if [ \$exit_status -ne 0 ]; then
                echo "TB-Profiler failed for sample ${sampleID} with exit status \$exit_status" >> ${sampleID}_tb_profiler.log
                echo "FAILED" > ${sampleID}_tb_profiler_status.txt
        else
                echo "TB-Profiler completed successfully for sample ${sampleID}" >> ${sampleID}_tb_profiler.log
                
                # Check if output files have content
                if [ -s bam/tbdb-${sampleID}.bam ] && [ -s vcf/tbdb-${sampleID}.targets.vcf.gz ] && [ -s results/tbdb-${sampleID}.results.json ] && [ -s results/tbdb-${sampleID}.results.txt ]; then
                echo "SUCCESS" > ${sampleID}_tb_profiler_status.txt
                else
                echo "TB-Profiler completed but some output files are empty. This may be expected for negative controls or samples with insufficient data." >> ${sampleID}_tb_profiler.log
                echo "LOW_DATA" > ${sampleID}_tb_profiler_status.txt
                fi
        fi

        # Ensure output files exist (even if empty) to satisfy Nextflow
        touch bam/tbdb-${sampleID}.bam
        touch vcf/tbdb-${sampleID}.targets.vcf.gz
        touch results/tbdb-${sampleID}.results.json
        touch results/tbdb-${sampleID}.results.txt

        # Always exit with status 0 to prevent pipeline failure
        exit 0
        """

}