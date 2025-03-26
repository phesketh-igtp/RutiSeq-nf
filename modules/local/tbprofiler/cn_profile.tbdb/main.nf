process CN_TBPROFILER_TBDB {

        tag "$sampleID"
        
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
                        path(reverse),
                        path(qc_results)

        output:
                path("bam/tbdb-${sampleID}.bam", optional: true)
                path("vcf/tbdb-${sampleID}.targets.vcf.gz", optional: true)
                path("results/tbdb-${sampleID}.results.json", optional: true)
                path("results/tbdb-${sampleID}.results.txt", optional: true)
                path("${sampleID}_tb_profiler.log", optional: true)

        script:
                def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

                """
                # Run TB-Proiler using TBDB database
                (
                        tb-profiler profile \\
                                -1 ${forward} \\
                                -2 ${reverse} \\
                                -p tbdb-${sampleID} \\
                                --txt --dir . \\
                                --db ${params.tbprofiler_tbdb} \\
                                --threads ${task.cpus} ${additional_args}
                ) > >(tee ${sampleID}_tb_profiler.log) 2>&1

                exit_status=\$?

                if [ \$exit_status -ne 0 ]; then
                        echo "TB-Profiler failed for sample ${sampleID} with exit status \$exit_status" >> ${sampleID}_tb_profiler.log
                else
                        echo "TB-Profiler completed successfully for sample ${sampleID}" >> ${sampleID}_tb_profiler.log
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