process CN_TBPROFILER_PROFILE_TBDB {

        tag "$sampleID"

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
                        path(reverse),
                        path(qc_results)

        output:
                path("bam/tbdb-${sampleID}.bam", optional: true)
                path("vcf/tbdb-${sampleID}.targets.vcf.gz", optional: true)
                path("results/tbdb-${sampleID}.results.json", optional: true)
                path("results/tbdb-${sampleID}.results.txt", optional: true)

        script:
                def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

                """
                # Run TB-Proiler using TBDB database
                        tb-profiler profile \\
                                -1 ${forward} \\
                                -2 ${reverse} \\
                                -p tbdb-${sampleID} \\
                                --txt --dir . \\
                                --db ${params.tbprofiler_tbdb} \\
                                --threads ${task.cpus} ${additional_args} \\
                1>>.command.out \\
                2>>.command.err || true

                touch bam/tbdb-${sampleID}.bam
                touch vcf/tbdb-${sampleID}.targets.vcf.gz
                touch results/tbdb-${sampleID}.results.json
                touch results/tbdb-${sampleID}.results.txt
                """
}