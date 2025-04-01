process TBPROFILER_PROFILE_TBDB {

        /*
                This module performs TB-Profiler using TBDB results to get MT lineage and 
                identify any potential contamination in the genome
        */

        tag "$sampleID"

        conda params.tbprofiler_env

        container { 
                if (workflow.containerEngine == 'singularity') return params.singularity_tbprofiler
                else if (workflow.containerEngine == 'docker') return params.docker_tbprofiler
                else if (workflow.containerEngine == 'apptainer') return params.apptainer_tbprofiler
                else return null
        }
        
        publishDir "${params.outdir}/bbdd/tbprofiler/", mode: 'copy'

        input:
                tuple val(sampleID), 
                        path(mtbc_forward), path(mtbc_reverse), path(mtbseq_class), 
                        path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                        path(tbdb_out), path(who_out), path(mtbseq_vcf)
                path(tbprofiler_update_handover)

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
                path(who_out), path(mtbseq_vcf),                            emit: updated_sample_ch2

        script:
                def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

                """
                # Run TB-Proiler using TBDB database
                        tb-profiler profile \\
                                -1 ${mtbc_forward} \\
                                -2 ${mtbc_reverse} \\
                                -p tbdb-${sampleID} \\
                                --txt --dir . \\
                                --threads ${task.cpus} ${additional_args}

                        touch bam/tbdb-${sampleID}.bam
                        touch vcf/tbdb-${sampleID}.targets.vcf.gz
                        touch results/tbdb-${sampleID}.results.json
                        touch results/tbdb-${sampleID}.results.txt


                # remove the published files from the previous module:
                        rm -f ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R1.fastq.gz
                        rm -f ${params.outdir}/bbdd/read-qc/mtbc_reads/${sampleID}_mtbc_R2.fastq.gz
                """
}

//--db ${params.tbprofiler_tbdb}