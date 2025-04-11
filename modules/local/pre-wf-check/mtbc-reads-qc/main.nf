process MTBC_READ_QC {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        This process was originally used to run the MTBC_READ_QC step of the pipeline.
        It has been modified to remove the kaiju step of the pipeline. And just generates
        sequencing statistics and down-samples the reads to 5,000,000 reads if there are
        more than 5,000,000 reads. The down-sampled reads are then used for the TBProfiler
        and MTBseq steps of the pipeline.
        TODO: Add kraken2 step back in to the pipeline (BUT! only to classify the reads and not
        to partition the reads. This is because the MTBseq pipeline generated much lower quality
        results when the reads were partitioned.)
*/
    
    tag "$sampleID"
    
    conda params.taxonomy_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_kaiju
            else if (workflow.containerEngine == 'docker') return params.docker_kaiju
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_kaiju
            else return null
        }

    publishDir "${params.outDir}/bbdd/read-qc/", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(forward), path(reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf)

    output:
        path("tables/${sampleID}.kaiju.out"),           optional: true
        path("tables/${sampleID}.kaiju_summary.tsv"),   optional: true

        // Emit ch for compiling read-QC
        tuple val(sampleID), path("tables/${sampleID}.qc.out"),      emit: qc_results, optional: true

        // Emit ch for the updated channel with all the outputs
        tuple val(sampleID), 
                path("mtbc_reads/${sampleID}_mtbc_R1.fastq.gz"), 
                path("mtbc_reads/${sampleID}_mtbc_R2.fastq.gz"), 
                path(mtbseq_class), path(mtbseq_stats), 
                path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf),     emit: updated_sample_ch1

    script:
        def additional_args_kaiju       = task.ext.additional_args_kaiju ?: ''
        def additional_args_kaiju2table = task.ext.additional_args_kaiju2table ?: ''

    """
        mkdir -p tables
        touch tables/${sampleID}.kaiju.out
        touch tables/${sampleID}.kaiju_summary.tsv
        touch tables/${sampleID}.qc.out

        read_count=\$(seqkit stats -abT -j ${task.cpus} ${forward} | sed "1d" | cut -f4)

        mkdir -p mtbc_reads

        if [[ \${read_count} > 5000000 ]]; then
            echo -e "Downsampling to 5,000,000 reads for TBProfiler/MTBseq"
            fastp --in1 ${forward} --in2 ${reverse} \\
                    --out1 mtbc_reads/${sampleID}_mtbc_R1.fastq.gz \\
                    --out2 mtbc_reads/${sampleID}_mtbc_R2.fastq.gz \\
                    --reads_to_process 5000000 --length_required 50
        else
            cp ${forward} mtbc_reads/${sampleID}_mtbc_R1.fastq.gz
            cp ${reverse} mtbc_reads/${sampleID}_mtbc_R2.fastq.gz
        fi

    """


}