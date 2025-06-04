process FASTP_READ_QC {

    /*

/*
    @author: Poppy J Hesketh Best
    @date: 2025-06-04
    @version: 2.0.0
    @description: 
        This QC step runs fastp for adaptor removal and down-sampling of reads.
        It takes the forward and reverse reads, performs quality control, and outputs the cleaned reads.
        It also emits the updated sample channel with the cleaned reads and other outputs.
    @changelog
        v1.0.0-2025-04-01: Initial version
        v2.0.0-2025-06-04: Changed this module to run fastp for adaptor removal and down-sampling of reads.
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
                path("fastp/${sampleID}_R1.fastq.gz"), 
                path("fastp/${sampleID}_R2.fastq.gz"), 
                path(mtbseq_class), path(mtbseq_stats), 
                path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf),     emit: updated_sample_ch1

    script:

    """
    mkdir -p fastp

    fastp --in1 ${forward} --in2 ${reverse} \\
        --out1 fastp/${sampleID}_R1.fastq.gz \\
        --out2 fastp/${sampleID}_R2.fastq.gz \\
        --reads_to_process 5000000 --length_required 50
    """


}