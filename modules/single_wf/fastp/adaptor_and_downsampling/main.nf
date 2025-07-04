process ADAPTORS_AND_DOWNSAMPLING {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-07-02
    @version: 1.1.0
    @description: 
        USe gfo fastp to remove any possible Illumina adapters from the reads and
        downsample to 5,000,000 reads for TBProfiler/MTBseq (if necessary).
        This will also halt the analysis of genomes with less than 500,000 reads, since
        these will likely fail the TBProfiler/MTBseq steps anyway.
    @changelog
        v1.0.0-2025-04-01: Initial version
        v1.0.1-2025-04-04: Added documentation and comments
        v1.1.0-2025-07-02: Updated half analysis of genomes with less than 500,000 reads
*/

    tag "$sampleID"

    conda params.taxonomy_env

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_kaiju
            else if (workflow.containerEngine == 'docker') return params.docker_kaiju
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_kaiju
            else return null
        }

    input:
        tuple val(sampleID), 
            path(forward), path(reverse), path(mtbseq_class), 
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf)

    output:
        // Emit ch for the updated channel with all the outputs
        tuple val(sampleID), 
            path("fastp/${sampleID}_R1.fastq.gz",), 
            path("fastp/${sampleID}_R2.fastq.gz",), 
            path(mtbseq_class), path(mtbseq_stats), 
            path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf), emit: updated_sample_ch1
    
        //path("failed_sample_entry.txt"), emit: failed_sample_entry

    script:

    """
    # Remove any possible illumina adapters from the reads and
    ## downsample to 5,000,000 reads for TBProfiler/MTBseq (if necessary)
        mkdir -p fastp
        fastp \\
            --in1 ${forward} \\
            --in2 ${reverse} \\
            --out1 fastp/${sampleID}_R1.fastq.gz \\
            --out2 fastp/${sampleID}_R2.fastq.gz \\
            --detect_adapter_for_pe \\
            --overrepresentation_analysis \\
            --reads_to_process ${params.fastp_max_reads} \\
            --length_required ${params.fastp_length_required}
    """
}