process ADAPTORS_AND_DOWNSAMPLING {

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

    input:
        tuple val(sampleID), 
            path(forward), path(reverse), path(mtbseq_class), 
            path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf)

    output:
        // Emit ch for the updated channel with all the outputs
        tuple val(sampleID), 
            path("fastp/${sampleID}_R1.fastq.gz"), 
            path("fastp/${sampleID}_R2.fastq.gz"), 
            path(mtbseq_class), path(mtbseq_stats), 
            path(mtbseq_pos), path(mtbseq_vars), 
            path(tbdb_out), path(who_out), path(mtbseq_vcf), emit: updated_sample_ch1

    script:

    """
    # Remove any possible illumina adapters from the reads and
    ## downsample to 5,000,000 reads for TBProfiler/MTBseq (if necessary)

        mkdir -p fastp

        fastp --in1 ${forward} --in2 ${reverse} \\
            --out1 fastp/${sampleID}_R1.fastq.gz \\
            --out2 fastp/${sampleID}_R2.fastq.gz \\
            --reads_to_process 5000000 \\
            --length_required 50
    """
}