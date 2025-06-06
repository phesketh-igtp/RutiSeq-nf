process SYLPH_CLASSIFICATION {

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

    conda params.taxonomy_env

    publishDir "${params.outDir}/bbdd/read-qc/sylph/${params.runID}/", mode: 'copy'

    input:
        val(runID)

    output:
        path("merged_sylph_sequence_abundance_file.tsv")

    script:

    """
        mkdir -p reads sylph
    
        # Example: iterate through all samples in the samplesheet
            while IFS="," read -r origID sampleID forward reverse type; do
                
                ln -s \${forward} \${sampleID}_R1.fastq.gz
                ln -s \${reverse} \${sampleID}_R2.fastq.gz

                sylph sketch \\
                    -1 \${sampleID}_R1.fastq.gz \\
                    -2 \${sampleID}_R2.fastq.gz \\
                    -d sylph/\${sampleID} \\
                    -t ${task.cpus}

            done < ${params.samplesheet}

    # Profile the sketches with Sylph
        sylph profile \\
            ${params.sylph_bbdd} \\
            sylph/*syldb \\
            --estimate-unknown \\
            --read-seq-id 0.99 \\
            -t ${task.cpus} \\
            -o sylph.tsv

    # Get taxonomy for the profiles
    sylph-tax download --download-to my_existing_folder/

    sylph-tax taxprof sylph.tsv \\
        -t ${params.sylph_bbdd_id} \\
        -o sylph/tax_

    # remove any empty files
        find sylph/ -type f -name 'tax_*.sylphmpa' -empty -delete

    sylph-tax merge \\
        sylph/tax_*.sylphmpa \\
        --column sequence_abundance \\
        -o merged_sylph_sequence_abundance_file.tsv
    """


}