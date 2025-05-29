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

    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_kaiju
            else if (workflow.containerEngine == 'docker') return params.docker_kaiju
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_kaiju
            else return null
        }

    publishDir "${params.outDir}/bbdd/read-qc/${params.runID}/", mode: 'copy'

    input:
        val all_samples

    output:
        path("merged_sylph_sequence_abundance_file.tsv")

    script:

    """
    # Access individual samples like:
    # all_samples[0] = first tuple [sampleID, forward, reverse, mtbseq_class, ...]
    # all_samples[1] = second tuple [sampleID, forward, reverse, mtbseq_class, ...]
    
        echo "Processing ${all_samples.size()} samples"

        mkdir -p reads sylph
    
        # Example: iterate through all samples

            for sample in ${all_samples}:
                echo "Sample ID: \${sample[0]}"
                echo "Forward reads: \${sample[1]}"
                echo "Reverse reads: \${sample[2]}"
                
                ln -s \${sample[1]}" reads/\${sample[0]}_R1.fastq.gz"
                ln -s \${sample[2]}" reads/\${sample[0]}_R2.fastq.gz"

                sylph sketch \\
                    -1 \${sample[1]} \\
                    -2 \${sample[2]}\\
                    -d sylph/\${sample[0]} \\
                    -t ${task.cpus}
            done

    # Profile the sketches with Sylph
        sylph profile \\
            ${params.sylph_bacterial_bbdd} \\
            sylph/*syldb \\
            --estimate-unknown \\
            --read-seq-id 0.99 \\
            -t ${task.cpus} \\
            -o sylph_bacterial.tsv

    # Get taxonomy for the profiles
    sylph-tax download --download-to my_existing_folder/

    sylph-tax taxprof sylph_bacterial.tsv \\
        -t ${params.sylph_bacterial_bbdd_id} \\
        -o sylph/tax_bacterial_

    # remove any empty files
        find sylph/ -type f -name 'tax_bacterial_*.sylphmpa' -empty -delete

    sylph-tax merge \\
        sylph/tax_*.sylphmpa \\
        --column sequence_abundance \\
        -o merged_sylph_sequence_abundance_file.tsv
    """


}