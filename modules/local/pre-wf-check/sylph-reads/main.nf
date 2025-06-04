process SYLPH_READ_CLASSIFICATION {

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
    # Create directory
        mkdir -p reads sylph
    
    # Parse the reads in the samplesheet
        while IFS=',' read -r originalID sampleID forward reverse type; do

            # Profile the sketches with Sylph
            sylph sketch \\
                -1 \${forward} \\
                -2 \${reverse}\\
                -d sylph/\${sampleID} \\
                -t ${task.cpus}

        done < ${params.samplesheet}

    # Profile the sketches
        sylph profile \\
                ${params.sylph_bbdd} \\
                sylph/*syldb \\
                --estimate-unknown \\
                --read-seq-id 0.99 \\
                -t ${task.cpus} \\
                -o sylph_bacterial.tsv

    # Get taxonomy for the profiles
        sylph-tax download --download-to my_existing_folder/
        sylph-tax taxprof sylph_bacterial.tsv \\
            -t ${params.sylph_bbdd_id} \\
            -o sylph/tax_

    # remove any empty files
        find sylph/ -type f -name 'tax_*.sylphmpa' -empty -delete

    # Merge the results
        sylph-tax merge \\
            sylph/tax_*.sylphmpa \\
            --column sequence_abundance \\
            -o merged_sylph_sequence_abundance_file.tsv
    """

}