
process SYLPH_CLASSIFICATION {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        Use sylph to classify reads from a sample.
        This process takes the sample ID and the Sylph database as input,
        and outputs a merged Sylph sequence abundance file.
    @changelog
        v1.0.0-2025-04-01: Initial version
*/

    conda params.taxonomy_env

    publishDir "${params.outDir}/bbdd/read-qc/", mode: 'copy'

    input:
        val(runID)

    output:
        path("${params.runID}.merged_sylph_sequence_abundance_file.tsv"), emit: sylph_res

    script:

    """
    mkdir -p sylph
    
    # Example: iterate through all samples in the samplesheet
        while IFS="," read -r origID sampleID forward reverse type; do

            ln -s "\${forward}" "\${sampleID}_R1.fastq.gz"
            ln -s "\${reverse}" "\${sampleID}_R2.fastq.gz"
            sylph sketch \\
                -1 "\${sampleID}_R1.fastq.gz" \\
                -2 "\${sampleID}_R2.fastq.gz" \\
                -d sylph/ \\
                -t ${task.cpus}

        done < <(sed '1d' ${params.samplesheet})


    # Profile the sketches with Sylph
        sylph profile \\
            ${params.sylph_bbdd} \\
            sylph/* \\
            --estimate-unknown \\
            --read-seq-id 0.99 \\
            -t ${task.cpus} \\
            -o sylph.tsv

    # Get taxonomy for the profiles
        mkdir -p my_existing_folder/
        sylph-tax download --download-to my_existing_folder/
        sylph-tax taxprof sylph.tsv \\
            -t ${params.sylph_bbdd_id} \\
            -o sylph/tax_
            
    # remove any empty files
            find sylph/ -type f -name 'tax_*.sylphmpa' -empty -delete
        sylph-tax merge \\
            sylph/tax_*.sylphmpa \\
            --column sequence_abundance \\
            -o ${params.runID}.merged_sylph_sequence_abundance_file.tsv
    """


}