
process SYLPH_CLASSIFICATION {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.1.0
    @description: 
        Use sylph to classify reads from a sample.
        This process takes the sample ID and the Sylph database as input,
        and outputs a merged Sylph sequence abundance file.
    @changelog
        v1.0.0-2025-04-01: Initial version
        v.1.1.0-2025-11-17: Modified for intergration with read QC report,
                            and for workflow download of database (GTRDB-R220 only)
*/

    conda params.readQC_env

    publishDir "${params.outDir}/db/qc/${params.runID}/", mode: 'copy'

    input:
        path(samplesheet)
        tuple path(sylph_db),
            path(sylph_tax)

    output:
        path("${params.runID}.sylph_sequence_abundance_file.tsv"), emit: sylph_out

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

        done < <(sed '1d' ${samplesheet})


    # Profile the sketches with Sylph
        sylph profile \\
            ${sylph_db} \\
            sylph/* \\
            --estimate-unknown \\
            --read-seq-id 0.98 \\
            -t ${task.cpus} \\
            -o sylph.tsv

    # Get taxonomy for the profiles
        sylph-tax taxprof sylph.tsv \\
            -t ${sylph_tax}/* \\
            -o sylph/tax_
            
    # remove any empty files
            find sylph/ -type f -name 'tax_*.sylphmpa' -empty -delete

        sylph-tax merge \\
            sylph/tax_*.sylphmpa \\
            --column sequence_abundance \\
            -o ${params.runID}.sylph_sequence_abundance_file.tsv
    """


}