
process SYLPH_CLASSIFICATION {

    conda params.readQC_env

    publishDir "${params.outDir}/db/qc/${params.runID}/"

    input:
        path(samplesheet)
        tuple path(sylph_db),
            path(sylph_tax)

    output:
        tuple path("${params.runID}.sylph_sequence_abundance.tsv"), 
            path("${params.runID}.sylph_relative_abundance.tsv"),
            path("${params.runID}.sylph_coverage.tsv"), emit: sylph_out

    script:

    """
    mkdir -p sylph/ sylph_tax/

    # Example: iterate through all samples in the samplesheet
    while IFS="," read -r origID sampleID forward reverse type; do

        # skip rows with missing FASTQs
        [[ -z \$forward || -z \$reverse ]] && continue

        sylph sketch \
            --out-name-db \$sampleID \
            -1 \$forward \
            -2 \$reverse \
            -d sylph/ \
            -t ${task.cpus}

    done < <(sed '1d' ${samplesheet})

    # Profile the sketches with Sylph
        sylph profile \\
            ${sylph_db} \\
            sylph/* \\
            --estimate-unknown \\
            --read-seq-id 0.95 \\
            --min-count-correct 3 \\
            --min-number-kmers 50 \\
            -t ${task.cpus} \\
            -o sylph.tsv

    # Get taxonomy for the profiles
        sylph-tax taxprof sylph.tsv \\
            -t ${sylph_tax}/* \\
            -o sylph_tax/tax_
            
    # remove any empty files
            find sylph_tax/ -type f -name 'tax_*.sylphmpa' -empty -delete

        sylph-tax merge \\
            sylph_tax/tax_*.sylphmpa \\
            --column sequence_abundance \\
            -o ${params.runID}.sylph_sequence_abundance.tsv

        sylph-tax merge \\
            sylph_tax/tax_*.sylphmpa \\
            --column relative_abundance \\
            -o ${params.runID}.sylph_relative_abundance.tsv
        
        sylph-tax merge \\
            sylph_tax/tax_*.sylphmpa \\
            --column relative_abundance \\
            -o ${params.runID}.sylph_coverage.tsv
    """

}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.2.0
    @description: 
        Use sylph to classify reads from a sample.
        This process takes the sample ID and the Sylph database as input,
        and outputs a merged Sylph sequence abundance file.
    @changelog
        v1.0.0-2025-04-01: Initial version
        v.1.1.0-2025-11-17: Modified for intergration with read QC report,
                            and for workflow download of database (GTRDB-R220 only)
        v.1.2.0-2025-11-26: Added sylph coverage and relative_abundance output.
                            Added sylph profile parameters for better accuracy 
                                included: --min-count-correct 3 (default)
                                            --min-number-kmers 50 (default)
        v1.3.0-2025-12-19: Modified to handle when the samlesheet has an empty line
*/