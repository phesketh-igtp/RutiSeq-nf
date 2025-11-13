process SNIPPY_CORE {

    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        path(mapping_stats)

    output:
        path("core.snps")
        path("core.vcf")
        path("core.aln")
        path("core.aln.full")

    script:

    """
    # Collect all the data from the snippy outputs
    for directory in "${params.outDir}/db/sample/*"; do
        sampleID=snippyDir_$(basename \$directory)
        mkdir -f \$directory
        ln -s \$directory/* \$sampleID/
    done

    paths=$(echo snippyDir_*)

    snippy-core \\
        --ref ${params.snippy_reference} \\
        \$paths

    # Housekeeping
    bgzip -9 core.vcf
    tabix -t core.vcf.bgz
    """
}