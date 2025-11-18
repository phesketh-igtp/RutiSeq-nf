process SNIPPY_CORE {

    conda params.snippy_env

    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        val(sampleID_list)

    output:
        path("core.snps")
        path("core.vcf.bgz")
        path("core.aln"), emit: snippy_out_ch
        path("core.aln.full")

    script:

    """
    # Collect all the paths for snippy core
        paths=\$(echo ${params.outDir}/db/samples/*/snippy)

    # Run snippy core
        snippy-core \\
            --ref ${params.snippy_reference} \\
            \$(echo ${params.outDir}/db/samples/*/snippy)

    # Housekeeping
        bgzip -9 core.vcf
        tabix -t core.vcf.bgz
    """
}