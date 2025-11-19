process SNIPPY_CORE {
    conda "bioconda::snippy=4.6.0"
    
    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        val(sampleID_list)  // Collection of VCF files from snippy runs
        path(vcf_files_ch)  // Reference genome file

    output:
        path("core.snps"), emit: snps
        path("core.vcf.bgz"), emit: vcf_bgz
        path("core.vcf.bgz.tbi"), emit: vcf_index
        path("core.aln"), emit: alignment
        path("core.aln.full"), emit: alignment_full
        path("core.tab")
        path("core.txt")

    script:

    def vcf_dirs = vcf_files_ch.collect { it.parent }.join(' ')

    """
    # Run snippy core with collected VCF directories
    snippy-core --ref ${params.snippy_reference} ${vcf_dirs}

    # Compress and index the VCF file
    bgzip -c core.vcf > core.vcf.bgz
    tabix -p vcf core.vcf.bgz
    """
}