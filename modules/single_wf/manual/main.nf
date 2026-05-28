process SINGLE_GENOME {


    tag "$sampleID"

    conda params.manual_tb_env

    input:
        tuple val(sampleID), 
            path(snippy_bam),
            path(snippy_bam_bai)

    output:
        // Emit ch for the updated channel with all the outputs
        tuple val(sampleID), 
            path("${sampleID}.vcf"),
            path("${sampleID}.LC_position.txt"), emit: lowFreq_ch

    script:

    """
    # Mpileup 
    samtools mpileup \\
        -q 20 \\
        -f ${params.snippy_reference} \\
        ${snippy_bam} \\
        -o ${sampleID}.mpileup

    # Low coverage positions 
    ## (from https://gitlab.com/LPCDRP/illumina-blindspots.pub, Modlin et al., 2020 Microbial Genomics, 
    ##   doi:10.1099/mgen.0.000465)
    python ${params.scriptDir}/LC_positions.py --mpileup ${sampleID}.mpileup

    # Variant calling (VarScan) 
    java -jar VarScan mpileup2cns "${sampleID}.mpileup" \\
        --min-avg-qual 20 \\
        --min-coverage 10 \\
        --variants \\
        --output-vcf 1 \\
        --strand-filter 0 \\
        > ${sampleID}.vcf

    # 
    """
}