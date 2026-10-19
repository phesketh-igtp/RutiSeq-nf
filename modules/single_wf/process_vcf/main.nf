process PROCESS_VCF {

    tag "$sampleID"

    conda params.snippy_env
    
    container { 
        if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
        else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
        else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
        else return null
    }
    publishDir "${params.outDir}/db/samples/${sampleID}/snippy/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(sampleID), 
            path(snippy_vcf)

    output:
        path("${sampleID}.QC.mask.filt.snps.consensus.fa.gz")
        path("${sampleID}.QC.mask.filt.vcf.gz")
        path("${sampleID}.QC.mask.filt.vcf.gz.tbi")
        path("${sampleID}.processed.pass.vcf.gz")
        path("${sampleID}.processed.pass.vcf.gz.tbi")

    script:
    """
    #--------------- Quality filtering ---------------#
    bcftools view \\
            -e 'ALT~"^<" || REF~"N"' \\
            -Ou \\
            "${snippy_vcf}" |
        bcftools norm \\
            -a \\
            -f "${params.snippy_reference_fa}" \\
            -m -any \\
            -d exact \\
            -Ou | \\
        bcftools filter \
            -e 'QUAL<100 || INFO/DP<8 || INFO/MQM<30 || INFO/SAF<${params.mtbseq_mincovf} || INFO/SAR<${params.mtbseq_mincovr} || INFO/AO/(INFO/AO+INFO/RO)<0.${params.mtbseq_minfreq}' \
            --soft-filter low_QUAL \\
            -Ou | \\
        bcftools filter \\
            --mask-file "${params.snippy_mask}" \\
            --soft-filter masked \
            -Ou | \\
        bcftools filter \\
            --SnpGap "${params.vcf_indel_distance}":indel \\
            --soft-filter SnpGap \\
            -Ou | \\
        bcftools reheader \\ 
            --header "${sampleID}" \\
            -Oz | \\
            -o ${sampleID}.QC.mask.filt.vcf.gz

    tabix -f ${sampleID}.QC.mask.filt.vcf.gz

    #--------------- Apply only SNPs to consensus genome ---------------#
    bcftools view \\
        -f PASS \\
        -Oz \\
        -o ${sampleID}.QC.mask.filt.pass.vcf.gz \\
        ${sampleID}.QC.mask.filt.vcf.gz

    tabix -f ${sampleID}.QC.mask.filt.pass.vcf.gz

    bcftools view \\
        -f PASS \\
        -v snps \\
        ${sampleID}.QC.mask.filt.vcf.gz \\
        -Oz \\
        -o ${sampleID}.QC.mask.filt.snps.vcf.gz

    tabix -f ${sampleID}.QC.mask.filt.snps.vcf.gz

    #--------------- Concensus  ---------------#
    bcftools consensus \\
        -f ${params.snippy_reference_fa} \\
        ${sampleID}.QC.mask.filt.snps.vcf.gz \\
        > ${sampleID}.QC.mask.filt.snps.consensus.fa

    sed -i "s@>NC_000962.3@>${sampleID}@g" \\
        ${sampleID}.QC.mask.filt.snps.consensus.fa

    pigz -p 8 --best ${sampleID}.QC.mask.filt.snps.consensus.fa
    """

}

/*
@author: Poppy J Hesketh Best
@date: 2026-10-19
@version: 1.0
@description: 
    Filters the snippy vcf file
@changelog:
    v1.0.0-2025-04-01: Initial version, using MTBSeq mpileup output
    v2.0.0-2025-06-10: Change from using MTBSeq outputs to using Snippy for SNP profiling
                        Added support for single-end reads
                        Extended available parameters for Snippy
    v2.1.0-2026-01-19: Updated to check for existing results to avoid re-running
*/