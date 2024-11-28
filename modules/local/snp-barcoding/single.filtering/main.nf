process SNP_FILTERING_SINGLE {

    tag "$sampleID"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/snp-profiling").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/snp-profiling" : "./modules/local/snp-barcoding/snp-profiling.yml" }

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f7/f7b4d880135db1027c7b0ead23ae6464da73fc498f84b2462ff36af4ac0f198b/data'

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'copy'

    input:
    tuple val(sampleID), path(mtbseq_vcf), path(mtbseq_vcf_index)

    output:
        tuple val(sampleID), path("${sampleID}.gatk.vcf.gz"), emit: mtbseq_vcf_annot

    script:

    // defined in the nextflow.config file
    def additional_args_snpeff = task.ext.additional_args_snpeff ?: '' // SnpEff arguments
    def additional_args_filt_snps = task.ext.additional_args_filt_snps ?: '' // filtering SNPs by strings

    """
    
    java -jar snpEff/snpEff.jar \\
        ${additional_args_snpeff} \\
        MTB_ancestor \\
        ${mtbseq_vcf} \\
        > ${sampleID}.gatk.annot.vcf.gz

    """

    stub:
    """

    """
}