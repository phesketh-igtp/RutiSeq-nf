process SNP_BARCODING_SINGLE {

    tag "$sampleID"

    conda './envs/conda/bcftools-env.yml'

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e1dc9586886df729616e7af235efe76cc2d31b5fa2a6afe0b0656efee6d983a/data'
            } else { 'community.wave.seqera.io/library/bcftools_snpeff_varscan_vcftools:3fa84761d1a9bed3' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/${sampleID}/SNP-Profiles", mode: 'link'

    input:
    tuple val(sampleID), path(mpileup)

    output:
    tuple val(sampleID), path("${sampleID}.gatk.vcf.gz"), emit: snp_barcoding_individual_vcf
    tuple val(sampleID), path("${sampleID}.gatk.vars"), emit: snp_barcoding_individual_vars
    tuple val(sampleID), path("${sampleID}.gatk.vcf.gz.tbi"), emit: snp_barcoding_individual_vcf_index

    script:
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    varscan mpileup2cns ${mpileup} \\
        ${additional_args} \\
        --output-vcf 1 \\
        --vcf-sample-list ${sampleID} \\
        > ${sampleID}.gatk.vcf

    varscan mpileup2cns ${mpileup} \\
        ${additional_args} \\
        > ${sampleID}.gatk.vars

    bgzip -c ${sampleID}.gatk.vcf > ${sampleID}.gatk.vcf.gz

    tabix -p vcf ${sampleID}.gatk.vcf.gz
    """

    stub:
    """
    touch ${sampleID}.gatk.vcf.gz
    touch ${sampleID}.gatk.vcf.gz.tbi
    touch ${sampleID}.gatk.vars
    """
}