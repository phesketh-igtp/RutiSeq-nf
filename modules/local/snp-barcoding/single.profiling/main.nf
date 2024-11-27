process SNP_PROFILING_SINGLE {

    tag "$sampleID"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/snp-profiling").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/snp-profiling" : "./modules/local/snp-barcoding/snp-profiling.yml" }

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/608c097132a7de8e156c452f40ea3b3fea6bf0a35b6988e4b2fe74d91524303f/data'

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