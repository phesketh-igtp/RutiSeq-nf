process SNP_PROFILING_SINGLE {

    tag "$sampleID"

    conda "../conda/snp-profiling.yml"

    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/608c097132a7de8e156c452f40ea3b3fea6bf0a35b6988e4b2fe74d91524303f/data'

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'copy'

    input:
    tuple val(sampleID), path(mtbseq_mpileup)

    output:
        tuple val(sampleID), path("${sampleID}.gatk.vcf.gz"), emit: mtbseq_vcf
        tuple val(sampleID), path("${sampleID}.gatk.vcf.gz.tbi"), emit: mtbseq_vcf_index

    script:
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    echo ${sampleID} > sample.list

    varscan mpileup2cns ${mtbseq_mpileup} \\
        ${additional_args} \\
        --output-vcf 1 \\
        --vcf-sample-list sample.list \\
        > ${sampleID}.gatk.vcf

    bgzip -c ${sampleID}.gatk.vcf > ${sampleID}.gatk.vcf.gz

    tabix -p vcf ${sampleID}.gatk.vcf.gz

    rm sample.list

    """

    stub:
    """
    touch ${sampleID}.gatk.vcf.gz
    touch ${sampleID}.gatk.vcf.gz.tbi
    """

}