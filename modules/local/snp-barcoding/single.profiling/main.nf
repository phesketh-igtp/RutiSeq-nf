process SNP_PROFILING_SINGLE {

    tag "$sampleID"

    conda './envs/conda/bcftools-env.yml'
    
    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e1dc9586886df729616e7af235efe76cc2d31b5fa2a6afe0b0656efee6d983a/data'
            } else { 'community.wave.seqera.io/library/bcftools_snpeff_varscan_vcftools:3fa84761d1a9bed3' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'copy'

    input:
    tuple val(sampleID), path(mtbseq_mpileup)

    output:
        tuple val(sampleID), path("${sampleID}.gatk.vcf.gz"),       emit: mtbseq_vcf
        tuple val(sampleID), path("${sampleID}.gatk.vcf.gz.tbi"),   emit: mtbseq_vcf_index

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

}