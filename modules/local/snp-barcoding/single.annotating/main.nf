process SNP_ANNOTATING_SINGLE {

    tag "$sampleID"

    conda 'bioconda::bcftools==1.21 bioconda::varscan==2.4.6 bioconda::snpeff==5.2'

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e1dc9586886df729616e7af235efe76cc2d31b5fa2a6afe0b0656efee6d983a/data'
            } else { 'community.wave.seqera.io/library/bcftools_snpeff_varscan_vcftools:3fa84761d1a9bed3' }
    }

    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbseq_vcf), path(mtbseq_vcf_index)
        path snpeff_ref_genome

    output:
        tuple val(sampleID), path("${sampleID}.gatk.annot.vcf.gz"), emit: mtbseq_vcf_annot
        tuple val(sampleID), path("${sampleID}.gatk.annot.vcf.gz.tbi"), emit: mtbseq_vcf_annot_index

    script:

    // defined in the nextflow.config file
    def additional_str_snpeff = task.ext.additional_str_snpeff ?: '' // SnpEff arguments

    """

    #Need to change the name of the vcf 'M.tuberculosis_H37Rv' to 'Mycobacterium_tuberculosis_h37rv'
    bgzip -d ${mtbseq_vcf} > ${sampleID}.gatk.vcf
    sed -i 's@${snpeff_ref_genome}@Chromosome@g' ${sampleID}.gatk.vcf
    bgzip ${sampleID}.gatk.vcf # recompress for annotation
    
    # Annotate SNPs
    snpEff ann \\
        ${task.additional_args_snpeff} \\
        ${snpeff_ref_genome} \\
        ${mtbseq_vcf} \\
        > ${sampleID}.gatk.annot.vcf

    # Recompress and index SNPs
    bgzip ${sampleID}.gatk.annot.vcf
    tabix ${sampleID}.gatk.annot.vcf

    """

    stub:
    """
    touch stub.gatk.annot.vcf.gz
    """
}