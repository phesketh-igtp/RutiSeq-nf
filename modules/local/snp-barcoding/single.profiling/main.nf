process SNP_PROFILING_SINGLE {

    tag "$sampleID"

    conda params.snp_profiling_env
    
    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5e1dc9586886df729616e7af235efe76cc2d31b5fa2a6afe0b0656efee6d983a/data'
            } else { 'community.wave.seqera.io/library/bcftools_snpeff_varscan_vcftools:3fa84761d1a9bed3' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'move', pattern: '*.{gatk.vcf.gz,gatk.vcf.gz.tbi}'

    input:
        tuple val(sampleID), 
                path(forward), path(reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf), path(mtbseq_mpileup)

    output:
        tuple val(sampleID), 
                path("${sampleID}.gatk.vcf.gz"), 
                path("${sampleID}.gatk.vcf.gz.tbi"),                 emit: gatk_vcf_ch

        // tuple for updating the sample ch
        tuple val(sampleID), path(forward), path(reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars),  
                path(tbdb_out), path(who_out), 
                path("${sampleID}.gatk.vcf.gz"),                     emit: updated_sample_ch5

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