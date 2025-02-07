process SNP_PROFILING_SINGLE {

    tag "$sampleID"

    conda params.snp_profiling_env
    
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
            else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
            else return null
        }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/SNP-Profiles/", mode: 'copy'

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


    # remove the published reads from the previous module:
        for file in \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${sampleID}_mtbc_R1.fastq.gz"" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/${sampleID}_mtbc_R2.fastq.gz"" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/tbdb-${sampleID}.results.txt" \\
            "${params.outdir}/bbdd/mtbseq/samples/${sampleID}/who-${sampleID}.results.txt" \\
        do
            if [ -f "\${file}" ] || [ -e "\${file}" ]; then
                rm "\${file}"
            fi
        done

    # Compress the outputs from MTBSeq mpileup
        gzip --force --best ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileup
        gzip --force --best ${params.outdir}/bbdd/mtbseq/samples/${sampleID}/Mpileup/${sampleID}.gatk.mpileuplog

    """

}