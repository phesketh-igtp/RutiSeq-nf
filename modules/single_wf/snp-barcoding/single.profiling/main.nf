process SNP_PROFILING_SINGLE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        This process runs the SNP profiling step using the MTBSeq Mpileup using VarScan2.
        It takes the mpileup output from the MTBSeq pipeline and generates a VCF file.
        The VCF file is then compressed and indexed using bgzip and tabix.
        The process also removes the published reads from the previous module to recover storage.
        The outputs are published to a specified directory.
*/

    tag "$sampleID"

    conda params.snp_profiling_env
    
    container { 
            if (workflow.containerEngine == 'singularity') return params.singularity_snp_profiling
            else if (workflow.containerEngine == 'docker') return params.docker_snp_profiling
            else if (workflow.containerEngine == 'apptainer') return params.apptainer_snp_profiling
            else return null
        }
    
    publishDir "${params.outDir}/db/samples/${sampleID}/snippy/", mode: 'copy'

    input:
        tuple val(sampleID), 
                path(fastq_1), path(fastq_2), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf), path(mtbseq_mpileup)

    output:
        // tuple for updating the sample ch
        tuple val(sampleID), path(fastq_1), path(fastq_2), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars),  
                path(tbdb_out), path(who_out), 
                path("${sampleID}.vcf"),                     emit: updated_sample_ch5

        path("${sampleID}.aligned.fa")
        path("${sampleID}.consensus.fa.gz")
        path("${sampleID}.vcf")

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # Check if R2 file has actual sequencing data
        if [[ \$(zcat ${fastq_2} | head -n 4 | wc -l) -eq 4 ]]; then
            echo "Using paired-end mode for ${sampleID}"
            
            snippy \\
                --sampleID ${sampleID} \\
                --outdir . \\
                --reference ${params.snippy_reference} \\
                --R1 ${fastq_1} --R2 ${fastq_2} \\
                --force
        else
            echo "Using single-end mode for ${sampleID}"
            
            snippy \\
                --sampleID ${sampleID} \\
                --outdir . \\
                --reference ${params.snippy_reference} \\
                --se ${fastq_1} \\
                --force
        fi

        gzip --best ${sampleID}.consensus.fa
    """

}