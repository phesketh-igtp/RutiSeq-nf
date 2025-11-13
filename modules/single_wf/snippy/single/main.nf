process SNIPPY_SINGLE {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        This process runs the SNP profiling step using a reference GBK file.
        It produces files that can used for snippy-core processing
        To generate large phylogenetic trees. Produces a consensus sequence.
    @changelog:
        v1.0.0-2025-04-01: Initial version, using MTBSeq mpileup output
        v2.0.0-2025-06-10: Change from using MTBSeq outputs to using Snippy for SNP profiling
                            Added support for single-end reads
                            Extended available parameters for Snippy
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
                path("${sampleID}.vcf"),                     emit: updated_sample_ch4

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
                --outdir . --force \\
                --R1 ${fastq_1} --R2 ${fastq_2} \\
                --ref ${params.snippy_reference} \\
                --mincov ${params.snippy_mincov} \\
                --minfrac ${params.snippy_minfrac} \\
                --mapqual ${params.snippy_mapqual} \\
                --minqual ${params.snippy_minqual} \\
                --basequal ${params.snippy_basequal} \\
                ${params.snippy_args}  
        else
            echo "Using single-end mode for ${sampleID}"
            snippy \\
                --sampleID ${sampleID} \\
                --outdir . --force \\
                --se ${fastq_1} \\
                --ref ${params.snippy_reference} \\
                --mincov ${params.snippy_mincov} \\
                --minfrac ${params.snippy_minfrac} \\
                --mapqual ${params.snippy_mapqual} \\
                --minqual ${params.snippy_minqual} \\
                --basequal ${params.snippy_basequal} \\
                ${params.snippy_args}  
        fi

        gzip --best ${sampleID}.consensus.fa
    """

}