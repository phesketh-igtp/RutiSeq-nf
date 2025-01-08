process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda params.mtbseq_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)
        
        tuple val(sampleID), path(forward), path(reverse), path(mtbseq_class), 
                path(mtbseq_stats), path(mtbseq_pos), path(mtbseq_vars), 
                path(tbdb_out), path(who_out), path(mtbseq_vcf)
                
    output:
        tuple val(sampleID), path("Bam/${sampleID}.bam")
        tuple val(sampleID), path("Bam/${sampleID}.bam.bai")
        tuple val(sampleID), path("Bam/${sampleID}.bamlog")
        tuple val(sampleID), path("Called/${sampleID}.gatk_position_uncovered_*.tab")
        tuple val(sampleID), path("Called/${sampleID}.gatk_position_variants_*.tab")
        tuple val(sampleID), path("Classification/Strain_Classification.tab")
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bam")
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bai")
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bamlog")
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.grp")
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.intervals")
        tuple val(sampleID), path("Mpileup/${sampleID}.gatk.mpileup"),                  emit: mtbseq_mpileup
        tuple val(sampleID), path("Mpileup/${sampleID}.gatk.mpileuplog")
        tuple val(sampleID), path("Position_Tables/${sampleID}.gatk_position_table.tab")
        tuple val(sampleID), path("Statistics/Mapping_and_Variant_Statistics.tab")

        // tuple for updating the sample ch
        tuple val(sampleID), path(forward), path(reverse), 
                path("Classification/Strain_Classification.tab"), 
                path("Statistics/Mapping_and_Variant_Statistics.tab"), 
                path("Position_Tables/${sampleID}.gatk_position_table.tab"), 
                path("Called/${sampleID}.gatk_position_variants_*.tab"),  
                path(tbdb_out), path(who_out), 
                path(mtbseq_vcf),                                                       emit: updated_sample_ch

    script:
    
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # remove the default symbolic links it does to prevent mtbseq using the reads
    unlink ${forward} 
    unlink ${reverse}

    # Run MTBseq for a single sample
    MTBseq --step TBfull \\
        --thread ${task.cpus} \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    ## --prefix ${sampleID} \\

    """
    
}