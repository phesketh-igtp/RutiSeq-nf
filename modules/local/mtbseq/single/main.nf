process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda "./envs/conda/mtbseq-env.yml"

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
    
    publishDir "${params.outdir}/bbdd/mtbseq/samples/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)

    output:
        tuple val(sampleID), path("Bam/${sampleID}.bam"),                                   emit: mtbseq_bam
        tuple val(sampleID), path("Bam/${sampleID}.bam.bai"),                               emit: mtbseq_bam_index
        tuple val(sampleID), path("Bam/${sampleID}.bamlog"),                                emit: mtbseq_bamlog
        tuple val(sampleID), path("Called/${sampleID}.gatk_position_uncovered_*.tab"),      emit: mtbseq_uncovered_positions
        tuple val(sampleID), path("Called/${sampleID}.gatk_position_variants_*.tab"),       emit: mtbseq_vars
        tuple val(sampleID), path("Classification/Strain_Classification.tab"),              emit: mtbseq_class
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bam"),                         emit: mtbseq_gatk_bam
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bai"),                         emit: mtbseq_gatk_bam_index
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.bamlog"),                      emit: mtbseq_gatk_bamlog
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.grp"),                         emit: mtbseq_gatk_grp
        tuple val(sampleID), path("GATK_Bam/${sampleID}.gatk.intervals"),                   emit: mtbseq_gatk_intervals
        tuple val(sampleID), path("Mpileup/${sampleID}.gatk.mpileup"),                      emit: mtbseq_mpileup
        tuple val(sampleID), path("Mpileup/${sampleID}.gatk.mpileuplog"),                   emit: mtbseq_mpileuplog
        tuple val(sampleID), path("Position_Tables/${sampleID}.gatk_position_table.tab"),   emit: mtbseq_pos
        tuple val(sampleID), path("Statistics/Mapping_and_Variant_Statistics.tab"),         emit: mtbseq_stats

    script:
    
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """

    # Run MTBseq for a single sample
    MTBseq --step TBfull \\
        --thread ${task.cpus} \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    ## --prefix ${sampleID} \\

    """
    
}