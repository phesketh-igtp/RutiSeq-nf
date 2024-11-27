process MTBSEQ_SINGLE {

    tag "$sampleID"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq" : "./modules/local/mtbseq/mtbseq.yml" }

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data"

    publishDir "${params.outdir}/bbdd/mtbseq/${sampleID}", mode: 'link'

    input:
        tuple val(sampleID), path(mtbc_forward), path(mtbc_reverse)

    output:
        path "Bam/", emit: bam_dir
        path "Bam/${sampleID}.bam", emit: bam
        path "Bam/${sampleID}.bam.bai", emit: bam_index
        path "Bam/${sampleID}.bamlog", emit: bamlog
        path "Position_Tables", emit: position_tables_dir
        path "Position_Tables/${sampleID}.gatk_position_table.tab", emit: position_tables
        path "Classification", emit: classification_dir
        path "Classification/Strain_Classification.tab", emit: classification
        path "Statistics", emit: statistics_dir
        path "Statistics/Mapping_and_Variant_Statistics.tab", emit: statistics
        path "Called/", emit: called_dir
        path "Called/*gatk_position_variants*.tab", emit: position_variants
        path "Mpileup/", emit: mpileup_dir
        path "Mpileup/${sampleID}*.gatk.mpileup", emit: mpileup

    script:
    
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """

    # Run MTBseq for a single sample
    MTBseq --step TBfull \\
        --thread ${task.cpus} \\
        --prefix ${sampleID} \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """

    stub:
    """
    mkdir -p Amend Position_Tables Classification Statistics Called
    touch Statistics/Mapping_and_Variant_Statistics.tab
    touch Classification/Strain_Classification.tab
    touch Called/gatk_position_variants.tab
    touch Position_Tables/gatk_position_table.tab
    """
}