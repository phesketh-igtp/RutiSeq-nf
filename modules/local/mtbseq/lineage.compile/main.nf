process MTBSEQ_LINEAGE_COMPILE {
    tag "cohort"
    label 'process_tbfull'
    label 'error_retry'

    conda "bioconda::mtbseq=1.1.0"

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }

    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${runID}", mode: 'copy'

    input:
        tuple val(old_name), val(sampleID), path(forward), path(reverse)

    output:
        path "Amend/",                                                  emit: amend_dir 
        path "Position_Tables/",                                        emit: position_tables_dir
        path "Classification/",                                         emit: classification_dir
        path "Statistics/",                                             emit: statistics_dir
        path "Statistics/Mapping_and_Variant_Statistics.tab",           emit: statistics
        path "Classification/Strain_Classification.tab",                emit: classification
        path "Called/*gatk_position_variants*.tab",                     emit: position_variants
        path "Position_Tables/*.gatk_position_table.tab",               emit: position_tables

    script:
        def args = task.ext.args ?: " --minbqual ${params.minbqual} --mincovf ${params.mincovf} --mincovr ${params.mincovr} --minphred ${params.minphred} --minfreq ${params.minfreq} --unambig ${params.unambig} --window ${params.window}"

        """

        """

}