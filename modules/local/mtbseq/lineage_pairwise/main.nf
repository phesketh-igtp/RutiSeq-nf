process MTBSEQ_LINEAGE_PAIRWISE {

    tag "${runID}_${lineage}"

    conda params.mtbseq_env

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
                
    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/${lineage}", mode: 'copy'

    input:
        val runID
        tuple val(lineage), path(samples_file)

    output:
        // Amend outputs
        tuple val(lineage), path("Amend/"),                                                  emit: amend_dir
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended.tab")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.fasta")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.plainIDs.fasta")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.tab")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.fasta")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.plainIDs.fasta")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*_removed.tab")
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.tab")
        // Join outputs
        tuple val(lineage), path("Join/"),                                                   emit: join_dir
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log")
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab")

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """

    mkdir Position_Tables/ Called/

    sed 's@\t@_@g' ${samples_file} samplesID.list

    while IFS=',' read -r sampleID; do
        ln -s ${params.outdir}/bbdd/mtbseq/pairwise/\${sampleID}/Position_Tables/* Position_Tables/
        ln -s ${params.outdir}/bbdd/mtbseq/pairwise/\${sampleID}/Called/* Called/
    done

    MTBseq --step TBjoin \\
        --thread ${task.cpus} \\
        --prefix ${lineage} \\
        --sample ${samples_file} \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    MTBseq --step TBamend \\
        --thread ${task.cpus} \\
        --prefix ${lineage} \\
        --samples ${samples_file} \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    """
}