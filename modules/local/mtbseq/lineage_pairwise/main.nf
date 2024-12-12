process MTBSEQ_LINEAGE_PAIRWISE {

    tag "${lineage}"

    conda "./envs/conda/mtbseq-env.yml"

    container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
            } else { 'quay.io/biocontainers/mtbseq' }
    }
                
    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${lineage}", mode: 'copy'

    input:
        val runID
        tuple val(lineage), val(sample_name), path(called_dir), path(position_tables_dir)

    output:
        // Amend outputs
        tuple val(lineage), path("Amend/"),                                                        emit: amended_dir
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended.tab"),                          emit: amended_tab
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.fasta"),               emit: phylo_fasta
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.plainIDs.fasta"),      emit: phylo_plainids_fasta
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo.tab"),                 emit: phylo_tab
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.fasta"),            emit: phylo_w_fasta
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.plainIDs.fasta"),   emit: phylo_w_plainids_fasta
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*_removed.tab"),      emit: phylo_w_removed_tab
        tuple val(lineage), path("Amend/${lineage}_joint_*_amended_u*_phylo_w*.tab"),              emit: phylo_w_tab
        // Join outputs
        tuple val(lineage), path("Join/"),                                                         emit: join_dir
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.log"),                                emit: join_log
        path("Join/${lineage}_joint_cf*_cr*_fr*_ph*_samples*.tab"),                                emit: join_tab

    script:

    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    # get the genomes that have >= "${params.mtbseq.minCov}"
    awk '{if (\$17 >= ${params.mtbseq.minCov}) print}' ${samples_file} | sed 's@_L@\tL@g' > ${lineage}.mtbseq.samples.txt

    MTBseq --step TBjoin \\
        --thread ${task.cpus} \\
        --prefix ${lineage} \\
        --sample ${lineage}.mtbseq.samples.txt \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq

    MTBseq --step TBamend \\
        --thread ${task.cpus} \\
        --prefix ${lineage} \\
        --samples ${lineage}.mtbseq.samples.txt \\
        ${additional_args} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq
    """
}