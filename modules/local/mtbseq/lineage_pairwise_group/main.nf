process MTBSEQ_LINEAGE_PAIRWISE_GROUP {

    tag "${lineage}_${snp_distance}"

    conda { file("/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq").exists() ? "/imppc/labs/emlab/phesketh/miniconda3/envs/mtbseq" : "./modules/local/mtbseq/mtbseq.yml" }

    container "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data"

    publishDir "${params.outdir}/bbdd/mtbseq/lineages/${lineage}", mode: 'copy'

    input:
        val(lineage)
        each snp_distance
        path(amended_dir)
        path(join_dir)
        path(samples_file)

    output:
        // Groups outputs
        tuple val(lineage), val(snp_distance), path("Groups/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u*_phylo_w*_d*.groups"), emit: groups_groups
        tuple val(lineage), val(snp_distance), path("Groups/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u*_phylo_w*_d*.list"), emit: groups_list
        tuple val(lineage), val(snp_distance), path("Groups/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u*_phylo_w*.matrix"), emit: groups_matrix 

    script:
    def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

    """
    MTBseq --step TBgroup \\
        --thread ${task.cpus} \\
        --prefix ${lineage} \\
        --sample ${samples_file} \\
        ${additional_args} \\
        --distance ${snp_distance} \\
        1>>.command.out \\
        2>>.command.err \\
        || true # NOTE This is a hack to overcome the exit status 1 thrown by mtbseq
    """
}