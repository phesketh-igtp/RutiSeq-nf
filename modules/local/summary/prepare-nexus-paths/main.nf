process PREPARE_NEXUS_PATHS{

    conda params.snp_profiling_env 

    tag "${lineage}"

    publishDir "${params.outdir}/bbdd/results/networks/${lineage}/", mode: 'copy'

    input:
        path pairwise_clusters
        tuple val(lineage), 
                path(joint_dir), 
                path(amend_dir)

    output:
        path("nexus.tuple.csv"),        emit: nexus_tuple

    script:

    """
        # Identify the unique clusters
        cut -f1 ${pairwise_clusters} | sort | uniq > unique.clusters.list

        for clusterID in \${cat unique.clusters.list}; do
            echo "${lineage},\${clusterID},${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.tab" >> nexus.tuple.csv
        done
    """
    
}