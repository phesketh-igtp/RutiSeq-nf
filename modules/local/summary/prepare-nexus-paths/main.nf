process PREPARE_NEXUS_PATHS{

    conda params.snp_profiling_env 

    tag "${lineage}"

    array 100

    publishDir "${params.outdir}/results/networks/${lineage}/", mode: 'copy'

    input:
        tuple val(lineage), 
            path(contree), 
            path(alignments)
        path pairwise_clusters


    output:
        path("nexus.tuple.csv"),    emit: nexus_tuple

    script:

    """
    # Identify the unique clusters
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f7 \\
            | sort \\
            | uniq \\
            | sed '1!{/^SampleID/d;}' \\
            | sed '1!{/^nX-/d;}' > unique.clusters.list

    # Remove clusters that are smaller than 5 genomes
        while read clusterID; do
            count=\$(grep -c "\${clusterID}" "${pairwise_clusters}")
            if [ "\${count}" -ge 5 ]; then
                echo "\${clusterID}" >> final_clusters.list
            fi
        done < unique.clusters.list

    # Create a CSV for generating a tuple of the paths
        for clusterID in `cat final_clusters.list`; do
            echo "${lineage},\${clusterID},${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.tab" >> nexus.tuple.csv
        done

    # touch the output incase the file is empty
        touch nexus.tuple.csv
    """

}