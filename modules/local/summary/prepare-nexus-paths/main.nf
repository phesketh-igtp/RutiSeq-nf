process PREPARE_NEXUS_PATHS{

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process prepares the paths for the NEXUS files for each cluster.
        It generates a CSV file with the paths to the NEXUS files and the
        corresponding tab files. The CSV file is used as input for the
        GENERATE_NEXUS process.
*/

    conda params.snp_profiling_env 

    tag "${lineage}"

    publishDir "${params.outDir}/results/networks/${lineage}/", mode: 'copy'

    input:
        tuple val(lineage), 
                path(contree), 
                path(alignments)
        path pairwise_clusters


    output:
        path("nexus.tuple.csv"), emit: nexus_tuple

    script:

    """
    # Identify the unique clusters
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f7 \\
            | sort \\
            | uniq \\
            | sed '1!{/^SampleID/d;}' \\
            | sed '1!{/^nX-/d;}' > unique.clusters.list

    # Remove clusters that are smaller than 3 genomes
        while read clusterID; do
            count=\$(grep -c "\${clusterID}" "${pairwise_clusters}")
            if [ "\${count}" -ge 3 ]; then
                echo "\${clusterID}" >> final_clusters.list
            fi
        done < unique.clusters.list

    # Create a CSV for generating a tuple of the paths
        for clusterID in `cat final_clusters.list`; do
            grep "\${clusterID}" ${pairwise_clusters} | cut -f1 | tr '\n' ';' > tmp-string
            echo "${lineage},\${clusterID},${params.outDir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outDir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab" >> nexus.tuple.csv
            rm tmp-string
        done

    # touch the output incase the file is empty
        touch nexus.tuple.csv
    """

}