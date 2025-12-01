process PREPARE_NEXUS_PATHS{

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0.2
    @description:
        This process prepares the paths for the NEXUS files for each cluster.
        It generates a CSV file with the paths to the NEXUS files and the
        corresponding tab files. The CSV file is used as input for the
        GENERATE_NEXUS process.
        The tuple has the following format:
        ["lineage", "clusterID", "fasta_path", "tab_path"]
    @changelog:
        v1.0.0-2025-04-01: Initial version
        v1.0.1-2025-04-04: Added filtering to remove clusters with less than 3 genomes
        v1.0.2-2025-12-01: Updated the paths for the new system.
*/

    conda params.snippy_env 

    tag "${lineage} (t=${distance})"

    publishDir "${params.outDir}/results/networks/${lineage}/", mode: 'copy'

    errorStrategy 'ignore'

    input:
        tuple val(lineage),
            val(distance),
            path(snp_fasta),
            path(snps_tab),
            path(clusters_tab)

    output:
        path("nexus.tuple.csv"), emit: nexus_tuple
        path(clusters_tab),      emit: clusters_tab

    script:

    def snp_fasta_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
    def snp_tab_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab"
    def clusters_tab_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/${lineage}_d${distance}.clusters.tsv"

    """
    lin=\$(echo "${lineage}" | sed 's/lineage/L/g')

    # Identify the unique clusters
        grep -w "${lineage}" ${clusters_tab} \\
            | cut -f4 \\
            | sort \\
            | uniq > unique.clusters.list

    # Remove clusters that are smaller than 3 genomes
        while read clusterID; do
            count=\$(grep -c "\${clusterID}" "${clusters_tab}")
            if [ "\${count}" -ge 3 ]; then
                echo "\${clusterID}" >> final_clusters.list
            fi
        done < unique.clusters.list

    # In the scenarip where the whole lineage is not present in the pairwise_clusters file
        if [ ! -s final_clusters.list ]; then

                echo "No clusters found for lineage ${lineage}. Creating empty nexus file."
                echo "" > nexus.tuple.csv

        else

            # Create a CSV for generating a tuple of the paths
                for clusterID in `cat final_clusters.list`; do
                    grep "\${clusterID}" ${clusters_tab} | cut -f1 | tr '\n' ';' > tmp-string
                    echo "${lineage},\${clusterID},${snp_fasta_path},${snp_tab_path},${clusters_tab_path}" >> nexus.tuple.csv
                    rm tmp-string
                done

            # Remove any unclustered clusters (nX-)
                cat nexus.tuple.csv | grep -v ",nX-" > tmp.nexus.tuple.csv
                mv tmp.nexus.tuple.csv nexus.tuple.csv

        fi
    """

}