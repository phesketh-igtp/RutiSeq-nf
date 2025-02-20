process CONCATENATE_CLUSTERS {

    publishDir "${params.outdir}/bbdd/results/main/", mode: 'copy'

    input:
        path(clusters)

    output:
        path("pairwise_clusters.tsv"),          emit: bbdd_clusters

    script:

        """
        # Create a header
            echo "lineage\tdistance\tgenomes\tgroup" > pairwise_clusters.tsv

        # Concatenate all files
            for file in ${params.outdir}/bbdd/mtbseq/pairwise/*/Groups/*clusters.tsv; do 
                cat \$file >> pairwise_clusters.tsv
            done

        """

}