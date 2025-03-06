process CONCATENATE_CLUSTERS {

    publishDir "${params.outdir}/bbdd/results/main/", mode: 'copy'

    input:
        path(clusters)

    output:
        path("unprocessed_clusters.tsv"),          emit: bbdd_clusters

    script:

        """
        # Create a header
            echo "lineage\tdistance\tgenomes\tgroup" > unprocessed_clusters.tsv

        # Concatenate all files
            for file in ${params.outdir}/bbdd/mtbseq/pairwise/*/Groups/*clusters.tsv; do 
                cat \$file >> unprocessed_clusters.tsv
            done

        """

}