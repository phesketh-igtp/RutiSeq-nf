process CONCATENATE_CLUSTERS {

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        path(clusters)

    output:
        path("pairwise_clusters.tsv"),          emit: bbdd_clusters

    script:

        """
        # Create a header
            echo "lineage,distance,genomes,group" > pairwise_clusters.tsv

        # Concatenate all files
            for file in ${clusters}; do cat \$file >> pairwise_clusters.tsv; done
        """  
}