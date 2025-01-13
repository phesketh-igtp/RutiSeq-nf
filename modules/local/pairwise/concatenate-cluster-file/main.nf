process CONCATENATE_CLUSTERS {

    publishDir "${params.outdir}/bbdd/mtbseq/pairwise/", mode: 'copy'

    input:
        path(bbdd_clusters)

    output:
        path("pairwise_clusters.tsv"),      emit: bbdd_clusters

    script:

        """

        # Concatenate all files, excluding the header from subsequent files

        for file in ${bbdd_clusters}; do

            cat \$file >> pairwise_clusters.tsv

        done
        """
}