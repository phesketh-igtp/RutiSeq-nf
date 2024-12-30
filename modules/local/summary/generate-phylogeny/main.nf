process GENERATE_PHYLOGENY {

    publishDir "${params.outdir}/bbdd/results/phylogeny", mode: 'copy'

    input:
        path pairwise_clusters
        path analysis_summary

    output:
        path("phylogeny/*"), emit: phylogeny_dir

    script:
        """
        """

}