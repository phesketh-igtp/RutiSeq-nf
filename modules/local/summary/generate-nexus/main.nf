process GENERATE_NEXUS {

    publishDir "${params.outdir}/bbdd/results/networks", mode: 'copy'

    input:
        path pairwise_clusters
        path analysis_summary

    output:

        path("mj-networks/nexus/*"),    emit: nexus_dir
        path("mj-networks/fasta/*")
        path("mj-networks/positions/*")

    script:
        """
        """

}