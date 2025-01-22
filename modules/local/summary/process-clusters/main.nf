process PROCESS_CLUSTERS {

    conda params.r_stats_env
    
    publishDir "${params.outdir}/bbdd/results/main/", mode: 'copy'

    input:
        val(runID)
        path(pairwise_clusters)
        path(analysis_summary)

    output:
        path("processed_clusters.tsv"),         emit: pairwise_clusters_processed
        path("${runID}_processed_clusters.tsv")

    script:
    """
    Rscript ${params.r_script_dir}/process_clusters.R \\
            --clusters ${pairwise_clusters} \\
            --summary ${analysis_summary}

    cp processed_clusters.tsv ${runID}_processed_clusters.tsv
    """

}