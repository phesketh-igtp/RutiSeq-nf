process PROCESS_CLUSTERS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process wrangles the cluster summary data and generates a processed clusters file.
        The resulting file has the clsuter IDs (numeric and not consistent across different runs)
        and asd the lineage ID to the numeric value to ensure that there are no duplicates. This 
        file is used for the summary report and the phylogeny plotting.
*/

    conda params.r_stats_env

    publishDir "${params.outDir}/results/${runID}/clusters/", mode: 'copy'

    input:
        val(runID)
        path(pairwise_clusters)
        path(analysis_summary)
        path(cluster_handover)

    output:
        path("processed_clusters.tsv"),         emit: pairwise_clusters_processed
        path("${runID}_processed_clusters.tsv")
        path("${runID}_singletons.tsv")

    script:
    
        """
        # Load the R script for processing clusters
            Rscript ${params.scriptDir}/R/process_clusters.R

        # Split the cluster into singletons and processed clusters
            grep 'nX-' processed_clusters.tsv > ${runID}_singletons.tsv
            grep -v 'nX-' processed_clusters.tsv > ${runID}_processed_clusters.tsv
            cp ${runID}_processed_clusters.tsv processed_clusters.tsv 
        """

}