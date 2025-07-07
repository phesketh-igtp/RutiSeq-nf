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

    output:
        path("processed_clusters.tsv"),         emit: pairwise_clusters_processed
        path("${runID}_processed_clusters.tsv")

    script:
    
        """
        # Load the R script for processing clusters
            Rscript ${params.r_script_dir}/process_clusters.R

        # Duplicate the processed clusters file with the runID
            cp processed_clusters.tsv ${runID}_processed_clusters.tsv

        # Copy the processed clusters file to the results directory
            cp processed_clusters.tsv ${params.outDir}/results/
            cp ${pairwise_clusters} ${params.outDir}/results/
        """

}