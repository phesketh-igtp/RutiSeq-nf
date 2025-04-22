process SPLIT_CLUSTER_GROUPS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-22
    @version: 1.0.0
    @description: 
        This process creates a CSV file with the following structure that is 
        converted into a datachannel within the main.nf of BARCODING_WF()
        lineage,clusterID,sampleID,vcf_path
    @changelog:
        v1.0.0-2025-04-22:  Initial version
*/

    tag "${runID}"

    input:
        val(runID)
        path(analysis_summary)
        path(pairwise_clusters)

    output:
        path("barcoding_tuple.csv"), emit: vcf_tuple

    script:
        """
        # create the tuple with an R-script
            Rscript -e "
            summary     <- read_delims(${analysis_summary})
            clusters    <- read_delims(${pairwise_clusters})

            tuple <- left_join(summary, clusters)

            tuple <-  tuple |>
                        mutate(vcf_path = paste0(""))
            "
        """

}