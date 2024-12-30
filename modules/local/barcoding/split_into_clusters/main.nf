process SPLIT_CLUSTER_GROUPS {


    tag "${runID}"

    input:
        val runID
        tuple analysis_summary
        tuple pairwise_clusters
        tuple mjn_positions

    output:


    script:
        """

        """

}