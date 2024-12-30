include { SPLIT_CLUSTER_GROUPS }        from './modules/barcoding/split_into_clusters/main.nf'
include { COMPUTE_FTS }        from './modules/barcoding/compute_Fts/main.nf'

workflow BARCODING_WORKFLOW {

    take:
        runID
        analysis_summary
        pairwise_clusters
        mjn_positions

    main:
        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        log.info """
        ${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ${color_red}Workflow: ${color_green}Barcoding analysis${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """

        /*
            Partition VCF files into clusters and lineages so that Fixation analysis can occur. 
            
            At end of this workflow I want generated a TSV file with all the cluster specific SNPs and fixation indexes
                i.e: genome,clusterID,lineage,position,allel,ref,Fts, emit: cluster_specific_snps

        */

            SPLIT_CLUSTER_GROUPS(runID,
                                    analysis_summary,
                                    pairwise_clusters,
                                    mjn_positions
                                )

            // Create a chanel with the following structure from the output of SPLIT_CLUSTERS_GROUPS
            /// samplesID, lineage, cluster, vcf_path, vcf.tbi_path

        /*
            Perform the analysis to calcualte the Fixation index
        */
            COMPUTE_FTS(vcf_ch)

        /* 
            Collect all the individual inputs into a single tsv file, and create a new channel with the input file
        */





}