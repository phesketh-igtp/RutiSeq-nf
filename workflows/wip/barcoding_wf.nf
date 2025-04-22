include { SPLIT_CLUSTER_GROUPS }    from '../../modules/local/barcoding/split_into_clusters/main.nf'
include { MERGE_VCFS }              from '../../modules/local/barcoding/merge_vcfs/main.nf'
include { COMPUTE_FTS }             from '../../modules/local/barcoding/compute_Fts/main.nf'

workflow BARCODING_WORKFLOW {

    take:
        runID
        analysis_summary
        pairwise_clusters

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
            SPLIT_CLUSTER_GROUPS( runID,
                                    analysis_summary,
                                    pairwise_clusters
                                )

            // Create a chanel with the following structure from the output of SPLIT_CLUSTERS_GROUPS
            /// [lineage: clusterID, sampleID, vcf_path]
                grouped_vcfs_lineage_ch = SPLIT_CLUSTER_GROUPS.out.vcf_tuple
                    .map { line -> 
                        def (lineage, clusterID, sampleID, vcf_path) = line.tokenize(',')
                        return [lineage, clusterID, sampleID, file(vcf_path)]
                    }
                    .groupTuple(by: 0)
                    .map { lineage, clusterIDs, sampleIDs, vcf_paths ->
                        return [lineage, clusterIDs.flatten(), sampleIDs.flatten(), vcf_paths.flatten()]
                    }

            // Merge the vcfs by the groupings
            MERGE_VCFS( grouped_vcfs_lineage_ch )

            /// group the resulting channel by the clsuterID
                merged_vcfs_cluster_grouped_ch = MERGE_VCFS.out.merged_vcs_tuple
                    .map { clusterID, lineage, sampleIDs, lineagePopList, mergedVcf, mergedVcfIndex ->
                        // Create a tuple with clusterID as the first element for grouping
                        [clusterID, sampleIDs, lineagePopList, mergedVcf, mergedVcfIndex]
                    }
                    .groupTuple(by: 0)
                    .map { clusterID, lineage, sampleIDsList, lineagePopLists, mergedVcfs, mergedVcfIndices ->
                        // Flatten nested lists if necessary and ensure single files for pop lists and VCFs
                        [
                            clusterID,
                            lineage,
                            sampleIDsList.flatten(),
                            lineagePopLists[0],  // Assuming all entries for a clusterID have the same lineage pop list
                            mergedVcfs[0],       // Assuming one merged VCF per clusterID
                            mergedVcfIndices[0]  // Assuming one merged VCF index per clusterID
                        ]
                    }

        /*
            Perform the analysis to calcualte the Fixation index
        */
            COMPUTE_FTS( merged_vcfs_cluster_grouped_ch )

        /* 
            Collect all the individual inputs into a single tsv file, and create a new channel with the input file
        */





}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-04
    @version: 1.0.1
    @description: 
        This is the barcoding workflow for the RutiSeq-nf pipeline. 
        WIP
    @changelog
        v1.0.0-2024-11-01: Initial version
        v1.0.1-2025-04-04: Added documentation and comments
*/