include { POST_SINGLE_BBDD_CLEANUP }    from '../modules/local/post-wf-cleaup/pairwise-bbdd-cleanup/main.nf'
include { PROCESS_CLUSTERS }            from '../modules/local/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT }     from '../modules/local/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY }         from '../modules/local/summary/plot-phylogeny/main.nf'
include { GENERATE_NEXUS }              from '../modules/local/summary/generate-nexus/main.nf'
//include { GENERATE_TIMETREES }          from '../modules/local/summary/generate-timetrees/main.nf'
//include { PLOT_TIMETREES }              from '../modules/local/summary/plot-timetrees/main.nf'


workflow SUMMARY_WF{

    take:
        runID
        pairwise_clusters
        analysis_summary
        who_resistance
        tbdb_resistance
        phylogeny_plotting_ch
        nexus_creation_ch

    main:

        // Process clusters
            PROCESS_CLUSTERS( runID, pairwise_clusters, analysis_summary )

        // Plot main ML phylogeny
            PLOT_MAIN_PHYLOGENY( PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        phylogeny_plotting_ch )

        // Generate base NEXUS files for each cluster
        /*    GENERATE_NEXUS( nexus_creation_ch, 
                            PROCESS_CLUSTERS.out.pairwise_clusters_processed
                            )
*/
        // Generate summary XLSX and CSV files for final results    
            GENERATE_SUMMARY_REPORT( runID,
                                        PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        analysis_summary,
                                        who_resistance,
                                        tbdb_resistance
                                    )

/*        
        if (params.metadata) {
            // Channel for metadata file
            ch_metadata = Channel.fromPath(params.metadata)

            GENERATE_TIMETREES( pairwise_clusters,
                                    analysis_summary,
                                    ch_metadata
                    )
            
            // Processes that depend on metadata
            PLOT_TIMETREES( pairwise_clusters,
                                    analysis_summary,
                                    ch_metadata
                                    )

            NEXUS_WITH_METADATA( pairwise_clusters,
                            analysis_summary,
                            ch_metadata
                            )

                } 
            }
*/

}