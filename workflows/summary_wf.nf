include { POST_SINGLE_BBDD_CLEANUP   }   from '../modules/local/post-wf-cleaup/pairwise-bbdd-cleanup/main.nf'
include { PROCESS_CLUSTERS           }   from '../modules/local/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT    }   from '../modules/local/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY        }   from '../modules/local/summary/plot-phylogeny/main.nf'
include { PREPARE_NEXUS_PATHS        }   from '../modules/local/summary/prepare-nexus-paths/main.nf'
include { GENERATE_NEXUS             }   from '../modules/local/summary/generate-nexus/main.nf'
include { TABULATE_VARIANT_SITES     }   from '../modules/local/summary/tabulate-variant-positions/main.nf'
include { CONCATENATED_VARIANT_FILES }   from '../modules/local/summary/concatenate-variant-positions/main.nf'
include { POST_SUMMARY_CLEANUP       }   from '../modules/local/summary/post-summary-cleanup-handover/main.nf'
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

    main:

        // Process clusters
            PROCESS_CLUSTERS( runID, pairwise_clusters, analysis_summary )

        // Generate summary XLSX and CSV files for final results    
            GENERATE_SUMMARY_REPORT( runID,
                                        PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        analysis_summary,
                                        who_resistance,
                                        tbdb_resistance
                                    )

        // Plot main ML phylogeny
            PLOT_MAIN_PHYLOGENY( phylogeny_plotting_ch,
                                    PROCESS_CLUSTERS.out.pairwise_clusters_processed )

        

        // Generate base NEXUS files for each cluster
            PREPARE_NEXUS_PATHS( phylogeny_plotting_ch,
                                    PROCESS_CLUSTERS.out.pairwise_clusters_processed )

                nexus_ch = PREPARE_NEXUS_PATHS.out.nexus_tuple
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, clusterID, fasta, tab) = row
                                    tuple(lineage, clusterID, file(fasta), file(tab))
                                    }
                // DEBUG: view the channel 
                ///nexus_ch.view()

            GENERATE_NEXUS( PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                            nexus_ch
                            )

            TABULATE_VARIANT_SITES( GENERATE_NEXUS.out.variant_sites_for_tabulation )

            CONCATENATED_VARIANT_FILES(
                    TABULATE_VARIANT_SITES.out.tabular_vars.collect(),
                    TABULATE_VARIANT_SITES.out.tabular_var_counts.collect()
                    )

        // Cleanup unwanted files
            POST_SUMMARY_CLEANUP( CONCATENATED_VARIANT_FILES.out.cleanup_handover )

}

/*        
        if (params.metadata) {
            // Channel for metadata file
            ch_metadata = Channel.fromPath(params.metadata)

            // Generate base NEXUS files for each cluster
                NEXUS_WITH_METADATA( GENERATE_NEXUS.out.nexus_w_no_metadata,
                                        ch_metadata
                                    )
            // Create timetrees
                GENERATE_TIMETREES( pairwise_clusters,
                                        analysis_summary,
                                        ch_metadata
                                    )

                PLOT_TIMETREES( pairwise_clusters,
                                        analysis_summary,
                                        ch_metadata
                                        )

                } 
            }
*/