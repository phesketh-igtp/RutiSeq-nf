include { POST_SINGLE_BBDD_CLEANUP  }   from '../modules/local/post-wf-cleaup/pairwise-bbdd-cleanup/main.nf'
include { PROCESS_CLUSTERS          }   from '../modules/local/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT   }   from '../modules/local/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY       }   from '../modules/local/summary/plot-phylogeny/main.nf'
include { PREPARE_NEXUS_PATHS       }   from '../modules/local/summary/prepare-nexus-paths/main.nf'
include { GENERATE_NEXUS            }   from '../modules/local/summary/generate-nexus/main.nf'
include { TABULATE_VARIANT_SITES    }   from '../modules/local/summary/tabulate-variant-positions/main.nf'
include { CONCATENATED_VARIANT_FILES}   from '../modules/local/summary/concatenate-variant-positions/main.nf'
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

        // Plot main ML phylogeny
            PLOT_MAIN_PHYLOGENY( PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        phylogeny_plotting_ch )

        // Generate base NEXUS files for each cluster
            lineage_ch = phylogeny_plotting_ch.first()
            lineage_ch.view()

            PREPARE_NEXUS_PATHS( lineage_ch, 
                                    PROCESS_CLUSTERS.out.pairwise_clusters_processed
                            )

                nexus_ch = PREPARE_NEXUS_PATHS.out.nexus_tuple
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, cluster, fasta, tab) = row
                                    tuple(lineage, cluster, fasta, tab)
                                }
                //nexus_ch.view()

            GENERATE_NEXUS( nexus_ch, 
                            PROCESS_CLUSTERS.out.pairwise_clusters_processed
                            )

            TABULATE_VARIANT_SITES( GENERATE_NEXUS.out.variant_sites_for_tabulation)

            CONCATENATED_VARIANT_FILES(
                    TABULATE_VARIANT_SITES.out.tabular_vars.collect(),
                    TABULATE_VARIANT_SITES.out.tabular_var_counts.collect()
                    )

        // Generate summary XLSX and CSV files for final results    
            GENERATE_SUMMARY_REPORT( runID,
                                        PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        analysis_summary,
                                        who_resistance,
                                        tbdb_resistance
                                    )

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