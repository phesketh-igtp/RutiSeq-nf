include { POST_SINGLE_BBDD_CLEANUP   }   from '../modules/local/post-wf-cleaup/pairwise-bbdd-cleanup/main.nf'
include { PROCESS_CLUSTERS           }   from '../modules/local/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT    }   from '../modules/local/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY        }   from '../modules/local/summary/plot-phylogeny/main.nf'
include { PREPARE_NEXUS_PATHS        }   from '../modules/local/summary/prepare-nexus-paths/main.nf'
include { GENERATE_NEXUS             }   from '../modules/local/summary/generate-nexus/main.nf'
include { TABULATE_VARIANT_SITES     }   from '../modules/local/summary/tabulate-variant-positions/main.nf'
include { CONCATENATED_VARIANT_FILES }   from '../modules/local/summary/concatenate-variant-positions/main.nf'
include { POST_SUMMARY_CLEANUP       }   from '../modules/local/summary/post-summary-cleanup-handover/main.nf'
include { GENERATE_TIMETREES         }   from '../modules/local/summary/generate-timetrees/main.nf'
include { PLOT_TIMETREES             }   from '../modules/local/summary/plot-timetrees/main.nf'
include { GENERATE_NEXUS_W_MRCA      }   from '../modules/local/summary/generate-nexus-with-ancestor/main.nf'
include { PREPARE_DATA_DELIVERY      }   from '../modules/local/summary/data-delivery/main.nf'
//include { GENERATE_NEXUS_W_METADATA  }   from '../modules/local/summary/generate-nexus-with-metadata/main.nf'

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
            PLOT_MAIN_PHYLOGENY( runID, phylogeny_plotting_ch,
                                    PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                    pairwise_clusters )

        // Generate base NEXUS files for each cluster
            PREPARE_NEXUS_PATHS( phylogeny_plotting_ch,
                                    PROCESS_CLUSTERS.out.pairwise_clusters_processed )

                clusters_ch = PROCESS_CLUSTERS.out.pairwise_clusters_processed

                nexus_ch = PREPARE_NEXUS_PATHS.out.nexus_tuple
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, clusterID, fasta, tab) = row
                                    tuple(lineage, clusterID, file(fasta), file(tab))
                                    }
                // DEBUG: view the channel 
                ///nexus_ch.view()

            GENERATE_NEXUS( runID, clusters_ch, nexus_ch )

            TABULATE_VARIANT_SITES( runID, GENERATE_NEXUS.out.variant_sites_for_tabulation )

            CONCATENATED_VARIANT_FILES(runID, 
                    TABULATE_VARIANT_SITES.out.tabular_vars.collect(),
                    TABULATE_VARIANT_SITES.out.tabular_var_counts.collect()
                    )

        // If metadata is provided then the following modules are run
            if (params.metadata) {
                log.info "Metadata provided. Generating time trees and ancestral sequences."
                
                // Channel for metadata file
                ch_metadata = Channel.fromPath(params.metadata)
                    .ifEmpty { error "Metadata file not found/empty: ${params.metadata}. Correct your metadata path/file and resume the analysis with '-resume'" }

                // Create timetrees
                GENERATE_TIMETREES(phylogeny_plotting_ch, ch_metadata)

                PLOT_TIMETREES( runID, GENERATE_TIMETREES.out.timetrees_ch, clusters_ch)
                
                timetree_ch = PLOT_TIMETREES.out.timetree_tuple
                        .splitCsv(header: false, sep: ',')
                        .map { row ->
                            def (lineage, clusterID, fasta, tab, ancestor) = row
                            tuple(lineage, clusterID, file(fasta), file(tab), file(ancestor))
                        }

                GENERATE_NEXUS_W_MRCA(runID, timetree_ch, clusters_ch)

                /*
                GENERATE_NEXUS_W_METADATA(GENERATE_NEXUS_W_ANCESTOR.out.nexus_w_no_metadata,
                                        ch_metadata, clusters_ch
                                    )
                */

            } else {

                log.info "No metadata provided. TimeTrees and ancestral sequence generation."
                
            }

        // Cleanup unwanted files
            //POST_SUMMARY_CLEANUP( CONCATENATED_VARIANT_FILES.out.cleanup_handover )

            PREPARE_DATA_DELIVERY( runID, CONCATENATED_VARIANT_FILES.out.cleanup_handover )



}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-04
    @version: 1.0.1
    @description: 
        This is the summary workflow for the RutiSeq-nf pipeline.
    @changelog
        v1.0.0-2024-11-01: Initial version
        v1.0.1-2025-04-04: Added documentation and comments
*/