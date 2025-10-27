//include { PROCESS_CLUSTERS           }   from '../modules/summary_wf/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT    }   from '../modules/summary_wf/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY        }   from '../modules/summary_wf/summary/plot-phylogeny/main.nf'
include { PREPARE_NEXUS_PATHS        }   from '../modules/summary_wf/summary/prepare-nexus-paths/main.nf'
include { GENERATE_NEXUS             }   from '../modules/summary_wf/summary/generate-nexus/main.nf'
include { TABULATE_VARIANT_SITES     }   from '../modules/summary_wf/summary/tabulate-variant-positions/main.nf'
include { CONCATENATED_VARIANT_FILES }   from '../modules/summary_wf/summary/concatenate-variant-positions/main.nf'
include { POST_SUMMARY_CLEANUP       }   from '../modules/summary_wf/summary/post-summary-cleanup-handover/main.nf'
include { GENERATE_TIMETREES         }   from '../modules/summary_wf/summary/generate-timetrees/main.nf'
include { PLOT_TIMETREES             }   from '../modules/summary_wf/summary/plot-timetrees/main.nf'
include { GENERATE_NEXUS_W_MRCA      }   from '../modules/summary_wf/summary/generate-nexus-with-ancestor/main.nf'
include { GENERATE_ANNOTATED_NEXUS   }   from '../modules/summary_wf/summary/generate-nexus-with-metadata/main.nf'
include { DATA_DELIVERY              }   from '../modules/summary_wf/summary/data-delivery/main.nf'
//include { GENERATE_REPORT            }   from '../modules/summary_wf/summary/generate-report/main.nf'

workflow SUMMARY_WF{

    take:
        runID
        processed_clusters
        unprocessed_clusters
        analysis_summary
        who_resistance
        tbdb_resistance
        phylogeny_plotting_ch
        nexus_creation_ch

    main:

    def red     = '\u001B[31m'
    def cyan    = '\u001B[36m'
    def no_col  = '\u001B[0m'

        // Process clusters

           // PROCESS_CLUSTERS( runID, pairwise_clusters, analysis_summary, cluster_handover )

        // Generate summary XLSX and CSV files for final results    
            GENERATE_SUMMARY_REPORT( runID,
                                        processed_clusters,
                                        //PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                        analysis_summary,
                                        who_resistance,
                                        tbdb_resistance
                                    )

        // Plot main ML phylogeny
            PLOT_MAIN_PHYLOGENY( runID, 
                                    phylogeny_plotting_ch,
                                    //PROCESS_CLUSTERS.out.pairwise_clusters_processed,
                                    processed_clusters,
                                    unprocessed_clusters 
                                )

        // Generate base NEXUS files for each cluster
            PREPARE_NEXUS_PATHS( nexus_creation_ch
                                    //phylogeny_plotting_ch,
                                    //PROCESS_CLUSTERS.out.pairwise_clusters_processed
                                    //processed_clusters 
                                )

                nexus_ch = PREPARE_NEXUS_PATHS.out.nexus_tuple
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, clusterID, fasta, tab, clusters_tab) = row
                                    tuple(lineage, clusterID, file(fasta), file(tab), file(clusters_tab))
                                    }

                // DEBUG: view the channel //nexus_ch.view()

            GENERATE_NEXUS( nexus_ch )

            TABULATE_VARIANT_SITES( runID, GENERATE_NEXUS.out.variant_sites_for_tabulation )

            CONCATENATED_VARIANT_FILES(runID, 
                    TABULATE_VARIANT_SITES.out.tabular_vars.collect(),
                    TABULATE_VARIANT_SITES.out.tabular_var_counts.collect()
                    )

        // If metadata is provided then the following modules are run
            if (params.metadata) {
                log.info "${cyan}Metadata provided. Generating time trees and ancestral sequences.${no_col}"
                
                // Channel for metadata file
                ch_metadata = Channel.fromPath(params.metadata)
                    .ifEmpty { error "${red}Metadata file not found/empty: ${params.metadata}. Correct your metadata path/file and resume the analysis with '-resume'${no_col}" }

                // Create timetrees
                GENERATE_TIMETREES( phylogeny_plotting_ch, ch_metadata )

                PLOT_TIMETREES( runID, GENERATE_TIMETREES.out.timetrees_ch, processed_clusters)
                
                timetree_ch = PLOT_TIMETREES.out.timetree_tuple
                        .splitCsv(header: false, sep: ',')
                        .map { row ->
                            def (lineage, clusterID, fasta, tab, ancestor) = row
                            tuple(lineage, clusterID, file(fasta), file(tab), file(ancestor))
                        }

                GENERATE_ANNOTATED_NEXUS( GENERATE_NEXUS.out.annotated_nexus_ch, 
                                            params.metadata,)

                /*
                GENERATE_NEXUS_W_MRCA(timetree_ch, 
                                        processed_clusters
                                    )
                */

            } else {
                log.info "${cyan}No metadata provided. TimeTrees and ancestral sequences will not be generated.${no_col}"
            }

        // Cleanup unwanted files
            //POST_SUMMARY_CLEANUP( CONCATENATED_VARIANT_FILES.out.cleanup_handover )

            DATA_DELIVERY( runID, CONCATENATED_VARIANT_FILES.out.cleanup_handover )

            //GENERATE_REPORT()

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