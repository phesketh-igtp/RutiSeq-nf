//include { PROCESS_CLUSTERS           }   from '../modules/summary_wf/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT    }   from '../modules/summary_wf/summary/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY        }   from '../modules/summary_wf/summary/plot-phylogeny/main.nf'
include { PREPARE_NEXUS_PATHS        }   from '../modules/summary_wf/summary/prepare-nexus-paths/main.nf'
include { GENERATE_NEXUS             }   from '../modules/summary_wf/summary/generate-nexus/main.nf'
include { POST_SUMMARY_CLEANUP       }   from '../modules/summary_wf/summary/post-summary-cleanup-handover/main.nf'
include { GENERATE_TIMETREES         }   from '../modules/summary_wf/summary/generate-timetrees/main.nf'
include { PLOT_TIMETREES             }   from '../modules/summary_wf/summary/plot-timetrees/main.nf'
include { GENERATE_NEXUS_W_MRCA      }   from '../modules/summary_wf/summary/generate-nexus-with-ancestor/main.nf'
include { GENERATE_ANNOTATED_NEXUS   }   from '../modules/summary_wf/summary/generate-nexus-with-metadata/main.nf'
include { DATA_DELIVERY              }   from '../modules/summary_wf/summary/data-delivery/main.nf'
//include { GENERATE_REPORT            }   from '../modules/summary_wf/summary/generate-report/main.nf'

workflow SUMMARY_WF{

    take:
        processed_clusters
        unprocessed_clusters
        analysis_summary
        who_resistance
        tbdb_resistance
        phylogeny_plotting_ch
        nexus_creation_ch
        sylph_results
        reads_taxonomy_qc_report_out

    main:

    def red     = '\u001B[31m'
    def cyan    = '\u001B[36m'
    def no_col  = '\u001B[0m'

        // Process clusters

        // Generate summary XLSX and CSV files for final results    
            GENERATE_SUMMARY_REPORT(
                                    processed_clusters,
                                    analysis_summary,
                                    who_resistance,
                                    tbdb_resistance,
                                    sylph_results,
                                    reads_taxonomy_qc_report_out
                                    )

        // Plot main ML phylogeny
            PLOT_MAIN_PHYLOGENY(
                                phylogeny_plotting_ch,
                                processed_clusters,
                                unprocessed_clusters 
                                )

        // Generate base NEXUS files for each cluster
            PREPARE_NEXUS_PATHS( 
                                nexus_creation_ch
                                )

                nexus_ch = PREPARE_NEXUS_PATHS.out.nexus_tuple
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, clusterID, fasta, tab, clusters_tab) = row
                                    tuple(lineage, clusterID, file(fasta), file(tab), file(clusters_tab))
                                    }

                // DEBUG: view the channel //nexus_ch.view()

        // Your current code
            GENERATE_NEXUS( nexus_ch )

        // Always get the base nexus channel
            base_nexus_ch = GENERATE_NEXUS.out.annotated_nexus_ch
                .mix(GENERATE_NEXUS.out.annotated_nexus_ch)

        // Collect all outputs before passing to DATA_DELIVERY
            collected_nexus_ch = base_nexus_ch.collect()

        // Conditionally mix with annotated nexus
        if (params.metadata) {
            log.info "${cyan}Metadata provided. Generating time trees and ancestral sequences.${no_col}"

        //    TODO: Repair timetree generation steps
            
            // Channel for metadata file
            ch_metadata = Channel.fromPath(params.metadata)
                .ifEmpty { error "${red}Metadata file not found/empty: ${params.metadata}. Correct your metadata path/file and resume the analysis with '-resume'${no_col}" }
/*
            // Create timetrees
            GENERATE_TIMETREES( phylogeny_plotting_ch, ch_metadata )

            PLOT_TIMETREES(
                        GENERATE_TIMETREES.out.timetrees_ch, 
                        processed_clusters
                        )
            
            timetree_ch = PLOT_TIMETREES.out.timetree_tuple
                    .splitCsv(header: false, sep: ',')
                    .map { row ->
                        def (lineage, clusterID, fasta, tab, ancestor, cluster_tab) = row
                        tuple(lineage, clusterID, file(fasta), file(tab), file(ancestor), file(cluster_tab))
                    }

            GENERATE_ANNOTATED_NEXUS( 
                                    GENERATE_NEXUS.out.annotated_nexus_ch, 
                                    params.metadata
                                    )

            GENERATE_NEXUS_W_MRCA( timetree_ch,
                                    PREPARE_NEXUS_PATHS.out.pairwise_clusters_processed
                                    )

            // Always get the base nexus channel
            base_nexus_ch = GENERATE_NEXUS_W_MRCA.out.annotated_nexus_ch
                .mix(GENERATEGENERATE_NEXUS_W_MRCA_NEXUS.out.annotated_nexus_ch)

        // Collect all outputs before passing to DATA_DELIVERY
            finish_handover = base_nexus_ch.collect()
*/
        } else {
            log.info "${cyan}No metadata provided. TimeTrees and ancestral sequences will not be generated.${no_col}"
            
            // Only use base nexus channel when no metadata
            finish_handover = collected_nexus_ch

        }

        // Cleanup unwanted files
            //POST_SUMMARY_CLEANUP( CONCATENATED_VARIANT_FILES.out.cleanup_handover )

            DATA_DELIVERY(
                        sylph_results,
                        reads_taxonomy_qc_report_out,
                        finish_handover
                        )

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
        v2.0.0-2025-11-15: Remove creation of variant sites tables from summary workflow
*/