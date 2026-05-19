//include { PROCESS_CLUSTERS           }   from '../modules/summary_wf/summary/process-clusters/main.nf'
include { GENERATE_SUMMARY_REPORT    }   from '../modules/summary_wf/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY        }   from '../modules/summary_wf/plot-phylogeny/main.nf'
include { NEXUS_GEN        }   from '../modules/summary_wf/prepare-nexus-paths/main.nf'
//include { GENERATE_NEXUS             }   from '../modules/summary_wf/generate-nexus/main.nf'
include { POST_SUMMARY_CLEANUP       }   from '../modules/summary_wf/post-summary-cleanup-handover/main.nf'
include { DATED_PHYLOGENY         }   from '../modules/summary_wf/generate-timetrees/main.nf'
include { PLOT_TIMETREES             }   from '../modules/summary_wf/plot-timetrees/main.nf'
//include { GENERATE_ANNOTATED_NEXUS   }   from '../modules/summary_wf/summary/generate-nexus-with-metadata/main.nf'
include { GENERATE_NEXUS_W_MRCA      }   from '../modules/summary_wf/generate-nexus-with-ancestor/main.nf'
include { DATA_DELIVERY              }   from '../modules/summary_wf/data-delivery/main.nf'

workflow SUMMARY_WF{

    take:
        processed_clusters
        unprocessed_clusters
        analysis_summary
        who_resistance
        tbdb_resistance
        phylogeny_plotting_ch
        nexus_creation_ch
        dated_phylogeny_ch
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
        NEXUS_GEN( nexus_creation_ch )

        // Collect all outputs before passing to DATA_DELIVERY
            //collected_handover_out = GENERATE_NEXUS.out.handover_out.collect()
            collected_handover_out = NEXUS_GEN.out.handover_out.collect()
        // Conditionally mix with annotated nexus - if metadata provided
        if (params.metadata) {
            // Check if metadata file exists
            def metadata_file = file(params.metadata)
            if (!metadata_file.exists() || metadata_file.isEmpty()) {
                log.warn "${red}Metadata file not found/empty: ${params.metadata}. Skipping time tree generation.${no_col}"
                
                // Fallback to base nexus
                finish_handover = collected_handover_out
            } else {
                log.info "${cyan}Metadata provided. Generating time trees and ancestral sequences.${no_col}"
                
/*              GENERATE_ANNOTATED_NEXUS( 
                                        GENERATE_NEXUS.out.annotated_nexus_ch, 
                                        params.metadata
                                        )
*/

                PLOT_TIMETREES(
                            dated_phylogeny_ch,
                            processed_clusters
                            )
                    timetree_ch = PLOT_TIMETREES.out.timetree_tuple
                            .splitCsv(header: false, sep: ',')
                            .map { row ->
                                def (lineage, clusterID, fasta, tab, ancestor, cluster_tab) = row
                                
                                // Debug: Print the row contents
                                println "DEBUG: Processing row: ${row}"
                                println "DEBUG: fasta=${fasta}, tab=${tab}, ancestor=${ancestor}, cluster_tab=${cluster_tab}"
                                
                                // Check for null values before creating file objects
                                if (!fasta || !tab || !ancestor || !cluster_tab) {
                                    error "One or more file paths are null in row: ${row}"
                                }
                                
                                tuple(lineage, clusterID, file(fasta), file(tab), file(ancestor), file(cluster_tab))
                            }
                GENERATE_NEXUS_W_MRCA(timetree_ch)
                // Mix both annotated nexus channels
                base_nexus_ch = GENERATE_NEXUS_W_MRCA.out.nexus_w_mrca_out
                // Collect all outputs before passing to DATA_DELIVERY
                finish_handover = base_nexus_ch.collect()
            }

        } else {
            log.info "${cyan}No metadata provided. TimeTrees and ancestral sequences will not be generated.${no_col}"
            
            // Only use base nexus channel when no metadata
            finish_handover = collected_handover_out
        }

        // Cleanup unwanted files
            //POST_SUMMARY_CLEANUP( CONCATENATED_VARIANT_FILES.out.cleanup_handover )
            DATA_DELIVERY(
                        sylph_results,
                        reads_taxonomy_qc_report_out,
                        finish_handover
                        )
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
    v2.1.0-2026-01-05: Added conditional handling for metadata input
    v2.2.0-2026-05-19: Merged nexus generation modules into a single one to reduce 
                        number of small processes launched.
*/