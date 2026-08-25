include { GENERATE_SUMMARY_REPORT }   from '../modules/summary_wf/summary-report/main.nf'
include { PLOT_MAIN_PHYLOGENY     }   from '../modules/summary_wf/plot-phylogeny/main.nf'
include { NEXUS_GEN               }   from '../modules/summary_wf/prepare-nexus-paths/main.nf'
include { POST_SUMMARY_CLEANUP    }   from '../modules/summary_wf/post-summary-cleanup-handover/main.nf'
include { DATA_DELIVERY           }   from '../modules/summary_wf/data-delivery/main.nf'

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
        NEXUS_GEN( nexus_creation_ch )

        // Collect all outputs before passing to DATA_DELIVERY
            finish_handover = NEXUS_GEN.out.handover_out.collect()

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
                        Remove TimeTree Compoenents due to complications and desire to
                            make the workflow much smaller and simpler for maintanance.
*/