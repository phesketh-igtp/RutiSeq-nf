include { GENERATE_SUMMARY_REPORT } from '../modules/local/summary/summary-report/main.nf'
include { PREPARE_PATHS }           from '../modules/local/summary/prepare-paths/main.nf'
include { PLOT_PHYLOGENY }          from '../modules/local/summary/generate-phylogeny/main.nf'
include { GENERATE_NEXUS }          from '../modules/local/summary/generate-nexus/main.nf'

workflow SUMMARY_WF{

    take:
        runID
        pairwise_clusters
        pairwise_matrix
        analysis_summary
        who_resistance
        tbdb_resistance

    main:

        GENERATE_SUMMARY_REPORT(runID,
                                pairwise_clusters,
                                pairwise_matrix,
                                analysis_summary,
                                who_resistance,
                                tbdb_resistance
                                )


        if (params.metadata) {
            // Channel for metadata file
            ch_metadata = Channel.fromPath(params.metadata)

            TIMETREES( pairwise_clusters,
                                    analysis_summary,
                                    ch_metadata
                    )
            
            // Processes that depend on metadata
            PLOT_PHYLOGENY_WITH_METADATA( pairwise_clusters,
                                    analysis_summary,
                                    ch_metadata
                                    )

            GENERATE_WITH_METADATA( pairwise_clusters,
                            analysis_summary,
                            ch_metadata
                            )

            } else {

                // Processes that don't require metadata
                PROCESS_WITHOUT_METADATA(input_ch)
        }

        PLOT_PHYLOGENY( pairwise_clusters,
                            analysis_summary,
                            metadata
                            )

        GENERATE_NEXUS( pairwise_clusters,
                        analysis_summary
                        )

/*
    emit:
        mjn_positions   =   GENERATE_NEXUS.out.mjn_positions
*/

}