include { GENERATE_SUMMARY_REPORT } from '../modules/local/summary/summary-report/main.nf'
include { PREPARE_PATHS }           from '../modules/local/summary/prepare-paths/main.nf'
include { GENERATE_PHYLOGENY }      from '../modules/local/summary/generate-phylogeny/main.nf'
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
        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        log.info """
        ${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ${color_red}Workflow: ${color_green}Analysis summary${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """

        GENERATE_SUMMARY_REPORT(runID,
                                pairwise_clusters,
                                pairwise_matrix,
                                analysis_summary,
                                who_resistance,
                                tbdb_resistance
                                )

        GENERATE_PHYLOGENY( pairwise_clusters,
                            analysis_summary
                            )

        GENERATE_NEXUS( pairwise_clusters,
                        analysis_summary
                        )

/*
    emit:
        mjn_positions   =   GENERATE_NEXUS.out.mjn_positions
*/

}