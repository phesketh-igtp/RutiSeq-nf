include { TBPROFILER_COMPILE_TBDB }         from '../modules/local/tbprofiler/compile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }          from '../modules/local/tbprofiler/compile.who/main.nf'
include { MTBSEQ_SAMPLE_FILTER }            from '../modules/local/mtbseq/sample_filter/main.nf'
include { MTBSEQ_LINEAGE_SPLITTING }        from '../modules/local/mtbseq/lineage_split/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE }         from '../modules/local/mtbseq/lineage_pairwise/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE_GROUP }   from '../modules/local/mtbseq/lineage_pairwise_group/main.nf'

workflow PAIRWISE_WF {
    
    take:
    runID
    workflow_complete

    main:
        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        log.info """
        ${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ${color_red}Workflow: ${color_green}Comparative/Pairwise genome analysis${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """

        // Compile TB-Profiler results
        TBPROFILER_COMPILE_TBDB(runID,
                                )

        TBPROFILER_COMPILE_WHO( runID,
                                )

        // Compile stats and classifications from MTBSeq
        MTBSEQ_SAMPLE_FILTER(params.mtbseq.min_cov,
                            runID,
                            TBPROFILER_COMPILE_TBDB.out.tbprofile_tdb_compile)
            
        MTBSEQ_SAMPLE_FILTER.out.mtbseq_join_paths
            .splitCsv(header: false)
            .map { row -> tuple(row[0], file(row[1]), file(row[2])) }
            .set { mtbseq_join_paths_ch }

        mtbseq_join_paths_ch.view { it -> "MTBseq join paths: $it" }

        // Uncomment and adjust the following sections as needed
        /*
        // Split the genomes into lineages based on params.
        MTBSEQ_LINEAGE_SPLITTING(runID,
                                TBPROFILER_COMPILE_TBDB.out.tbprofiler_tbdb_compile,
                                Channel.fromList(params.lineage_pairwise))
    
        // Create a channel for all the mtbseq paths needed for compiling
        MTBSEQ_LINEAGE_SPLITTING.out.mtbseq_paths
                .splitCsv(header: false)
                .map { row -> tuple(row[0], row[1], file(row[2]), file(row[3])) }
                .set { mtbseq_paths_ch }

        // Split the genomes into lineages based on params.
        MTBSEQ_LINEAGE_PAIRWISE(runID, 
                                mtbseq_paths_ch,
                                mtbseq_called_results,
                                mtbseq_pos_var_results)

        MTBSEQ_LINEAGE_PAIRWISE_GROUP(runID, 
                                    MTBSEQ_LINEAGE_PAIRWISE.out,
                                    Channel.fromList(params.mtbseq.snp_distance))
        */
    /*
    emit:
        tbprofiler_tbdb_results = TBPROFILER_COMPILE_TBDB.out.tbprof_tbdb_res
        tbprofiler_who_results = TBPROFILER_COMPILE_WHO.out.tbprof_who_res
        mtbseq_filtered_results = MTBSEQ_SAMPLE_FILTER.out.mtbseq_join_paths
    */
        // Uncomment the following lines if you uncomment the corresponding processes above
        /*
        lineage_splitting_results = MTBSEQ_LINEAGE_SPLITTING.out.lineages_ch
        lineage_pairwise_results = MTBSEQ_LINEAGE_PAIRWISE.out.lineage_results
        lineage_pairwise_group_results = MTBSEQ_LINEAGE_PAIRWISE_GROUP.out.group_results
        */
}