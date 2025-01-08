include { TBPROFILER_COMPILE_TBDB }     from '../modules/local/tbprofiler/compile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }      from '../modules/local/tbprofiler/compile.who/main.nf'
include { MTBSEQ_STATS_COMPILE }        from '../modules/local/mtbseq/stats-compile/main.nf'
include { COMPILE_SEQUENCING_STATS }    from '../modules/local/filtering/compile-sequencing-stats/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE }     from '../modules/local/mtbseq/lineage_pairwise/main.nf'
include { MTBSEQ_LINEAGE_GROUPING }     from '../modules/local/mtbseq/lineage-grouping/main.nf'
include { MTBSEQ_PAIRWISE_RESULTS  }    from '../modules/local/mtbseq/compile-pairwise/main.nf'

workflow PAIRWISE_WF {
    
    take:
        runID
        single_updated_samples_ch

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

        // Create channels of just the necessary outputs contained within the tuple
            mtbseq_stats_files = single_updated_samples_ch.map { tuple ->
                    def (sampleID, forward, reverse, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return mtbseq_stats
                }

            mtbseq_class_files = single_updated_samples_ch.map { tuple ->
                    def (sampleID, forward, reverse, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return mtbseq_class
                }

            tbdb_out_files = single_updated_samples_ch.map { tuple ->
                    def (sampleID, forward, reverse,  mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return tbdb_out
                }

            who_out_files = single_updated_samples_ch.map { tuple ->
                    def (sampleID, forward, reverse,  mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return who_out
                }

            // make the channels
            mtbseq_stats_ch     =       mtbseq_stats_files.collect()
            mtbseq_class_ch     =       mtbseq_class_files.collect()
            tbdb_out_ch         =       tbdb_out_files.collect()
            who_out_ch          =       who_out_files.collect()

        // Compile TB-Profiler results
            TBPROFILER_COMPILE_TBDB( runID, 
                                        tbdb_out_ch
                                    )

            TBPROFILER_COMPILE_WHO( runID, 
                                        who_out_ch
                                    )

            // DEBUG:
                TBPROFILER_COMPILE_TBDB.out.tbdb_results.view()
                TBPROFILER_COMPILE_WHO.out.who_results.view()

        // Compile stats and classifications from MTBSeq
            MTBSEQ_STATS_COMPILE( mtbseq_stats_ch, 
                                    mtbseq_class_ch
                                )
            
        // Determine infection type (Mixed vs Clonal using both tbprofiler and mtbseq outputs)
        //// and filter genomes based on quality parameters (min coverage)
            COMPILE_SEQUENCING_STATS(   runID, 
                                        TBPROFILER_COMPILE_TBDB.out.tbdb_results,
                                        TBPROFILER_COMPILE_WHO.out.who_results,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_strains,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_map_stats
                                    )

            // Create tuple and data channel from lineage_samples_paths.csv
                lineage_samples_ch = COMPILE_SEQUENCING_STATS.out.lineage_sample_path
                    .splitCsv(header:false)
                    .map { row -> tuple(row[0], file(row[1])) }

            // Now you can use lineage_samples_tuple in subsequent processes
                lineage_samples_ch.view()

        // Run the pairwise analysis by lineages
            MTBSEQ_LINEAGE_PAIRWISE(    
                                        runID, 
                                        lineage_samples_ch
                                    )

            MTBSEQ_LINEAGE_GROUPING(
                                        runID, 
                                        lineage_samples_ch,
                                        MTBSEQ_LINEAGE_PAIRWISE.out.amended_dir,
                                        MTBSEQ_LINEAGE_PAIRWISE.out.join_dir
                                    )

            // Collect all cluster and matrix outputs
            all_clusters = MTBSEQ_LINEAGE_GROUPING.out.clusters.collect()
            all_matrices = MTBSEQ_LINEAGE_GROUPING.out.matrices.collect()

        // Compile the pairwise analysis results into a single cluster file
            MTBSEQ_PAIRWISE_RESULTS(all_clusters, all_matrices)


    emit:
        pairwise_clusters       =       MTBSEQ_PAIRWISE_RESULTS.out.master_clusters
        pairwise_matrix         =       MTBSEQ_PAIRWISE_RESULTS.out.master_matrices
        analysis_summary        =       COMPILE_SEQUENCING_STATS.out.analysis_summary
        who_resistance          =       COMPILE_SEQUENCING_STATS.out.who_resistance
        tbdb_resistance         =       COMPILE_SEQUENCING_STATS.out.tbdb_resistance

}