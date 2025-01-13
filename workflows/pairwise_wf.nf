include { TBPROFILER_COMPILE_TBDB }                 from '../modules/local/tbprofiler/compile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }                  from '../modules/local/tbprofiler/compile.who/main.nf'
include { MTBSEQ_STATS_COMPILE }                    from '../modules/local/mtbseq/stats-compile/main.nf'
include { COMPILE_SEQUENCING_STATS }                from '../modules/local/filtering/compile-sequencing-stats/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE }                 from '../modules/local/mtbseq/lineage_pairwise/main.nf'
include { CONCATENATED_VARIABLE_REGION_PHYLOGENY }  from '../modules/local/phylogeny/concatenated_snp_phylogeny-nf'
include { CONCATENATE_CLUSTERS }                    from '../modules/local/pairwise/concatenate-cluster-file/main.nf'

workflow PAIRWISE_WF {
    
    take:
        runID
        mtbseq_stats_ch
        mtbseq_class_ch
        tbdb_out_ch
        who_out_ch

    main:
    
        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        // Compile TB-Profiler results
            TBPROFILER_COMPILE_TBDB( runID, tbdb_out_ch )
            // DEBUG: 
            ///TBPROFILER_COMPILE_TBDB.out.tbdb_results.view { file -> 
            ///"Content of ${file.name}:\n${file.text}" }

            TBPROFILER_COMPILE_WHO( runID, who_out_ch )
            // DEBUG: 
            ///TBPROFILER_COMPILE_WHO.out.who_results.view { file -> 
            ///"Content of ${file.name}:\n${file.text}" }

            /*
            // DEBUG:
                TBPROFILER_COMPILE_TBDB.out.tbdb_results.view()
                TBPROFILER_COMPILE_WHO.out.who_results.view()
            */

        // Compile stats and classifications from MTBSeq
            MTBSEQ_STATS_COMPILE( mtbseq_stats_ch, mtbseq_class_ch )
            
        // Determine infection type (Mixed vs Clonal using both tbprofiler and mtbseq outputs)
        //// and filter genomes based on quality parameters (min coverage)
            COMPILE_SEQUENCING_STATS(   runID,
                                        TBPROFILER_COMPILE_TBDB.out.tbdb_results,
                                        TBPROFILER_COMPILE_TBDB.out.lineage_fractions,
                                        TBPROFILER_COMPILE_WHO.out.who_results,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_strains,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_map_stats
                                    )

            // Create tuple and data channel from lineage_samples_paths.csv
            /// the channel needs to be grouped by the lineage
                lineage_samples_ch = COMPILE_SEQUENCING_STATS.out.lineage_sample_tuple
                    .splitCsv(header: false)
                    .map { row -> tuple(row[0], row[1]) }
                    .groupTuple()

                // DEBUG: View the grouped channel
                lineage_samples_ch.view()

        // Run the pairwise analysis by lineages
            MTBSEQ_LINEAGE_PAIRWISE( runID, lineage_samples_ch )

            CONCATENATED_VARIABLE_REGION_PHYLOGENY( runID, 
                                    MTBSEQ_LINEAGE_PAIRWISE.out.snp_phylogeny_ch )

            // Collect all cluster and matrix outputs
            bbdd_clusters = MTBSEQ_LINEAGE_PAIRWISE.out.clusters.collect()

            CONCATENATE_CLUSTERS(bbdd_clusters)

    emit:
        pairwise_clusters       =   CONCATENATE_CLUSTERS.out.bbdd_clusters
        analysis_summary        =   COMPILE_SEQUENCING_STATS.out.analysis_summary
        who_resistance          =   COMPILE_SEQUENCING_STATS.out.who_resistance
        tbdb_resistance         =   COMPILE_SEQUENCING_STATS.out.tbdb_resistance
        phylogeny_plotting_ch   =   CONCATENATED_VARIABLE_REGION_PHYLOGENY.out.phylogeny_plotting_ch


}