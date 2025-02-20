include { TBPROFILER_COMPILE_TBDB }                 from '../modules/local/tbprofiler/compile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }                  from '../modules/local/tbprofiler/compile.who/main.nf'
include { MTBSEQ_STATS_COMPILE }                    from '../modules/local/mtbseq/stats-compile/main.nf'
include { COMPILE_SEQUENCING_STATS1 }                from '../modules/local/filtering/compile-sequencing-stats/main.nf'
include { COMPILE_SEQUENCING_STATS2 }                from '../modules/local/filtering/compile-sequencing-stats/main2.nf'
include { MTBSEQ_LINEAGE_JOINT_AMEND }              from '../modules/local/mtbseq/lineage_joint-amend/main.nf'
include { MTBSEQ_LINEAGE_GROUP }                    from '../modules/local/mtbseq/lineage_group/main.nf'
include { CONCATENATED_VARIABLE_REGION_PHYLOGENY }  from '../modules/local/phylogeny/concatenated_snp_phylogeny-nf'
include { CONCATENATE_CLUSTERS }                    from '../modules/local/pairwise/concatenate-cluster-file/main.nf'

workflow PAIRWISE_WF {
    
    take:
        runID
        mtbseq_stats_ch
        mtbseq_class_ch
        tbdb_out_ch
        who_out_ch
        sampleID_list

    main:
    
        def purple  = '\u001B[35m'
        def green   = '\u001B[32m'
        def red     = '\u001B[31m'
        def cyan    = '\u001B[36m'
        def no_col  = '\u001B[0m'

        // Compile TB-Profiler results
            TBPROFILER_COMPILE_TBDB( runID, tbdb_out_ch )
            TBPROFILER_COMPILE_WHO( runID, who_out_ch )

        // Compile stats and classifications from MTBSeq
            MTBSEQ_STATS_COMPILE( mtbseq_stats_ch, mtbseq_class_ch )

        // Determine infection type (Mixed vs Clonal using both tbprofiler and mtbseq outputs)
        //// and filter genomes based on quality parameters (min coverage)
            COMPILE_SEQUENCING_STATS1(   runID,
                                        TBPROFILER_COMPILE_TBDB.out.tbdb_results,
                                        TBPROFILER_COMPILE_TBDB.out.lineage_fractions,
                                        TBPROFILER_COMPILE_WHO.out.who_results,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_strains,
                                        MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_map_stats,
                                        sampleID_list
                                    )

            COMPILE_SEQUENCING_STATS2(   runID,
                                        COMPILE_SEQUENCING_STATS1.out.pairwise_analysis_list,
                                        sampleID_list
                                    )

            // Create tuple and data channel from lineage_samples_paths.csv
            /// the channel needs to be grouped by the lineage
                lineage_samples_ch = COMPILE_SEQUENCING_STATS.out.lineage_sample_tuple
                    .splitCsv(header: false)
                    .map { row -> tuple(row[0], row[1]) }
                    .groupTuple()

                lineage_samples_ch.view { lineage, samples -> 
                    "${cyan}Clustering - ${green}Lineage: ${red}${lineage}${green} || Genomes: ${red}${samples.size()}${no_col}" }

                skipped_lineages_ch = COMPILE_SEQUENCING_STATS.out.skipped_lineages
                    .splitCsv(header: false)
                    .map { row -> tuple(row[0], row[1]) }
                    .groupTuple()

                skipped_lineages_ch.view { lineage, samples -> 
                    "${purple}Skipping clustering - Lineage: ${lineage} || Genomes: ${samples.size()}${no_col}" }

                // DEBUG: View the grouped channel
                //lineage_samples_ch.view()

        // Run the pairwise analysis by lineages
            MTBSEQ_LINEAGE_JOINT_AMEND( runID, lineage_samples_ch )

            // row[0] lineage, distance, join_dir, amend_dir, samples_txt
            mtbseq_group_ch = MTBSEQ_LINEAGE_JOINT_AMEND.out.mtbseq_group_tuple_csv
                            .splitCsv(header: false, sep: ',')
                            .map { row ->
                                def (lineage, distance, joint_path, amend_path, samples_path) = row
                                tuple(lineage, distance, joint_path, amend_path, samples_path)
                            }
                // DEBUG: View the channel
                //mtbseq_group_ch.view()

            MTBSEQ_LINEAGE_GROUP( runID, mtbseq_group_ch )

            CONCATENATED_VARIABLE_REGION_PHYLOGENY( runID, MTBSEQ_LINEAGE_JOINT_AMEND.out.snp_phylogeny_ch )

        // Collect all cluster and matrix outputs
            bbdd_clusters = MTBSEQ_LINEAGE_GROUP.out.clusters.collect()
            CONCATENATE_CLUSTERS(bbdd_clusters)

        // Create a channel to emit for the nexus generation
            nexus_creation_ch = MTBSEQ_LINEAGE_JOINT_AMEND.out.mtbseq_group_tuple_csv
                            .splitCsv(header: false, sep: ',')
                            .map { row ->
                                def (lineage, distance, joint_path, amend_path, samples_path) = row
                                tuple(lineage, joint_path, amend_path)
                            }
                // DEBUG: View the channel
                //nexus_creation_ch.view()

    emit:
        pairwise_clusters       =   CONCATENATE_CLUSTERS.out.bbdd_clusters
        analysis_summary        =   COMPILE_SEQUENCING_STATS1.out.analysis_summary
        who_resistance          =   COMPILE_SEQUENCING_STATS1.out.who_resistance
        tbdb_resistance         =   COMPILE_SEQUENCING_STATS1.out.tbdb_resistance
        phylogeny_plotting_ch   =   CONCATENATED_VARIABLE_REGION_PHYLOGENY.out.phylogeny_plotting_ch
        nexus_creation_ch       =   nexus_creation_ch

}