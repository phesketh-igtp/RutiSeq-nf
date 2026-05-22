// ASSEMBLE RESULTS
include { TBPROFILER_COMPILE }         from '../modules/pairwise_wf/tbprofiler/compile/main.nf'
include { MTBSEQ_STATS_COMPILE }       from '../modules/pairwise_wf/mtbseq/stats-compile/main.nf'
include { COMPILE_SEQUENCING_STATS }   from '../modules/pairwise_wf/filtering/compile-sequencing-stats/main.nf'
//SNIPPY
include { SNIPPY_LINEAGE_CORE }           from '../modules/pairwise_wf/pairwise/snippy-lineage-core/main.nf'
include { SNIPPY_LINEAGE_CORE_GUBBINS }   from '../modules/pairwise_wf/pairwise/snippy-lineage-core-gubbins/main.nf'
include { SNIPPY_LINEAGE_CORE_PHYLOGENY } from '../modules/pairwise_wf/pairwise/snippy-lineage-core-phylogeny/main.nf'
include { SNIPPY_CORE }                   from '../modules/pairwise_wf/pairwise/snippy-core/main.nf'
include { SNIPPY_PHYLOGENY }              from '../modules/pairwise_wf/pairwise/snippy-phylogeny/main.nf'
// MTBSeq
include { ASSESS_SAMPLES }             from '../modules/pairwise_wf/filtering/assess_samples/main.nf'
include { MTBSEQ_LINEAGE_JOINT_AMEND } from '../modules/pairwise_wf/mtbseq/lineage_joint-amend/main.nf'
include { MTBSEQ_LINEAGE_GROUP }       from '../modules/pairwise_wf/mtbseq/lineage_group/main.nf'
// Work in progress to re-write the inneficient MTBseq Join module in python
    //include { MTBSEQ_LINEAGE_JOINT }       from '../modules/pairwise_wf/mtbseq/lineage_joint/main.nf'
    //include { MTBSEQ_LINEAGE_AMEND }       from '../modules/pairwise_wf/mtbseq/lineage_amend/main.nf'
// Phylogeny
include { SNP_PHYLOGENY }              from '../modules/pairwise_wf/phylogeny/concatenated_snp_phylogeny/main.nf'
// Process Cluster
include { PREPROCESS_CLUSTER }         from '../modules/pairwise_wf/pairwise/preprocess_clusters/main.nf'
include { CONCATENATE_CLUSTERS }       from '../modules/pairwise_wf/pairwise/concatenate-cluster-file/main.nf'

workflow PAIRWISE_WF {

    take:
        sampleID_list

    main:

        def purple  = '\u001B[35m'
        def green   = '\u001B[32m'
        def red     = '\u001B[31m'
        def cyan    = '\u001B[36m'
        def no_col  = '\u001B[0m'

        // Compile TB-Profiler results
        TBPROFILER_COMPILE( 
                            sampleID_list, 
                            )
        // Compile stats and classifications from MTBSeq
        MTBSEQ_STATS_COMPILE(                            
                            sampleID_list, 
                            )
        // Determine infection type (Mixed vs Clonal using both tbprofiler and mtbseq outputs)
        //// and filter genomes based on quality parameters (min coverage)
        COMPILE_SEQUENCING_STATS(
                                TBPROFILER_COMPILE.out.tbdb_out,
                                MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_strains,
                                MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_map_stats,
                                )
        // Assess samples and prepare PW tuples
        ASSESS_SAMPLES(
                        COMPILE_SEQUENCING_STATS.out.pairwise_analysis_list,
                        sampleID_list
                        )
        // Create tuple and data channel from lineage_samples_paths.csv
        // Create tuple and data channel from lineage_samples_paths.csv
        lineage_samples_ch = ASSESS_SAMPLES.out.lineage_sample_tuple
            .splitCsv(header: false)
            .map { row -> tuple(row[0], row[1]) }
            .groupTuple()
            .map { lineage, sampleIDs -> 
                tuple(lineage, sampleIDs, sampleIDs.size())
            }

        skipped_lineages_ch = ASSESS_SAMPLES.out.skipped_lineages_tuple
            .splitCsv(header: false)
            .map { row -> tuple(row[0], row[1]) }
            .groupTuple()

        // Count samples per lineage and total counts
        // FIX: Change to handle 3-element tuple
        lineage_counts_ch = lineage_samples_ch
            .map { lineage, samples, count -> count }  // ← Use count directly
            .reduce(0) { acc, count -> acc + count }

        skipped_lineage_counts_ch = skipped_lineages_ch
            .map { lineage, samples -> samples.size() }
            .reduce(0) { acc, count -> acc + count }

        // Count number of lineages
        // FIX: Change to handle 3-element tuple
        lineage_groups_count = lineage_samples_ch
            .map { lineage, samples, count -> 1 }  // ← Add count parameter
            .reduce(0) { acc, val -> acc + val }

        skipped_lineage_groups_count = skipped_lineages_ch
            .map { lineage, samples -> 1 }
            .reduce(0) { acc, val -> acc + val }

        // Create summary strings for lineages
        // FIX: Change to handle 3-element tuple
        lineage_summary_ch = lineage_samples_ch
            .map { lineage, samples, count -> 
                "  ${cyan}${lineage} (n = ${count})${no_col}"  // ← Use count parameter
            }
            .collect()
            .map { list -> list.join('\n') }

        skipped_lineage_summary_ch = skipped_lineages_ch
            .map { lineage, samples -> 
                "  ${purple}${lineage} (n = ${samples.size()})${no_col}"
            }
            .collect()
            .map { list -> list.join('\n') }

            // Combine and log lineage summary
            lineage_counts_ch
                .combine(skipped_lineage_counts_ch)
                .combine(lineage_groups_count)
                .combine(skipped_lineage_groups_count)
                .combine(lineage_summary_ch)
                .combine(skipped_lineage_summary_ch)
                .subscribe { total_genomes, skipped_genomes, active_lineages, skipped_lineages, lineage_details, skipped_details ->
                    log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                    log.info "${green}PAIRWISE_WF() WORKFLOW SUMMARY${no_col}"
                    log.info "${purple}Clustering analysis grouped by lineage for pairwise comparisons${no_col}"
                    log.info "${purple}only lineages present in the current analysis is performed.${no_col}"
                    log.info "${green}--pairwise_split ${red}${params.pairwise_split}${no_col}"
                    log.info "${green}runID: ${red}${params.runID}${green}"
                    log.info "${green}  Processing : ${red}${active_lineages}${green} new genomes plus database (n = ${red}${total_genomes}${green})"
                    log.info "${green}  Ignored : ${red}${skipped_lineages}${green} lineages (including ${red}${skipped_genomes}${green} genomes from database)${no_col}"
                    log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                    log.info "${green}Active lineages for clustering:${no_col}"
                    log.info "${lineage_details}"
                        if (skipped_genomes > 0) {
                            log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                            log.info "${green}Ignored lineages:${no_col}"
                            log.info "${skipped_details}"
                        }
                    log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                        }
            // DEBUG: View the grouped channel
            //lineage_samples_ch.view()
        // Run the pairwise analysis by lineages
            MTBSEQ_LINEAGE_JOINT_AMEND( lineage_samples_ch )
            SNIPPY_LINEAGE_CORE( lineage_samples_ch )

                // row[0] lineage, distance, join_dir, amend_dir, samples_txt
                mtbseq_group_ch = MTBSEQ_LINEAGE_JOINT_AMEND.out.mtbseq_group_tuple_csv
                                .splitCsv(header: false, sep: ',')
                                .map { row ->
                                    def (lineage, distance, joint_path, amend_path, samples_path, sampleID_count) = row
                                    tuple(lineage, distance, joint_path, amend_path, samples_path, sampleID_count)
                                }
                    // DEBUG: View the channel
                    //mtbseq_group_ch.view()
            MTBSEQ_LINEAGE_GROUP( mtbseq_group_ch )
        // Collect all cluster and matrix outputs
            PREPROCESS_CLUSTER( MTBSEQ_LINEAGE_GROUP.out.clusters )
            processed_clusters_collected = PREPROCESS_CLUSTER.out.pairwise_clusters_processed.collect()
            CONCATENATE_CLUSTERS(
                                processed_clusters_collected, 
                                COMPILE_SEQUENCING_STATS.out.analysis_summary
                                )

        // Assemble all the variable region phylogenies
            SNP_PHYLOGENY( MTBSEQ_LINEAGE_JOINT_AMEND.out.snp_phylogeny_ch )
            SNIPPY_LINEAGE_CORE_PHYLOGENY( SNP_PHYLOGENY.out.snippy_lin_phylo_ch )

        // SNIPPY_CORE and SNIPPY_PHYLOGENY
        if ( params.snippy_core ) { // Only run if enabled in params (--snippy_core 'true')
            
            log.info "${green}Including SNIPPY_CORE analysis${no_col}"

            SNIPPY_CORE( 
                        sampleID_list,
                        COMPILE_SEQUENCING_STATS.out.analysis_summary
                        )

            SNIPPY_PHYLOGENY( 
                            SNIPPY_CORE.out.snippy_core_phylo_alignment
                            )
        /*
            // Generate dated phylogenies if metadata is provided
            if (params.metadata) {
                SNIPPY_DATED_PHYLOGENY( 
                                        SNIPPY_PHYLOGENY.out.snippy_core_dated_phylogeny,
                                        params.metadata
                                    )
        */

        } else {
            log.info "${purple}Skipping SNIPPY_CORE and SNIPPY_PHYLOGENY analysis (params.snippy_core = ${params.snippy_core})${no_col}"
        }

/*
        // Create a channel to emit for the nexus generation
            nexus_creation_ch = MTBSEQ_LINEAGE_JOINT_AMEND.out.mtbseq_group_tuple_csv
                            .splitCsv(header: false, sep: ',')
                            .map { row ->
                                def (lineage, distance, joint_path, amend_path, samples_path) = row
                                tuple(lineage, joint_path, amend_path)
                            }
*/

    emit:
        processed_clusters    = CONCATENATE_CLUSTERS.out.pairwise_clusters_processed
        unprocessed_clusters  = CONCATENATE_CLUSTERS.out.pairwise_clusters_unprocessed
        analysis_summary      = COMPILE_SEQUENCING_STATS.out.analysis_summary
        who_resistance        = COMPILE_SEQUENCING_STATS.out.who_resistance
        tbdb_resistance       = COMPILE_SEQUENCING_STATS.out.tbdb_resistance
        phylogeny_plotting_ch = SNP_PHYLOGENY.out.main_phylogeny_out
        nexus_creation_ch     = PREPROCESS_CLUSTER.out.nexus_ch
}

/*
@author: Poppy J Hesketh Best
@date: 2026-01-19
@version: 1.3.0
@description: 
    This is the pairwise genome workflow for the RutiSeq-nf pipeline.
@changelog
    v1.0.0-2024-11-01: Initial version
    v1.0.1-2025-04-04: Added documentation and comments
    v1.1.0-2026-01-05: Added SNIPPY_DATED_PHYLOGENY module for dated phylogenies from snippy core alignments
    v1.2.0-2026-01-19: Cleaner handling of dated phylogeny emission based on metadata presence
    v1.3.0-2026-05-19: Timetree generations to simplify the entire workflow (was overcomplicating the wf)
*/