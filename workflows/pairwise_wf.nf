include { TBPROFILER_COMPILE_TBDB }         from '../modules/local/tbprofiler/compile.tbdb/main.nf'
include { TBPROFILER_COMPILE_WHO }          from '../modules/local/tbprofiler/compile.who/main.nf'
include { MTBSEQ_SAMPLE_FILTER }            from '../modules/local/mtbseq/sample_filter/main.nf'
include { MTBSEQ_LINEAGE_SPLITTING }        from '../modules/local/mtbseq/lineage_split/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE }         from '../modules/local/mtbseq/lineage_pairwise/main.nf'
include { MTBSEQ_LINEAGE_PAIRWISE_GROUP }   from '../modules/local/mtbseq/lineage_pairwise_group/main.nf'
include { MTBSEQ_STATS_COMPILE }            from '../modules/local/mtbseq/stats-compile/main.nf'

workflow PAIRWISE_WF {
    
    take:
        runID
        pairwise_input // This will be a flattened list of all files for all samples
    
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

        // Group files by type
            def strain_classifications = pairwise_input.findAll { it.toString().contains("Strain_Classification.tab") }
            def mapping_statistics = pairwise_input.findAll { it.toString().contains("Mapping_and_Variant_Statistics.tab") }
            def position_tables = pairwise_input.findAll { it.toString().contains("Position_Table/") }
            def variant_tables = pairwise_input.findAll { it.toString().contains("Called/") }
            def tbprofiler_results = pairwise_input.findAll { it.toString().contains("tbprofiler/results") }
            def tbprofiler_who_results = pairwise_input.findAll { it.toString().contains("tbprofiler/who-only/results") }


        // Compile TB-Profiler results
            TBPROFILER_COMPILE_TBDB(runID, tbprofiler_results)
                TBPROFILER_COMPILE_TBDB.out.sample_lineage
                    .splitCsv(sep: '\t')
                    .map { sampleID, lineage -> [sampleID, lineage] }
                    .set { ch_sample_lineage }

            TBPROFILER_COMPILE_WHO(runID, tbprofiler_who_results)

        // Compile stats and classifications from MTBSeq
            MTBSEQ_STATS_COMPILE(runID, strain_classifications, mapping_statistics)
                MTBSEQ_STATS_COMPILE.out.sample_coverage
                    .splitCsv(sep: '\t')
                    .map { sampleID, coverage -> [sampleID, coverage] }
                    .set { ch_sample_coverage }

        // Merge the channels using join
            merged_channel = ch_sample_lineage
                .join(ch_sample_coverage, by: 0)  // Join based on the first element (sampleID)
                .join(pairwise_input, by: 0)      // Join the result with pairwise_input

            // View the merged channel (optional)
            merged_channel.view()

        /*
            // Filter the merged channel based on coverage
            filtered_channel = merged_channel.filter { tuple ->
                // Assuming coverage is the third element (index 2) in the tuple
                // Adjust this index if the coverage is at a different position
                def coverage = tuple[2]
                coverage >= params.mtbseq_min_cov}


        // Filter the pairwise_input to include only genomes with the minimum coverage (params.mtbseq.min_cov)
        // and add the lineage to the tuple for splitting afterwards

                // Prepare input for MTBSEQ_SAMPLE_FILTER
                def mtbseq_sample_tables = position_tables.combine(variant_tables, by: 0)
                    .map { sampleID, position_table, variant_table -> 
                        tuple(sampleID, position_table, variant_table)
                    }

            MTBSEQ_SAMPLE_FILTER(runID, 
                                params.mtbseq.min_cov, params.lineage_pairwise
                                mtbseq_sample_tables,
                                MTBSEQ_STATS_COMPILE.out.mtbseq_compiled_map_stats,
                                TBPROFILER_COMPILE_TBDB.out.tbprofile_tdb_compile)
        */

}