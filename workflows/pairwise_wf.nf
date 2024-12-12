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
        pairwise_input_ch // This will be a flattened list of all files for all samples
    
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
            mtbseq_stats_files = pairwise_input_ch.map { tuple ->
                    def (sampleID, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return mtbseq_stats}.collect()
            mtbseq_class_files = pairwise_input_ch.map { tuple ->
                    def (sampleID, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return mtbseq_class}.collect()
            tbdb_out_files = pairwise_input_ch.map { tuple ->
                    def (sampleID, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return tbdb_out}.collect()
            who_out_files = pairwise_input_ch.map { tuple ->
                    def (sampleID, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = tuple
                    return mtbseq_class
                    }.collect()
        
        // Compile TB-Profiler results
            TBPROFILER_COMPILE_TBDB(runID, tbdb_out_files)
                TBPROFILER_COMPILE_TBDB.out.sample_lineage
                    .splitCsv(sep: '\t')
                    .map { sampleID, lineage -> [sampleID, lineage] }
                    .set { ch_sample_lineage }

            TBPROFILER_COMPILE_WHO(runID, who_out_files)

        // Compile stats and classifications from MTBSeq
            MTBSEQ_STATS_COMPILE(runID, 
                                mtbseq_stats_files, 
                                mtbseq_class_files)
            
            MTBSEQ_STATS_COMPILE.out.sample_coverage
                    .splitCsv(sep: '\t')
                    .map { sampleID, coverage -> [sampleID, coverage] }
                    .set { ch_sample_coverage }

}