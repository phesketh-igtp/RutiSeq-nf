include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { COMBINE_QC_RESULTS }        from '../modules/local/pre-wf-check/combine-qc-results/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_ANNOTATING_SINGLE }     from '../modules/local/snp-barcoding/single.annotating/main.nf'

workflow SINGLE_WF {

    /*
        Define the inputs from main.nf
    */

    take:
        single_samples_ch
        kaiju_names
        kaiju_nodes
        kaiju_fmi
        tbprofiler_tbdb
        tbprofiler_who

    main:

        /*
            Opening message for workflow
        */

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        log.info """
        ${color_purple}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ${color_red}Workflow: ${color_green}Single genome analysis${color_purple}
        2024-12-09
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~${no_color}
        """          

        // Run MTBC_READ_QC on filtered samples
            MTBC_READ_QC(single_samples_ch,
                        kaiju_names,
                        kaiju_nodes,
                        kaiju_fmi
                        )

        // Collect all QC results
            all_qc_results = MTBC_READ_QC.out.qc_results.map { it[1] }.collect()

        // Combine QC results : NEEDS REPAIRING
            COMBINE_QC_RESULTS( // produce TSV of read QC results
                            all_qc_results, 
                            params.runID)

        // Explicitly capture the mtbc_reads output
            mtbc_reads_ch = MTBC_READ_QC.out.mtbc_reads

            mtbc_reads_ch.view()

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE_TBDB(mtbc_reads_ch, tbprofiler_tbdb)

            TBPROFILER_PROFILE_WHO(mtbc_reads_ch, tbprofiler_who)
    
        // Run MTBSEQ_SINGLE
            MTBSEQ_SINGLE(mtbc_reads_ch)

        // Run SNP_PROFILING_SINGLE using the mpileup output
            SNP_PROFILING_SINGLE(MTBSEQ_SINGLE.out.mtbseq_mpileup)

        // Collect all the output paths from the single analysis and create tuple that is emitted for the final output    
            collected_outputs = MTBSEQ_SINGLE.out.mtbseq_class
                            .mix(MTBSEQ_SINGLE.out.mtbseq_stats)
                            .mix(MTBSEQ_SINGLE.out.mtbseq_pos)
                            .mix(MTBSEQ_SINGLE.out.mtbseq_vars)
                            .mix(TBPROFILER_PROFILE_TBDB.out.tbdb_out)
                            .mix(TBPROFILER_PROFILE_WHO.out.who_out)                          
                            .mix(SNP_PROFILING_SINGLE.out.mtbseq_vcf)
                            .groupTuple(by: 0)

    collected_outputs.view()

    emit:
        analyzed_single_samples_ch = collected_outputs

}