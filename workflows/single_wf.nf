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

        // Check if the input channel is empty
        def collected_outputs
        def all_qc_results
        def mtbc_reads_ch

        single_samples_ch
            .ifEmpty { 
                log.info "No single samples to process. Creating empty 'collected_outputs' channel."
                collected_outputs = Channel.empty()
                return Channel.empty()
            }
            .set { samples_to_process }

        // Process samples if the channel is not empty
            samples_to_process.branch {
                has_samples: it != null
                no_samples: it == null
            }
            .set { branched_samples }

            branched_samples.has_samples.ifEmpty { Channel.empty() }
                
            // Taxonomically classify and partition the MTBC reads
                MTBC_READ_QC(
                            samples_to_process, 
                            kaiju_names, 
                            kaiju_nodes, 
                            kaiju_fmi
                            )

            // Collect all QC results
                all_qc_results = MTBC_READ_QC.out.qc_results.map { it[1] }.collect()

            // Combine QC results
                COMBINE_QC_RESULTS(all_qc_results, params.runID)

            // Explicitly capture the mtbc_reads output
                mtbc_reads_ch = MTBC_READ_QC.out.mtbc_reads

            // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
                TBPROFILER_PROFILE_TBDB(mtbc_reads_ch)

                TBPROFILER_PROFILE_WHO(mtbc_reads_ch)

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
                                    .groupTuple()
                                    .toList()


    emit:
        analyzed_single_samples_ch = branched_samples.has_samples
            ? collected_outputs
            : Channel.empty()

}