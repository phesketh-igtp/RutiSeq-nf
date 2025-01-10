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
        comp_samples_ch

    main:

        /*
            Opening message for workflow
        */ 

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        // Check if the input channel is empty

        // Use the branch operator to split the channel
            branched_channel = comp_samples_ch.branch {
                with_reads: it[1] != [] && it[2] != [] // zero-indexed so [1] is the second value in the tuple, ect
                without_reads: it[1] == [] || it[2] == [] }

            sample_ch_skip = branched_channel.without_reads

            /*
            // DEBUG:: View the results
            branched_channel.with_reads.view { "With reads: $it" }
            */

            // Taxonomically classify and partition the MTBC reads
                MTBC_READ_QC(
                            branched_channel.with_reads,
                            )

            // Collect all QC results
                all_qc_results = MTBC_READ_QC.out.qc_results.map { it[1] }.collect()

            // Combine QC results
                COMBINE_QC_RESULTS(all_qc_results, params.runID)

            // Explicitly capture the mtbc_reads output
                mtbc_reads_ch = MTBC_READ_QC.out.mtbc_reads

            // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
                TBPROFILER_PROFILE_TBDB(mtbc_reads_ch,
                                        MTBC_READ_QC.out.updated_sample_ch)

                TBPROFILER_PROFILE_WHO(mtbc_reads_ch,
                                        TBPROFILER_PROFILE_TBDB.out.updated_sample_ch)

            // Run MTBSEQ_SINGLE
                MTBSEQ_SINGLE(mtbc_reads_ch,
                                TBPROFILER_PROFILE_WHO.out.updated_sample_ch)

            // Run SNP_PROFILING_SINGLE using the mpileup output
                SNP_PROFILING_SINGLE(MTBSEQ_SINGLE.out.mtbseq_mpileup,
                                    MTBSEQ_SINGLE.out.updated_sample_ch)

                // create updated channel
                branched_channel_with_reads_updated = SNP_PROFILING_SINGLE.out.updated_sample_ch

            // Merge the processed samples with the samples without reads
            /// TODO: not a 100% certain that the channels are being merged, but will come back and troubleshoot this
                final_updated_sample_ch = branched_channel_with_reads_updated.mix(sample_ch_skip)

            /*
            // View the merged results
                final_updated_sample_ch.view { "Final channel: $it" }
            */

    emit:
        single_updated_samples_ch = final_updated_sample_ch

}