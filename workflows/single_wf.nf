include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
include { COMBINE_QC_RESULTS }        from '../modules/local/pre-wf-check/combine-qc-results/main.nf'
include { TBPROFILER_PROFILE_TBDB }   from '../modules/local/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }    from '../modules/local/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }             from '../modules/local/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }      from '../modules/local/snp-barcoding/single.profiling/main.nf'
include { SNP_ANNOTATING_SINGLE }     from '../modules/local/snp-barcoding/single.annotating/main.nf'
include { POST_SINGLE_BBDD_CLEANUP }  from '../modules/local/post-wf-cleaup/single-bbdd-cleanup/main.nf'

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

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE_TBDB( MTBC_READ_QC.out.updated_sample_ch1 )

            TBPROFILER_PROFILE_WHO( TBPROFILER_PROFILE_TBDB.out.updated_sample_ch2 )

        // Run MTBSEQ_SINGLE
            MTBSEQ_SINGLE( TBPROFILER_PROFILE_WHO.out.updated_sample_ch3 )

        // Run SNP_PROFILING_SINGLE using the mpileup output
            SNP_PROFILING_SINGLE( MTBSEQ_SINGLE.out.updated_sample_ch4 )

            // create updated channel
            branched_channel_with_reads_updated = SNP_PROFILING_SINGLE.out.updated_sample_ch5

        // Merge the processed samples with the samples without reads
            final_updated_sample_ch = branched_channel_with_reads_updated.mix(sample_ch_skip)

            // DEBUG:: 
                final_updated_sample_ch.view { "Final channel: $it" }
            

        // Cleanup to reduce storage usage in the publish directory
            sampleid_list_ch = branched_channel_with_reads_updated.map { it[0] }
            POST_SINGLE_BBDD_CLEANUP(sampleid_list_ch)

    emit:
        single_updated_samples_ch = final_updated_sample_ch

}