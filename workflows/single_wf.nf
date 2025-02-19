include { MTBC_READ_QC }              from '../modules/local/pre-wf-check/mtbc-reads-qc/main.nf'
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
        runID
        comp_samples_ch

    main:

        /*
            Opening message for workflow
        */ 

        def purple  = '\u001B[35m'
        def green   = '\u001B[32m'
        def red     = '\u001B[31m'
        def cyan    = '\u001B[36m'
        def no_col  = '\u001B[0m'

        // Check if the input channel is empty

        // Use the branch operator to split the channel
            branched_channel = comp_samples_ch.branch {
                with_reads: it[1] != [] && it[2] != [] // zero-indexed so [1] is the second value in the tuple, ect
                without_reads: it[1] == [] || it[2] == [] }

            sample_ch_skip = branched_channel.without_reads

            // Count the number of samples in each channel
                with_reads_count = branched_channel.with_reads.count()
                without_reads_count = sample_ch_skip.count()

            // Combine the counts and log the message
                with_reads_count.combine(without_reads_count)
                    .map { with_reads, without_reads -> 
                        log.info "${green}runID: ${red}${runID}${green} || For ${cyan}SINGLE_WF()${green} : ${red}${with_reads}${green} samples || Skipped until ${cyan}PAIRWISE()${green}: ${red}${without_reads}${green} samples${no_col}"
                    }

        /*
        // DEBUG:: View the results
            branched_channel.with_reads.view { "With reads: $it" }
        */

        // Taxonomically classify and partition the MTBC reads
            MTBC_READ_QC( branched_channel.with_reads )

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
                //final_updated_sample_ch.view { "Final channel: $it" }

        // Cleanup to reduce storage usage in the publish directory (all of these should be deleted, this is just ensuring they are properly gone)
            sampleid_list_ch = branched_channel_with_reads_updated.map { it[0] }

                sampleid_list_ch = branched_channel_with_reads_updated
                    .map { it[0] }
                    .tap { count -> println "${green}Current completed ${cyan}SINGLE_WF()${green} samples: ${red}${count.count()}${no_col}" }

    emit:
        single_updated_samples_ch = final_updated_sample_ch

}