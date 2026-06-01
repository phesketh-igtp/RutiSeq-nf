include { ADAPTORS_AND_DOWNSAMPLING }       from '../modules/single_wf/fastp/adaptor_and_downsampling/main.nf'
include { TBPROFILER_PROFILE }              from '../modules/single_wf/tbprofiler/profiler/main.nf'
include { MTBSEQ_SINGLE }                   from '../modules/single_wf/mtbseq/main.nf'
include { SNIPPY_SINGLE }                   from '../modules/single_wf/snippy/main.nf'
include { VARSCAN_SINGLE }                  from '../modules/single_wf/varscan/main.nf'
include { POST_SINGLE_DB_CLEANUP }          from '../modules/single_wf/post-wf-cleaup/single-db-cleanup/main.nf'

workflow SINGLE_WF {

    /*
        Define the inputs from main.nf
    */

    take:
        comp_samples_ch
        tbprofiler_db

    main:

        /*
            Opening message for workflow
        */ 

        //def purple  = '\u001B[35m'
        def green   = '\u001B[32m'
        def red     = '\u001B[31m'
        def cyan    = '\u001B[36m'
        def purple  = '\u001B[35m'
        def no_col  = '\u001B[0m'

        // Check if the input channel is empty

        // Use the branch operator to split the channel
            branched_channel = comp_samples_ch.branch {
                with_reads: it[1] != []     // Has at least one read file
                without_reads: it[1] == []  // No reads at all
            }

            skipped_samples_ch = branched_channel.without_reads

            // Count the number of samples in each channel
                nonskipped_samples_count = branched_channel.with_reads.count()
                skipped_samples_count = skipped_samples_ch.count()

        // Extract [0] index and collect into a list, then join with commas
        skipped_samples_names = skipped_samples_ch.map { it[0] }.collect().map { it.join(', ') }

        // Combine and log
        nonskipped_samples_count
            .combine(skipped_samples_count)
            .combine(skipped_samples_names)
            .subscribe { with_reads, without_reads, skipped_names_string ->
                log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                log.info "${green}SINGLE_WF WORKFLOW SUMMARY${no_col}"
                log.info "${purple}This workflow is additive, meaning if a sampleID already exists in the central database${no_col}"
                log.info "${purple}and will not be re-analised. The sampleID must be manually removed/moved to another location.${no_col}"
                log.info "${green}runID: ${red}${params.runID}${green} || For ${cyan}SINGLE_WF()${green} : ${red}${with_reads}${green} samples || Skipped${green}: ${red}${without_reads}${green} samples${no_col}"
                log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
                log.info "${green}Skipped samples are as follows:${no_col}"
                log.info "${purple}${skipped_names_string}${no_col}"
                log.info "${red}⚠️  If you expected a samples to be analysed, rename you sampleID or remove duplicated sampleIDs in:${no_col}"
                log.info "      - ${purple}${params.outDir}/db/samples/${no_col}"
                log.info "${green}-----------------------------------------------------------------------------------------${no_col}"
            }

        /*
        // DEBUG:: View the results
            branched_channel.with_reads.view { "With reads: $it" }
        */

        // Sylph read classification as a quality control step
        /// SYLPH_READ_CLASSIFICATION( params.samplesheet )

        // Remove and remaining Illumina adapters and downsample the reads (if necessary)
            ADAPTORS_AND_DOWNSAMPLING( branched_channel.with_reads )
            //ADAPTORS_AND_DOWNSAMPLING_SE(branched_channel.with_reads_se )

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE( 
                                ADAPTORS_AND_DOWNSAMPLING.out.updated_sample_ch1,
                                tbprofiler_db
                                )

        // Run MTBSEQ_SINGLE
            MTBSEQ_SINGLE( TBPROFILER_PROFILE.out.updated_sample_ch2 )

        // Run SNP_PROFILING_SINGLE using the mpileup output
            SNIPPY_SINGLE( MTBSEQ_SINGLE.out.updated_sample_ch3 )
            //VARSCAN_SINGLE( SNIPPY_SINGLE.out.snippy_bam )

            // create updated channel
            branched_channel_with_reads_updated = SNIPPY_SINGLE.out.updated_sample_ch4

        // Merge the processed samples with the samples without reads
            final_updated_sample_ch = branched_channel_with_reads_updated.mix(skipped_samples_ch)

            // DEBUG:: 
                //final_updated_sample_ch.view { "Final channel: $it" }

        // TODO: change all the modules to just be input of fastp reads, then mix all the output channels into a specific tuple structure for the final output/emit
        // Cleanup to reduce storage usage in the publish directory (all of these should be deleted, 
        // this is just ensuring they are properly gone)
            sampleid_list_ch = branched_channel_with_reads_updated.map { it[0] }
                //sampleid_list_ch.view() // check the channel is as you would expect

                POST_SINGLE_DB_CLEANUP(sampleid_list_ch)

    emit:
        single_updated_samples_ch   = final_updated_sample_ch

}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-04
@version: 1.0.1
@description: 
    This is the single genome workflow for the RutiSeq-nf pipeline.
@changelog
    v1.0.0-2024-11-01: Initial version
    v1.0.1-2025-04-04: Added documentation and comments
    v2.0.0-2025-11-13: Restructured workflow to merge TB-Profiler into a single step
                    Changed the channels outputs with the removal of one module.
                    Extended warning message for skipped samples, added path of samples db.
*/