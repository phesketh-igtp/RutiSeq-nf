include { SYLPH_CLASSIFICATION }            from '../modules/single_wf/sylph/read_classification/main.nf'
include { SYLPH_CLASSIFICATION_INSPECTION } from '../modules/single_wf/sylph/read_classification_inspection/main.nf'
include { ADAPTORS_AND_DOWNSAMPLING }       from '../modules/single_wf/fastp/adaptor_and_downsampling/main.nf'
include { TBPROFILER_PROFILE_TBDB }         from '../modules/single_wf/tbprofiler/profile.tbdb/main.nf'
include { TBPROFILER_PROFILE_WHO }          from '../modules/single_wf/tbprofiler/profile.who/main.nf'
include { MTBSEQ_SINGLE }                   from '../modules/single_wf/mtbseq/single/main.nf'
include { SNP_PROFILING_SINGLE }            from '../modules/single_wf/snp-barcoding/single.profiling/main.nf'
include { SNP_ANNOTATING_SINGLE }           from '../modules/single_wf/snp-barcoding/single.annotating/main.nf'
include { POST_SINGLE_BBDD_CLEANUP }        from '../modules/single_wf/post-wf-cleaup/single-bbdd-cleanup/main.nf'

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

        //def purple  = '\u001B[35m'
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

        // Sylph read classification as a quality control step
        /// SYLPH_READ_CLASSIFICATION( params.samplesheet )

        // Remove and remaining Illumina adapters and downsample the reads (if necessary)
            ADAPTORS_AND_DOWNSAMPLING( branched_channel.with_reads )

            // Collect failed samples using collectFile
            failed_samples_report = ADAPTORS_AND_DOWNSAMPLING.out.failed_sample_entry
                .collectFile(
                    name: "${params.runID}_failed_samples.txt",
                    storeDir: "${params.outDir}/bbdd/read-qc/",
                    keepHeader: false,
                    skip: { it.size() == 0 }  // Skip empty files
                ) { file ->
                    if (file.size() > 0) {
                        return file.text
                    } else {
                        return ""
                    }
                }

        // Run SYLPH_CLASSIFICATION to classify the reads
            SYLPH_CLASSIFICATION( runID )

                // Identify reads that have mixed taxonomy
                SYLPH_CLASSIFICATION_INSPECTION( SYLPH_CLASSIFICATION.out.sylph_res, 
                                                    failed_samples_report )

        // Run TBPROFILER_PROFILE_TBDB after MTBC_READ_QC is done
            TBPROFILER_PROFILE_TBDB( ADAPTORS_AND_DOWNSAMPLING.out.updated_sample_ch1 )

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

        // Cleanup to reduce storage usage in the publish directory (all of these should be deleted, 
        // this is just ensuring they are properly gone)
            sampleid_list_ch = branched_channel_with_reads_updated.map { it[0] }
                //sampleid_list_ch.view() // check the channel is as you would expect

                POST_SINGLE_BBDD_CLEANUP(sampleid_list_ch)

    emit:
        single_updated_samples_ch   = final_updated_sample_ch
        sylph_results               = SYLPH_CLASSIFICATION.out.sylph_res 
        failed_samples_report       = failed_samples_report

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
*/