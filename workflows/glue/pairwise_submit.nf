// Include the SINGLE_WF workflow
include { PAIRWISE_WF }               from '../pairwise_wf.nf'

workflow PAIRWISE_WF_SUBMIT {

    take:
        runID
        pairwise_samples_ch
        single_results // this is empty but forces PAIRWISE_WF to wait until SINGLE_WF has been resolved

    main:

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

        // Check if the channel is empty
        pairwise_samples_ch
            .ifEmpty { 
                log.info "${color_red}Workflow execution:${color_green} No pairwise samples to process. ${color_red}PAIRWISE_WF will not be executed.${no_color}"
                Channel.empty()
            }
            .set { samples_to_process }

        // Only run SINGLE_WF if there are samples to process
            samples_to_process
                .branch {
                    has_samples: it
                    no_samples: Channel.empty()
                }
                .set { pairwise_branched_samples }

        // Finally run the PAIRWISE_WF
            PAIRWISE_WF(runID, 
                        pairwise_branched_samples.has_samples
                    )

    emit:
        single_results      = single_results
        pairwise_clusters   = PAIRWISE_WF.out.pairwise_clusters
        pairwise_matrix     = PAIRWISE_WF.out.pairwise_matrix
        analysis_summary    = PAIRWISE_WF.out.analysis_summary
        who_resistance      = PAIRWISE_WF.out.who_resistance
        tbdb_resistance     = PAIRWISE_WF.out.tbdb_resistance
        //snp_vcf       = PAIRWISE_WF.out.snp_vcf

}