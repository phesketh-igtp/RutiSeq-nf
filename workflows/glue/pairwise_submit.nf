// Include the SINGLE_WF workflow
include { PAIRWISE_WF }               from '../pairwise_wf.nf'

workflow PAIRWISE_WF_SUBMIT {
    take:
        runID
        final_pairwise_samples_ch

    main:

        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

    // Check if the channel is empty
    final_pairwise_samples_ch
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
        .set { branched_samples }

    PAIRWISE_WF(params.runID, 
                branched_samples
            )

    emit:
        pairwise_results = PAIRWISE_WF.out
}