// Include the SINGLE_WF workflow
include { SINGLE_WF } from '../single_wf.nf'

workflow SINGLE_WF_SUBMIT {
    take:
        single_samples_ch

    main:
        def color_purple = '\u001B[35m'
        def color_green = '\u001B[32m'
        def color_red = '\u001B[31m'
        def color_cyan = '\u001B[36m'
        def no_color = '\u001B[0m'

    // Check if the channel is empty
        single_samples_ch
            .ifEmpty { 
                log.info "${color_red}Workflow execution:${color_green} No single samples to process. ${color_red}SINGLE_WF will not be executed.${no_color}"
                Channel.empty()
            }
            .set { samples_to_process }

    // Only run SINGLE_WF if there are samples to process
        single_results = samples_to_process
            .branch {
                has_samples: it
                no_samples: Channel.empty()
            }
            .map { 
                has_samples = it.has_samples
                if (has_samples) {
                    SINGLE_WF(
                        has_samples,
                        file(params.kaiju_names),
                        file(params.kaiju_nodes),
                        file(params.kaiju_fmi)
                    )
                } else {
                    // Create empty channels to match SINGLE_WF output structure
                    // Adjust the number and names of these channels based on SINGLE_WF's actual output
                    [
                        Channel.empty(),
                        Channel.empty(),
                        Channel.empty()
                    ]
                }
            }

    emit:
        single_results
}