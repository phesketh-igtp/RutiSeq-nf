#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SINGLE_WF }               from './workflows/single_wf.nf'
include { FILE_CHECK }              from './modules/local/file-checks/main.nf'

workflow {
    
    def color_purple = '\u001B[35m'
    def color_green = '\u001B[32m'
    def color_red = '\u001B[31m'
    def color_reset = '\u001B[0m'

    log.info """
    ${color_purple}
    ╔════════════════════════════════════════════════════════════════════════╗
    ║  ██████╗ ██╗   ██╗████████╗██╗███████╗███████╗ ██████╗                 ║
    ║  ██╔══██╗██║   ██║╚══██╔══╝██║██╔════╝██╔════╝██╔═══██╗                ║
    ║  ██████╔╝██║   ██║   ██║   ██║███████╗█████╗  ██║   ██║                ║
    ║  ██╔══██╗██║   ██║   ██║   ██║╚════██║██╔══╝  ██║▄▄ ██║                ║
    ║  ██║  ██║╚██████╔╝   ██║   ██║███████║███████╗╚██████╔╝                ║
    ║  ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚═╝╚══════╝╚══════╝ ╚══▀▀═╝  v.1.0.0-alpha  ║
    ╚════════════════════════════════════════════════════════════════════════╝
    ${color_reset}
    """

    /*

    Define all the expected argument to be provided at time of running 
        nextflow at CLI

        --samplesheet /path/to/sample-sheet
        --runID [a-zA-Z0-9]
        --workflow [full, single, pairwise, summary, barcoding]

    */

    // Create channel from sample sheet
        if (params.samplesheet == null) {
            error "Please provide a samplesheet CSV file with --samplesheet"
        }

    // Create channel from sample sheet
        Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true, sep: ',')
            .map { row ->
                if (row.sampleID == null || row.forward_path == null || row.reverse_path == null) {
                    error "Missing required column in samplesheet: ${row}"
                }
                tuple(row.sampleID, file(row.forward_path, checkIfExists: true), file(row.reverse_path, checkIfExists: true))
            }
            .set { samples_ch }

        // Report the samples part of the samplesheet
            log.info "${color_purple}Input samples:${color_reset}"
            samples_ch.view { sampleID, forward, reverse ->
                "${color_red}Sample: ${color_green}$sampleID${color_red} | Forward: ${color_green}$forward${color_red} | Reverse: ${color_green}$reverse${color_reset}"
            }

        // Check if the genome has previously been analyzed
        FILE_CHECK(samples_ch)



        // Call the SINGLE_WORKFLOW only for samples missing files
        SINGLE_WF(
            samples_ch, 
            file(params.kaiju_names),
            file(params.kaiju_nodes),
            file(params.kaiju_fmi),
            file(params.tbprofiler_db)
                )

}
