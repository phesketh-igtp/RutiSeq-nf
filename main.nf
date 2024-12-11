#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SINGLE_WF }               from './workflows/single_wf.nf'
include { PAIRWISE_WF }             from './workflows/pairwise_wf.nf'
include { FILE_CHECK }              from './modules/local/file-checks/main.nf'
//include { SUMMARY_WORKFLOW }        from './workflows/summary_wf.nf'
//include { BARCODING_WORKFLOW }      from './workflows/barcoding_wf.nf'
//include { REMOVE_SAMPLE_WORKFLOW }  from './workflows/remove-sample_wf.nf'

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

    // Check if RunID is provided
        if (params.runID == null ) {
            error "Please provide a RunID using --runID"
        }

    /* UNDER CONSTRUCTION!
    // Validate workflow option
        def valid_workflows = ['full', 'single', 'pairwise', 'summary', 'barcoding']
        if (!valid_workflows.contains(params.workflow)) {
            error "Invalid workflow option. Please choose from: ${valid_workflows.join(', ')}"
        }   
    */ 

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

/*
    // Report the samples part of the samplesheet
        log.info "${color_purple}Input samples:${color_reset}"
        samples_ch.view { sampleID, forward, reverse ->
            "${color_red}Sample: ${color_green}$sampleID${color_red} | Forward: ${color_green}$forward${color_red} | Reverse: ${color_green}$reverse${color_reset}"
        }
*/

    // Run FILE_CHECK process
        FILE_CHECK(samples_ch)

    // Collect all pairwise samples
        pairwise_samples = FILE_CHECK.out.pairwise_input.collect()

    // Process pairwise samples
        pairwise_input = pairwise_samples
            .filter { it.size() > 0 }
            .map { samples -> 
                def all_files = samples.flatten()
                return all_files
            }

    // Log pairwise input
        pairwise_input.subscribe { files ->
            log.info "Pairwise input files: $files"
        }

    // Call the SINGLE_WORKFLOW only for samples missing files
        SINGLE_WF(
            FILE_CHECK.out.single_input, 
            file(params.kaiju_names),
            file(params.kaiju_nodes),
            file(params.kaiju_fmi),
            file(params.tbprofiler_db)
        )

    // Call PAIRWISE_WF only if there are pairwise samples
        pairwise_input
            .filter { it.size() > 0 }
            .ifEmpty { log.warn "No pairwise samples found. Skipping PAIRWISE_WF." }
            .set { final_pairwise_input }

        PAIRWISE_WF(params.runID, final_pairwise_input)
        
}

    /*
    //
    SUMMARY_WORKFLOW()

    //
    BARCODING_WORKFLOW()

    //
    Channel
        .fromPath(params.removelist)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            if (row.sampleID == null ) {
                error "Missing required column in samplesheet: ${row}"
            }
            tuple(row.sampleID)
        }
        .set { remove_ch }
    REMOVE_SAMPLE_WORKFLOW(remove_ch)
    */