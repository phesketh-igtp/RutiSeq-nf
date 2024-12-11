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

        // Collect and parse the pairwise samples into the desired structure
        single_samples = FILE_CHECK.out.single_input
                                        .collectFile(name: 'all_single_samples.txt', newLine: true)

        single_samples
            .splitCsv()
            .map { row -> 
                def (forward, reverse) = row
                def sampleID = forward.tokenize('/')[-1].split('_')[0]  // Assuming sampleID is the first part of the filename
                tuple(sampleID, 
                    file(forward, checkIfExists: true), 
                    file(reverse, checkIfExists: true))
            }
            .set { single_samples_ch }

        // Collect and parse the pairwise samples into the desired structure
        pairwise_samples = FILE_CHECK.out.pairwise_input
                                        .collectFile(name: 'all_pairwise_samples.txt', newLine: true)

            // Parse the pairwise samples into the desired structure
            pairwise_samples
                .splitCsv()
                .map { row -> 
                    def (mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out) = row
                    def sampleID = mtbseq_class.tokenize('/')[-3]  // Assuming sampleID is the third-to-last part of the path
                    tuple(sampleID, 
                        file(mtbseq_class, checkIfExists: true),
                        file(mtbseq_stats, checkIfExists: true),
                        file(mtbseq_pos, checkIfExists: true),
                        file(mtbseq_vars, checkIfExists: true),
                        file(tbdb_out, checkIfExists: true),
                        file(who_out, checkIfExists: true))
                }
                .set { pairwise_samples_ch }

        // Extract sampleIDs from pairwise_samples_ch
            pairwise_sampleIDs = pairwise_samples_ch.map { it[0] }.collect()

        // Check if all samples from samples_ch are in pairwise_samples_ch
            all_samples_pairwise = samples_ch
                .map { it[0] }  // Extract sampleID from samples_ch
                .collect()
                .map { samples -> 
                    pairwise_sampleIDs.val.containsAll(samples)
                }

}

/*
    // Use a conditional workflow execution
        all_samples_pairwise.branch {
        true: { PAIRWISE_WF(params.runID, pairwise_samples_ch) }
        false: { log.warn "Not all samples have pairwise data. Skipping PAIRWISE_WF." }
        }

*/
        // Now you can access each file type separately
        ///parsed_pairwise_samples.strain_classification.view { "Strain Classification: $it" }
        ///parsed_pairwise_samples.mapping_statistics.view { "Mapping Statistics: $it" }
        ///parsed_pairwise_samples.position_table.view { "Position Table: $it" }
        ///parsed_pairwise_samples.variant_table.view { "Variant Table: $it" }
        ///parsed_pairwise_samples.tbprofiler_results.view { "TB Profiler Results: $it" }
        ///parsed_pairwise_samples.tbprofiler_who_results.view { "TB Profiler WHO Results: $it" }

        // For single samples, you can split them into forward and reverse reads if needed

/*
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
            single_samples, 
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
        */

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