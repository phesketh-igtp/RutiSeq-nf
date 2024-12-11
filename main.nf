#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SINGLE_WF }               from './workflows/single_wf.nf'
include { FILE_CHECK }              from './modules/local/file-checks/main.nf'
include { PAIRWISE_WF }             from './workflows/pairwise_wf.nf'

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
            error "Please provide a samplesheet CSV file with --samplesheet (csv)"
        }

    // Create channel from sample sheet
        if (params.runID == null) {
            error "Please provide a runID file with --runID (chr)"
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

        FILE_CHECK.out.single_input.view()

        // Collect and parse the pairwise samples into the desired structure
        single_samples = FILE_CHECK.out.single_input
                                        .collectFile(name: 'all_single_samples.txt', newLine: true)
        single_samples.view()                                

        single_samples
            .splitCsv()
            .map { row -> 
                def (sampleID, forward, reverse) = row
                tuple(sampleID, 
                    file(forward, checkIfExists: true), 
                    file(reverse, checkIfExists: true))
            }
            .set { single_samples_ch }

        single_samples_ch.view()

        // Call the SINGLE_WORKFLOW only for samples missing files
        SINGLE_WF(
            single_samples_ch, 
            file(params.kaiju_names),
            file(params.kaiju_nodes),
            file(params.kaiju_fmi),
            file(params.tbprofiler_db)
                )

        
        // Collect and parse the pairwise samples into the desired structure
        pairwise_samples = FILE_CHECK.out.pairwise_input
                                        .collectFile(name: 'all_pairwise_samples.txt', newLine: true)

            // Parse the pairwise samples into the desired structure
            pairwise_samples
                .splitCsv()
                .map { row -> 
                    def (sampleID,mtbseq_class,mtbseq_stats,mtbseq_pos,mtbseq_vars,tbdb_out,who_out,mtbseq_vcf) = row
                    tuple(sampleID, 
                        file(mtbseq_class),
                        file(mtbseq_stats),
                        file(mtbseq_pos),
                        file(mtbseq_vars),
                        file(tbdb_out),
                        file(who_out),
                        file(mtbseq_vcf))
                }
                .set { pairwise_samples_ch }
        
        pairwise_samples.view()

        // Create a merged channel that has all the paths for the outputs needed for the pairwise analysis
        // since the wf will wait for the output from SINGLE_WF, this is how im thinking i get around
        // the expanding BBDD issue and how nextflow isnt really intended for this kind of processes.
        // When merging the WFs, this channel will be emitted and fed into the pairwise workflow
        pairwise_input_channel = pairwise_samples_ch
            .mix(SINGLE_WF.out.analyzed_single_samples_ch.map { sampleID, files ->
                tuple(
                    sampleID,
                    files[0], // mtbseq_class
                    files[1], // mtbseq_stats
                    files[2], // mtbseq_pos
                    files[3], // mtbseq_vars
                    files[4], // tbdb_out
                    files[5], // who_out
                    files[6]  // mtbseq_vcf
                )
            })

        pairwise_input_channel.view()


    PAIRWISE_WF(pairwise_input_channel,
                params.runID)

}