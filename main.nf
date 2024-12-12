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
            DEFINE INPUT ARGUMENTS: expected argument to be provided at time of running 
                nextflow at CLI

            nextflow run main.nf \
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

        /*
            CREATE sample_ch FROM SAMPLESHEETS
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

        // Report the samples part of the samplesheet
            log.info "${color_purple}Input samples:${color_reset}"
            samples_ch.view { sampleID, forward, reverse ->
                "${color_red}Sample: ${color_green}$sampleID${color_red} | Forward: ${color_green}$forward${color_red} | Reverse: ${color_green}$reverse${color_reset}"
            }

        /*
            INSPECT BBDD FOR INTERMEDIATE FILES (i.e. sample has been previously analyzed)
        */

        // Check if the genome has previously been analyzed
            FILE_CHECK(samples_ch)

        /*
            SINGLE genome analysis
        */

        // Collect and parse the pairwise samples into the desired structure
            single_samples = FILE_CHECK.out.single_input
                                            .collectFile(name: 'all_single_samples.txt', newLine: true)

            single_samples
                .splitCsv()
                .map { row -> 
                    def (sampleID, forward, reverse) = row
                    tuple(sampleID, 
                        file(forward, checkIfExists: true), 
                        file(reverse, checkIfExists: true))
                }
                .set { single_samples_ch }

        // Call the SINGLE_WORKFLOW only for samples missing files
            SINGLE_WF(
                single_samples_ch, 
                file(params.kaiju_names),
                file(params.kaiju_nodes),
                file(params.kaiju_fmi)
                    )

        /*
            PAIRWISE genome analysis
        */

        // Collect and parse the pairwise samples into the desired structure
            pairwise_samples = FILE_CHECK.out.pairwise_input
                .collectFile(name: 'all_pairwise_samples.txt', newLine: true)

            pairwise_samples.view()

        // Parse the pairwise samples into the desired structure
            pairwise_samples_ch = pairwise_samples
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

        // Create a merged channel that has all the paths for the outputs needed for the pairwise analysis
            pairwise_input_ch = pairwise_samples_ch.mix(
                SINGLE_WF.out.analyzed_single_samples_ch.map { sampleID, files ->
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
                }
            )

        // DEV: Inspect the resulting channel has the expected structure (tuple: sampleID,mtbseq_class,mtbseq_stats,mtbseq_pos,tbdb_out,who_out,mtbseq_vcf)
        pairwise_input_ch.view()

/*
            PAIRWISE_WF(pairwise_input_ch, params.runID)
*/

}
