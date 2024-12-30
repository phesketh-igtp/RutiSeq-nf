#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FILE_CHECK }              from './modules/local/file-checks/main.nf'
include { SINGLE_WF_SUBMIT }        from './workflows/glue/single_submit.nf'
include { SINGLE_WF }               from './workflows/single_wf.nf'
include { PAIRWISE_WF }             from './workflows/pairwise_wf.nf'
//include { SUMMARY_WF }              from './workflows/summary_wf.nf'
//include { BARCODING_WF }            from './workflows/barcoding_wf.nf'

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
    ║ ${color_green}Pre-release testing version${color_purple}                                            ║    
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
                    tuple(row.sampleID, file(row.forward_path, checkIfExists: true), 
                            file(row.reverse_path, checkIfExists: true))
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

            // After the FILE_CHECK process
            pairwise_samples = FILE_CHECK.out.pairwise_input
                .collectFile(name: 'all_pairwise_samples.txt', newLine: true, storeDir: params.outdir)
                .ifEmpty { file("${params.outdir}/empty_pairwise_samples.txt") }

            single_samples = FILE_CHECK.out.single_input
                .collectFile(name: 'all_single_samples.txt', newLine: true, storeDir: params.outdir)
                .ifEmpty { file("${params.outdir}/empty_single_samples.txt") }

            // Debug: Print the content of all_pairwise_samples.txt and all_single_samples.txt
            pairwise_samples.view   { "DEBUG - Content of all_pairwise_samples.txt:\n${it.text}" }
            single_samples.view     { "DEBUG - Content of all_single_samples.txt:\n${it.text}" }

            // Parse the pairwise samples into the desired structure
            pairwise_samples_ch = pairwise_samples
                .splitCsv()
                .map { row -> 
                    log.debug "DEBUG - Processing pairwise row: $row"
                    if (row.size() == 8) {
                        def (sampleID, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = row
                        return tuple(
                            sampleID, 
                            file(mtbseq_class.trim()),
                            file(mtbseq_stats.trim()),
                            file(mtbseq_pos.trim()),
                            file(mtbseq_vars.trim()),
                            file(tbdb_out.trim()),
                            file(who_out.trim()),
                            file(mtbseq_vcf.trim())
                        )   
                    } else {
                        log.warn "Skipping pairwise row with incorrect number of elements: $row"
                        return null
                    }
                }
                .filter { it != null }

        // Parse the single samples into the desired structure
        single_samples_ch = single_samples
            .splitCsv()
            .map { row -> 
                log.debug "DEBUG - Processing single row: $row"
                if (row.size() == 3) {
                    def (sampleID, forward, reverse) = row
                    return tuple(sampleID, file(forward.trim()), file(reverse.trim()))
                } else {
                    log.warn "Skipping single row with incorrect number of elements: $row"
                    return null
                }
            }
            .filter { it != null }

            // Debug: Print the content of pairwise_samples_ch and single_samples_ch
            pairwise_samples_ch.view { sample -> "DEBUG - Pairwise sample: $sample" }
            single_samples_ch.view { sample -> "DEBUG - Single sample: $sample" }

        /*
            SINGLE sample analysis
        */
        // Call the SINGLE_WORKFLOW only for samples missing necessary files
            SINGLE_WF(
                        single_samples_ch,
                        params.kaiju_names,
                        params.kaiju_nodes,
                        params.kaiju_fmi
                    )

        /*
            PAIRWISE sample analysis that have all the intermediate documents OR 
            post SINGLE_WF analysis (depends on if intermediate files were present in the BBDD)
        */

        // Call the PAIRWISE_WF con
            PAIRWISE_WF(
                            params.runID, 
                            pairwise_samples_ch,
                            SINGLE_WF_SUBMIT.out.single_results
                        )

        /*
            SUMMARY_WF to generate the XCEL summary tables, produce ML phylogenetic trees 
            and visualise them, and generate MJN files for visualisation in PopArt, and
            interactive python networks
        */
/*
            SUMMARY_WF(     params.runID,
                            PAIRWISE_WF.out.pairwise_clusters,
                            PAIRWISE_WF.out.pairwise_matrix,
                            PAIRWISE_WF.out.analysis_summary,
                            PAIRWISE_WF.out.who_resistance,
                            PAIRWISE_WF.out.tbdb_resistance
                    )
*/

        /*
            BARCODING_WF to perform barcoding analysis of the VCF files generated from the single workflow.
            This analysis has a much lower priority
        */
/*
            BARCODING_WF(
                            params.runID,
                            PAIRWISE_WF.out.pairwise_clusters
                            PAIRWISE_WF.out.analysis_summary
                            SUMMARY_WF.out.mjn_positions
                        )
*/

}