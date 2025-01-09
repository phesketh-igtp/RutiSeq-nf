#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FILE_CHECK }              from './modules/local/file-checks/main.nf'
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
            Channel.fromPath(params.samplesheet)
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                        if (row.sampleID == null || row.forward_path == null || row.reverse_path == null || row.type == null) {
                            error "Missing required column in samplesheet: ${row}"
                        }
                        tuple( row.sampleID, 
                                file(row.forward_path, checkIfExists: true), 
                                file(row.reverse_path, checkIfExists: true), 
                                row.type
                                )
                    }
                    .branch {
                        sample: it[3] == 'sample'
                        control: it[3] == 'control'
                    }
                    .set { branched_samples }

        // Create the sample_ch and control_ch
            samples_ch = branched_samples.sample
            controls_ch = branched_samples.control

        // Report the samples part of the samplesheet
            log.info "${color_green}Input samples:${color_reset}"
            samples_ch.view { sampleID, forward, reverse ->
                "${color_red}Sample: ${color_green}$sampleID${color_red} | Forward: ${color_green}$forward${color_red} | Reverse: ${color_green}$reverse${color_reset}"
            }

        /*
            INSPECT BBDD FOR INTERMEDIATE FILES (i.e. sample has been previously analyzed)
        */

        // Check if the genome has previously been analyzed
            FILE_CHECK(samples_ch)

            // After the FILE_CHECK process
            verified_samples_ch = FILE_CHECK.out.sample_paths
                .collectFile(name: 'all_sample_paths.txt', newLine: true, storeDir: params.outdir)
                .ifEmpty { file("${params.outdir}/empty_all_sample_paths.txt") }

            // Parse the samples into the desired tuple structure
            comp_samples_ch = verified_samples_ch
                .splitCsv()
                .map { row -> 
                    log.debug "DEBUG - Processing sample row: $row"
                    if (row.size() == 10) {
                        def (sampleID, forward, reverse, mtbseq_class, mtbseq_stats, mtbseq_pos, mtbseq_vars, tbdb_out, who_out, mtbseq_vcf) = row
                        tuple(
                            sampleID,
                            forward ? file(forward.trim()) : [],
                            reverse ? file(reverse.trim()) : [],
                            mtbseq_class ? file(mtbseq_class.trim()) : [],
                            mtbseq_stats ? file(mtbseq_stats.trim()) : [],
                            mtbseq_pos ? file(mtbseq_pos.trim()) : [],
                            mtbseq_vars ? file(mtbseq_vars.trim()) : [],
                            tbdb_out ? file(tbdb_out.trim()) : [],
                            who_out ? file(who_out.trim()) : [],
                            mtbseq_vcf ? file(mtbseq_vcf.trim()) : []
                        )   
                    } else {
                        log.warn "Error with channel: $row"
                        null
                    }
                }
                .filter { it != null }

            // Demonstrate the content of the channel
            comp_samples_ch.view { sample -> "Sample: $sample" }


        /*
            SINGLE sample analysis
        */

        // Call the SINGLE_WORKFLOW only for samples missing necessary files
        /// this should also update the channel to contain the paths to any missing single-analysis 
        /// results required for the pairwise comparison
            SINGLE_WF( comp_samples_ch,
                    )
                // Demonstrate the content of the channel
                SINGLE_WF.out.single_updated_samples_ch.view { sample -> "Sample: $sample" }

        /*
            PAIRWISE sample analysis
        

        // Call the PAIRWISE_WF con
            PAIRWISE_WF( SINGLE_WF.out.single_updated_samples_ch,
                            params.runID
                        )
*/
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