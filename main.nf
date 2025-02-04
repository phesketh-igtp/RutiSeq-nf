#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FILE_CHECK }                  from './modules/local/file-checks/main.nf'
include { SINGLE_WF }                   from './workflows/single_wf.nf'
include { PAIRWISE_WF }                 from './workflows/pairwise_wf.nf'
include { NEGATIVE_CONTROL_WF }         from './workflows/negative_ctrl_wf.nf'
include { SUMMARY_WF }                  from './workflows/summary_wf.nf'
//include { BARCODING_WF }                from './workflows/barcoding_wf.nf'

/* 
    Help Message
*/
def helpMessage() {
    log.info"""
    Usage:

    Mandatory arguments:
        --samplesheet           [CSV]   Path to input data (must be surrounded with quotes)
        --outdir                [path]  The output directory where the results will be saved

    Optional arguments:
        --metadata              [CSV]   Metadata file containing the sampleID,sampling_data;loc data.

    Additional parameters:

        MTBseq optional arguments
            --mtbseq_minbqual   [num]   Defines minimum positional mapping quality during variant calling (default: 20).
            --mtbseq_mincovf    [num]   Defines minimum forward read coverage for a putative variant position (default: 4).
            --mtbseq_mincovr    [num]   Defines minimum reverse read coverage for a putative variant position (default: 4).
            --mtbseq_minphred20 [num]   Defines the minimum number of reads having a phred score above or equal 20 to be considered as a putative variant (default: 4).
            --mtbseq_minfreq    [num]   Defines minimum allele frequency for majority allele (default: 75).
            --mtbseq_unambig    [num]   Defines minimum percentage of samples having unambigous base call in TBamend analysi (default: 95).
            --mtbseq_window     [num]   Defines window for SNP cluster look up. Reduces putative false positives in TBamend (default: 10).
            --mtbseq_distance   [num]   Defines SNP distance for the single linkage clustering in TBgroups. Tuple containing a range of values (default: "[5, 10, 15]").

        TBProfiler optional arguments
            --tbprof_           
            --tbprof_

        IQ-Tree optional arguments
            --iqtree_bootstraps [num]   Defines the number of bootstraps used by IQ-Tree for variant positions phylogeny (default: 1000).
            --iqtree_model      [chr]   Defines the maximum-likelihood model used by ID-Tree for variant position phylogeny (default: GTR+G4).

    """
}

/* 
    MAIN WORKFLOW
*/

workflow {
    
    def color_purple = '\u001B[35m'
    def color_green  = '\u001B[32m'
    def color_red    = '\u001B[31m'
    def color_reset  = '\u001B[0m'
    def color_cyan   = '\u001B[36m'

    log.info """
    ${color_cyan}
    ╔════════════════════════════════════════════════════════════════════════╗
    ║  ██████╗ ██╗   ██╗████████╗██╗███████╗███████╗ ██████╗                 ║
    ║  ██╔══██╗██║   ██║╚══██╔══╝██║██╔════╝██╔════╝██╔═══██╗                ║
    ║  ██████╔╝██║   ██║   ██║   ██║███████╗█████╗  ██║   ██║                ║
    ║  ██╔══██╗██║   ██║   ██║   ██║╚════██║██╔══╝  ██║▄▄ ██║                ║
    ║  ██║  ██║╚██████╔╝   ██║   ██║███████║███████╗╚██████╔╝                ║
    ║  ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚═╝╚══════╝╚══════╝ ╚══▀▀═╝  v.1.0.0-beta   ║
    ║ ${color_green}Pre-release development version${color_cyan}                                        ║    
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
        ······································································································
            CREATION OF CHANNELS
                The section creates the samples_ch and the controls_ch from the samplesheet. First the 
                sample sheet is imported and split by the 'type' column into sample or control, and these
                two sets of samples are directed into seperate workflows.
        ······································································································
        */

            // Create channel from sample sheet
                Channel
                    .fromPath(params.samplesheet)
                    .ifEmpty { error "Sample sheet file '${params.samplesheet}' not found or empty" }
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                        def requiredColumns = ['originalID', 'sampleID', 'forward_path', 'reverse_path', 'type']
                        def missingColumns = requiredColumns.findAll { !row.containsKey(it) }
                        if (missingColumns) {
                            error "Missing required column(s) in samplesheet: ${missingColumns.join(', ')}"
                        }
                        
                        // Check for empty paths
                        if (!row.forward_path.trim() || !row.reverse_path.trim()) {
                            error "Empty file path found for sample ${row.sampleID}. Both forward and reverse paths must be provided."
                        }
                        
                        // Use the file function with error checking
                        def forwardFile = file(row.forward_path.trim(), checkIfExists: true)
                        def reverseFile = file(row.reverse_path.trim(), checkIfExists: true)
                        
                        tuple(row.sampleID.trim(), 
                            forwardFile, 
                            reverseFile, 
                            row.type.trim()
                        )
                    }
                    .branch {
                        sample: it[3] == 'sample'
                        control: it[3] == 'control'
                    }
                    .set { branched_samples_by_type }

            // Remove the 'type' from the tuples and ensure only 3 elements
            samples_ch = branched_samples_by_type.sample.map { it -> 
                tuple(it[0], it[1], it[2]) // keep only the sampleID, forward and reverse reads
            }

            controls_ch = branched_samples_by_type.control.map { it -> 
                tuple(it[0], it[1], it[2]) // keep only the sampleID, forward and reverse reads
            }

            // Report the samples part of the samplesheet
                samples_ch.view { sampleID, forward, reverse ->
                    "${color_cyan}Sample: ${color_green}$sampleID${color_reset} | ${color_cyan}Forward: ${color_green}$forward${color_reset} | ${color_cyan}Reverse: ${color_green}$reverse${color_reset}"
                }

             // Report the controls part of the samplesheet
                controls_ch.view { sampleID, forward, reverse ->
                    "${color_red}Control: ${color_green}$sampleID${color_reset} | ${color_red}Forward: ${color_green}$forward${color_reset} | ${color_red}Reverse: ${color_green}$reverse${color_reset}"
                }

        /*
        ······································································································
            NEGATIVE CONTROL WORKFLOW (NEGATIVE_CONTROL_WF)

                - From the controls_ch, the samples are taxonomically classified with Kaiju
                - Taxonomically classified sample reads and produces a summary of the reads
        ······································································································
        */

        // Call the workflow
        // TODO: need to figure out if this is working as intended and correct the channel to not have that empty index [4]
            ///NEGATIVE_CONTROL_WF(controls_ch)

        /*
        ······································································································
            INSPECT BBDD FOR INTERMEDIATE FILES
                - Inspects the BBDD for the sampleID and SINGLE_WF outputs and creates a channel containins paths
        ······································································································
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

            /*
            // DEBUG:: Demonstrate the content of the channel
            comp_samples_ch.view { sample -> "Sample: $sample" }
            */

        /*
        ······································································································
            SINGLE SAMPLE ANALYSIS (SINGLE_WF):
                
                - Taxonomically classified sample reads and produces a summary read stats
                - Paritions MTBc reads for downstream analysis
                - Performs TB-Profiler analysis on MTBc reads (using both WHO and TBDB databases)
                - Performs  MTBSeq analysis on MTBc reads
                - Creates correctly formatted VCF files using MTBseq mpileup files
        ······································································································
        */

            SINGLE_WF( comp_samples_ch )
                    
                // DEBUG: Demonstrate the content of the channel
                ///     SINGLE_WF.out.single_updated_samples_ch.view { sample -> "Sample: $sampleID" }

            // prepare the channel for the pairwise analysis
            pairwise_samples_ch = SINGLE_WF.out.single_updated_samples_ch

            // pairwise_samples_ch.subscribe { println "Debug - pairwise_samples_ch: $it" }

        /*
        ······································································································
            PAIRWISE SAMPLE ANALYSIS (PAIRWISE_WF):

                - PAIRWISE_WF performs comparative analysis of
                - Performs pairwise SNP clustering
                - Creates concatenated variant positions alignments
                - Produces ML phylogenetic trees (with reference and ancestral reconstruction)
                - Produces molecular clock, with ML trees (with reference and ancestral reconstruction)
        ······································································································
        */

            // Structure of the channel : pairwise_samples_ch
            /// [0] sampleID        [1] forward         [2] reverse     [3] mtbseq_class    [4] mtbseq_stats
            /// [5] mtbseq_pos      [6] mtbseq_vars     [7] tbdb_out    [8] who_out         [9] mtbseq_vcf

            // filter channels of just the necessary output files contained within the tuple (by calling the index)
                mtbseq_class_files  =   pairwise_samples_ch.map { it -> it[3] ?: null }
                mtbseq_stats_files  =   pairwise_samples_ch.map { it -> it[4] ?: null }
                tbdb_out_files      =   pairwise_samples_ch.map { it -> it[7] ?: null }
                who_out_files       =   pairwise_samples_ch.map { it -> it[8] ?: null }

                // make the channels
                mtbseq_stats_ch     =   mtbseq_stats_files.collect()
                mtbseq_class_ch     =   mtbseq_class_files.collect()
                tbdb_out_ch         =   tbdb_out_files.collect()
                who_out_ch          =   who_out_files.collect()

            PAIRWISE_WF( params.runID,
                            mtbseq_stats_ch,
                            mtbseq_class_ch,
                            tbdb_out_ch,
                            who_out_ch
                        )

        /*
        ······································································································
            SUMMARY WORKFLOW (SUMMARU_WF):

                - Produces the EXCEL summary tables
                - Visualise phylogenetic trees 
                - Generate MJN files for visualisation in PopArt
        ······································································································
        */

            SUMMARY_WF( params.runID,
                        PAIRWISE_WF.out.pairwise_clusters,
                        PAIRWISE_WF.out.analysis_summary,
                        PAIRWISE_WF.out.who_resistance,
                        PAIRWISE_WF.out.tbdb_resistance,
                        PAIRWISE_WF.out.phylogeny_plotting_ch
                    )

        /*
        ······································································································
            BARCODING ANALYSIS (BARCODING_WF)
                Perform barcoding analysis of the VCF files generated from the single workflow.
                This analysis has a much lower priority
        ······································································································
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