#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.version = '1.2.0-beta'

include { PREPARE_SAMPLES_WF } from './workflows/prepare-samples_wf.nf'
include { DATABASE_ACCESS_WF } from './workflows/database.access_wf.nf'
include { CONTROL_WF }         from './workflows/control_wf.nf'
include { SINGLE_WF }          from './workflows/single_wf.nf'
include { PAIRWISE_WF }        from './workflows/pairwise_wf.nf'
include { SUMMARY_WF }         from './workflows/summary_wf.nf'

/* 
    Help Message
*/

def helpMessage() {

    log.info """
    Usage:

    Mandatory arguments:
        --samplesheet           [CSV]   Path to input data (must be surrounded with quotes)
        --outDir                [path]  The output directory where the results will be saved
        --workDir               [path]  The temporary work directory for intermediate files (can be deleted when 
                                            analysis is complete to recovered storage space)

    Optional arguments:
        --metadata              [CSV]   Metadata file containing the sampleID,sampling_data;loc data.
        --help                          Print this help message.

    Additional parameters:

        MTBseq optional arguments
            --mtbseq_minbqual   [num]   Defines minimum positional mapping quality during variant calling (default: 20).
            --mtbseq_mincovf    [num]   Defines minimum forward read coverage for a putative variant position (default: 4).
            --mtbseq_mincovr    [num]   Defines minimum reverse read coverage for a putative variant position (default: 4).
            --mtbseq_minphred20 [num]   Defines the minimum number of reads having a Phred score above or equal 20 to be 
                                            considered as a putative variant (default: 4).
            --mtbseq_minfreq    [num]   Defines minimum allele frequency for majority allele (default: 75).
            --mtbseq_unambig    [num]   Defines minimum percentage of samples having unambiguous base call in TBamend 
                                            analysis (default: 95).
            --mtbseq_window     [num]   Defines window for SNP cluster look up. Reduces putative false positives in 
                                            TBamend (default: 10).
            --mtbseq_distance   [num]   Defines SNP distance for the single linkage clustering in TBgroups. Tuple 
                                            containing a range of values (default: "[5, 10, 15]").

        IQ-Tree optional arguments
            --iqtree_bootstraps [num]   Defines the number of bootstraps used by IQ-Tree for variant positions 
                                            phylogeny (default: 1000).
            --iqtree_model      [chr]   Defines the maximum-likelihood model used by ID-Tree for variant position 
                                            phylogeny (default: GTR+G4).
    """
    exit 0
}

/* 
    MAIN WORKFLOW
*/

workflow {

    //def color_purple = '\u001B[35m'
    def green   = '\u001B[32m'
    def red     = '\u001B[31m'
    def no_col  = '\u001B[0m'
    def cyan    = '\u001B[36m'
    def purple  = '\u001B[35m'

    if (params.help) { helpMessage() }

    def missingParams = []
        if (params.samplesheet == null) missingParams << "samplesheet"
        if (params.runID == null) missingParams << "runID"
        if (params.outDir == null) missingParams << "outDir"
        if (params.workDir == null) missingParams << "workDir"

        if (missingParams.size() > 0) {
            error "The following required parameters are missing: ${missingParams.join(', ')}. Please provide them with the appropriate flags."
            helpMessage()
        }

        /*
            DEFINE INPUT ARGUMENTS: expected argument to be provided at time of running 
                nextflow at CLI
            nextflow run main.nf \
                --samplesheet /path/to/sample-sheet
                --runID [a-zA-Z0-9]
                --workflow [full, single, pairwise, summary, barcoding]
        */

        if (params.samplesheet == null) { error "Please provide a samplesheet CSV file with --samplesheet (csv)"; helpMessage() }
        if (params.runID == null) { error "Please provide a runID file with --runID (chr)"; helpMessage() }
        if (params.outDir == null) { error "Please provide a results/database directory for the RutiSeq db (location where new or past results will be) with --outDir (path)"; helpMessage() }
        if (params.workDir == null) { error "Please provide a work directory for the temporary intermediate files --workDir (path)"; helpMessage() }

        /*
        ······································································································
            CREATION OF CHANNELS
                The section creates the samples_ch and the controls_ch from the samplesheet. First the 
                sample sheet is imported and split by the 'type' column into sample or control, and these
                two sets of samples are directed into separate workflows.
        ······································································································
        */

        // Call the subworkflow
            PREPARE_SAMPLES_WF(
                params.samplesheet,
                )
            
            // DEBUG the outputs
            //PREPARE_SAMPLES.out.samples.view { "Sample output: ${it}" }
            //PREPARE_SAMPLES.out.controls.view { "Control output: ${it}" }
            //PREPARE_SAMPLES_WF.out.all_samples.view { "All samples output: ${it}" }

        /*
        ······································································································
            DATABASE_ACCESS (DATABASE_ACCESS_WF):
                - Download following databaset
                    - Sylph GDR220
                    - TB-Profiler (tbdb and who)
        ······································································································
        */

            DATABASE_ACCESS_WF(
                            params.runID
                            )

        /*
        ······································································································
            CONTROLS CHECKS (CONTROL_WF):
                - Taxonomically classifies control samples, and confirms no reads belonging to the 
                    Mycobacterium species specified
        ······································································································
        */

            CONTROL_WF( 
                    params.samplesheet,
                    DATABASE_ACCESS_WF.out.sylph_db
                    )


        /*
        ······································································································
            SINGLE SAMPLE ANALYSIS (SINGLE_WF):
                - Taxonomically classified sample reads and produces a summary read stats
                - Partitions MTBc reads for downstream analysis
                - Performs TB-Profiler analysis on MTBc reads (using both WHO and TBDB databases)
                - Performs  MTBSeq analysis on MTBc reads
                - Creates correctly formatted VCF files using MTBseq mpileup files
        ······································································································
        */

            SINGLE_WF( 
                    PREPARE_SAMPLES_WF.out.all_samples,
                    DATABASE_ACCESS_WF.out.tbprofiler_db
                    )
                    
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
            /// [5] mtbseq_pos      [6] mtbseq_vars     [7] tbdb_out    [8] who_out         [9] snippy_vcf
            /// filter channels of just the necessary output files contained within the tuple (by calling the index)
                sampleID_dump       =   pairwise_samples_ch.map { it -> it[0] ?: null }
                sampleID_list       =   sampleID_dump.collect()
                //  sampleID_list.view()

            PAIRWISE_WF( 
                        sampleID_list
                        )

        /*
        ······································································································
            SUMMARY WORKFLOW (SUMMARY_WF):
                - Produces the EXCEL summary tables
                - Visualise phylogenetic trees 
                - Generate MJN files for visualisation in PopArt
        ······································································································
        */

            SUMMARY_WF(
                        PAIRWISE_WF.out.processed_clusters,
                        PAIRWISE_WF.out.unprocessed_clusters,
                        PAIRWISE_WF.out.analysis_summary,
                        PAIRWISE_WF.out.who_resistance,
                        PAIRWISE_WF.out.tbdb_resistance,
                        PAIRWISE_WF.out.phylogeny_plotting_ch,
                        PAIRWISE_WF.out.nexus_creation_ch,
                        PAIRWISE_WF.out.dated_phylogeny_ch,
                        CONTROL_WF.out.sylph_results,
                        CONTROL_WF.out.read_qc_report
                        
                    )


    /* 
    // Report workflow parameters
    */
    log.info """
    ${green}════════════════════════════════════════════════════════════════════════${cyan}
        ██████╗ ██╗   ██╗████████╗██╗███████╗███████╗ ██████╗                 
        ██╔══██╗██║   ██║╚══██╔══╝██║██╔════╝██╔════╝██╔═══██╗                
        ██████╔╝██║   ██║   ██║   ██║███████╗█████╗  ██║   ██║                
        ██╔══██╗██║   ██║   ██║   ██║╚════██║██╔══╝  ██║▄▄ ██║                
        ██║  ██║╚██████╔╝   ██║   ██║███████║███████╗╚██████╔╝                
        ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚═╝╚══════╝╚══════╝ ╚══▀▀═╝                 
        ${green}Pre-release development version${no_col}
    ${green}════════════════════════════════════════════════════════════════════════
    ${cyan}    RutiSeq.nf: Main workflow for the RutiSeq pipeline (${params.version})${no_col}
    ${green}════════════════════════════════════════════════════════════════════════
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${green}PIPELINE PARAMETER SUMMARY${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}Run Configuration:${no_col}
    ${cyan}   Run ID:${no_col} ${red}${params.runID}${no_col}
    ${cyan}   Output Directory:${no_col} ${red}${params.outDir}${no_col}
    ${cyan}   Input Directory:${no_col} ${red}${params.inputDir}${no_col}
    ${cyan}   Work Directory:${no_col} ${red}${workflow.workDir}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}Pairwise Analysis Configuration:${no_col}
    ${cyan}   Pairwise Split:${no_col} ${red}${params.pairwise_split}${no_col}
    ${cyan}   Main Lineages:${no_col} ${red}${params.lineage_pairwise_main.join(', ')}${no_col}
    ${cyan}   Sub Lineages:${no_col} ${red}${params.lineage_pairwise_sub.join(', ')}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}fastp Parameters:${no_col}
    ${cyan}   Min length:${no_col} ${red}${params.fastp_length_required}${no_col}
    ${cyan}   Max Reads:${no_col} ${red}${params.fastp_max_reads}${no_col}
    ${cyan}   Min Reads:${no_col} ${red}${params.fastp_min_reads}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}MTBseq Parameters:${no_col}
    ${cyan}   Reference:${no_col} ${red}${params.mtbseq_reference}${no_col}
    ${cyan}   Min Base Quality:${no_col} ${red}${params.mtbseq_minbqual}${no_col}
    ${cyan}   Min Coverage Forward:${no_col} ${red}${params.mtbseq_mincovf}${no_col}
    ${cyan}   Min Coverage Reverse:${no_col} ${red}${params.mtbseq_mincovr}${no_col}
    ${cyan}   Min Phred20:${no_col} ${red}${params.mtbseq_minphred20}${no_col}
    ${cyan}   Min Frequency:${no_col} ${red}${params.mtbseq_minfreq}${no_col}
    ${cyan}   Unambiguous:${no_col} ${red}${params.mtbseq_unambig}${no_col}
    ${cyan}   Window:${no_col} ${red}${params.mtbseq_window}${no_col}
    ${cyan}   SNP Distances:${no_col} ${red}${params.mtbseq_snp_distance.join(', ')}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.mtbseq_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}Snippy Parameters:${no_col}
    ${cyan}   Reference:${no_col} ${red}${params.snippy_reference}${no_col}
    ${cyan}   Min Coverage:${no_col} ${red}${params.snippy_mincov}${no_col}
    ${cyan}   Min Fraction:${no_col} ${red}${params.snippy_minfrac}${no_col}
    ${cyan}   Map Quality:${no_col} ${red}${params.snippy_mapqual}${no_col}
    ${cyan}   Min Quality:${no_col} ${red}${params.snippy_minqual}${no_col}
    ${cyan}   Base Quality:${no_col} ${red}${params.snippy_basequal}${no_col}
    ${purple}Core Analysis:${no_col} ${red}${params.snippy_core ? 'Yes' : 'No'}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.snippy_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}IQ-Tree Parameters:${no_col}
    ${cyan}   Bootstraps:${no_col} ${red}${params.iqtree_bootstraps}${no_col}
    ${cyan}   Model:${no_col} ${red}${params.iqtree_model}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.snippy_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}Genome Filtering Parameters:${no_col}
    ${cyan}   Min Coverage:${no_col} ${red}${params.filt_min_cov}${no_col}
    ${cyan}   Min Depth:${no_col} ${red}${params.filt_min_depth}${no_col}
    ${cyan}   Min Reads:${no_col} ${red}${params.filt_min_reads}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.snippy_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}ONT Parameters (WIP):${no_col}
    ${cyan}   Min Coverage:${no_col} ${red}${params.ont_filt_min_cov}${no_col}
    ${cyan}   Min Depth:${no_col} ${red}${params.ont_filt_min_depth}${no_col}
    ${cyan}   Min Reads:${no_col} ${red}${params.ont_filt_min_reads}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.snippy_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    ${purple}SKA Parameters:${no_col}
    ${cyan}   K-mer Size:${no_col} ${red}${params.ska_kmer}${no_col}
    ${cyan}   Min Count:${no_col} ${red}${params.ska_min_count}${no_col}
    ${cyan}   Proportion Reads:${no_col} ${red}${params.ska_proportion_reads}${no_col}
    ${cyan}   Quality Filter:${no_col} ${red}${params.ska_qual_filter}${no_col}
    ${cyan}   Min Quality:${no_col} ${red}${params.ska_min_qual}${no_col}
    ${cyan}   Min Frequency:${no_col} ${red}${params.ska_min_freq}${no_col}
    ${cyan}   Build Args:${no_col} ${red}${params.ska_build_args ?: 'none'}${no_col}
    ${cyan}   Distance Args:${no_col} ${red}${params.ska_distance_args ?: 'none'}${no_col}
    ${cyan}   Additional Args:${no_col} ${red}${params.snippy_args ?: 'none'}${no_col}
    ${green}════════════════════════════════════════════════════════════════════════${no_col}
    """
}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-04
    @version: 1.2.0-beta
    @description: 
        This is the main workflow for the RutiSeq-nf pipeline. It is designed to be run with Nextflow and 
        takes a samplesheet as input. The workflow performs the following steps:
            - Update the TBProfiler database
            - Perform negative control analysis
            - Perform single sample analysis
            - Perform pairwise sample analysis
            - Produce summary tables and visualisations
            - Perform barcoding analysis (optional-WIP)
    @changelog
        - 2024-11-01: Initial version
*/