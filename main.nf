#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SINGLE_GENOME_ANALYSIS }  from './workflows/single_genome_analysis.nf'
//include { PAIRIWISE_GENOME_ANALYSIS }      from './workflows/pairwise_analysis.nf'

workflow {
    def color_purple = '\u001B[35m'
    def color_green = '\u001B[32m'
    def color_red = '\u001B[31m'
    def color_reset = '\u001B[0m'

    log.info """
    ${color_purple}
    ╔════════════════════════════════════════════════════════════════════════════╗
    ║                                                                            ║
    ║     ██████╗ ██╗   ██╗████████╗██╗███████╗███████╗ ██████╗                  ║
    ║     ██╔══██╗██║   ██║╚══██╔══╝██║██╔════╝██╔════╝██╔═══██╗                 ║
    ║     ██████╔╝██║   ██║   ██║   ██║███████╗█████╗  ██║   ██║                 ║
    ║     ██╔══██╗██║   ██║   ██║   ██║╚════██║██╔══╝  ██║▄▄ ██║                 ║
    ║     ██║  ██║╚██████╔╝   ██║   ██║███████║███████╗╚██████╔╝                 ║
    ║     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚═╝╚══════╝╚══════╝ ╚══▀▀═╝  v.1.0.0-alpha   ║
    ║                                                                            ║
    ╚════════════════════════════════════════════════════════════════════════════╝
    ${color_reset}
    """

    // Create channel from sample sheet
    if (params.samplesheet == null) {
        error "Please provide a samplesheet CSV file with --samplesheet"
    }

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

    // Call the subworkflow

    SINGLE_GENOME_ANALYSIS(
        samples_ch,
        file(params.kaiju_names),
        file(params.kaiju_nodes),
        file(params.kaiju_fmi),
        file(params.tbprofiler_db)
    )

    // You can now use the outputs from the subworkflow if needed
    SINGLE_GENOME_ANALYSIS.out.qc_results.view()

    // Combine all the MTBSeq results into a single channel
    mtbseq_results = SINGLE_GENOME_ANALYSIS.out.mtbseq_bam_dir
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_bam)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_bam_index)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_bamlog)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_position_tables_dir)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_position_tables)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_classification_dir)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_classification)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_statistics_dir)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_statistics)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_called_dir)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_position_variants)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_mpileup_dir)
                .combine(SINGLE_GENOME_ANALYSIS.out.mtbseq_mpileup)



}