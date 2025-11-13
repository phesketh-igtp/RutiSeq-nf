process READ_TAXONOMY_QC_REPORT {

    publishDir "${params.outdir}/RutiSeq/${params.runID}/reports", mode: 'move', overwrite: true

    input:
        path(sylph_output, stageAs: 'sylph_classification.csv')
        path(ska_input, stageAs: 'ska_profiling.csv')
        path(read_stats, stageAs: 'read_stats.csv')

    output:
    tuple
        path("${params.runID}_reads-controls-report.qmd"),
        path("${params.runID}_warnings.out"), emit: reads_taxonomy_qc_report_out
    
    script:
    
    """
    # Copy the quarto document to the work directory
    cp ${params.scriptDir}/quarto/read_taxonomy_qc_report.qmd ${params.runID}_reads-controls-report.qmd

    # Render the report
        quarto render ${params.runID}_reads-controls-report.qmd \
            --to html
    """

}