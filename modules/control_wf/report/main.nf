process READ_TAXONOMY_QC_REPORT {

/*
    @author: Poppy J hesketh Best
    @date: 2025-11-17
    @version: v1.0.0
    @description:
        Produces read statistics report from the outputs of Sylph, SKA2, seqkit.
    @changelog:
        v1.0.0-2025-11-17: Functioning module created.
*/

    conda params.r_stats_env

    publishDir "${params.outDir}/db/qc/${params.runID}/", mode: 'copy'

    input:
        path(sylph_output, stageAs: 'sylph_classification.tsv')
        path(ska_input, stageAs: 'ska_distance.tsv')
        path(read_stats, stageAs: 'read_stats.tsv')
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
    tuple path("${params.runID}_reads-controls-report.qmd", optional: true),
        path("${params.runID}_warnings.out", optional: true), emit: html
    
    script:
    
    """
    # Create an empty output incase the file is empty
    > ${params.runID}_warnings.out
    
    # Copy the quarto document to the work directory
    cp ${params.scriptDir}/quarto/reads-controls-report.qmd ${params.runID}_reads-controls-report.qmd

    # Render the report
        quarto render ${params.runID}_reads-controls-report.qmd \
            --to html

    
    """

}