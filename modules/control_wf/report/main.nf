process READ_TAXONOMY_QC_REPORT {

    conda params.r_stats_env

    publishDir "${params.outDir}/results/${params.runID}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple path(sylph_sequence_abundance, stageAs: 'sylph_sequence_abundance.tsv'),
            path(sylph_relative_abundance, stageAs: 'sylph_relative_abundance.tsv'),
            path(sylph_coverage, stageAs: 'sylph_coverage.tsv')

        path(ska_input, stageAs: 'ska_distance.tsv')
        path(read_stats, stageAs: 'read_stats.tsv')
        path(samplesheet, stageAs: 'samplesheet.csv')

    output:
    tuple path("${params.runID}_reads-controls-report.html", optional: true),
        path("${params.runID}_warnings.out", optional: true), emit: html
    
    script:
    
    """
    if [[ -s "${params.outDir}/results/${params.runID}/${params.runID}_reads-controls-report.html" ]]; then
        echo -e "Already generated, copying from results file"
        cp ${params.outDir}/results/${params.runID}/${params.runID}_reads-controls-report.html .
        echo "Already generated, copying from results file: ${params.outDir}/results/${params.runID}/${params.runID}_reads-controls-report.html" > ${params.runID}_warnings.out
    else
        # Create an empty output incase the file is empty
        > ${params.runID}_warnings.out

        # Copy the quarto document to the work directory
        cp ${params.scriptDir}/quarto/reads-controls-report.qmd ${params.runID}_reads-controls-report.qmd

        # Render the report
        quarto render ${params.runID}_reads-controls-report.qmd --to html
    fi
    """

}


/*
@author: Poppy J hesketh Best
@date: 2025-11-17
@version: v1.0.0
@description:
    Produces read statistics report from the outputs of Sylph, SKA2, seqkit.
@changelog:
    v1.0.0-2025-11-17: Functioning module created.
*/
