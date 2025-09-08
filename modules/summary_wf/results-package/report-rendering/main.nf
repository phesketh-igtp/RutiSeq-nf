process RENDER_REPORT {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-08-07
    @version: 1.0.0
    @description:
        This process prepares the data delivery for the final results. It moves the relevant files to the appropriate directories and cleans up any unnecessary files.
        It is called after the summary report generation and before the final data delivery.
    @changelog
        v1.0.0-2025-08-07: Initial version
*/

    tag "${params.runID}"

    publishDir "${params.outDir}/results/${params.runID}/", mode: 'copy'

    input:
        path(analysisSummary)
        path(processedClusters)
        path(whoResistance)
        path(tbdbResistance)

    output:
        path("report.html"), emit: report

    script:

        """
        cp ${params.scriptDir}/quarto/report.qmd .

        quarto render report.qmd --execute \\
            --params runID:${params.runID} \\
            --params analysisSummary:'${analysisSummary}' \\
            --params processedClusters:'${processedClusters}' \\
            --params whoResistance:'${whoResistance}' \\
            --params tbdbResistance:'${tbdbResistance}' \\
            --params samplesheet:'${params.samplesheet}' \\
            --params metadata:'${params.metadata}' \\
            --params minCov:${params.minCov} \\
            --params minCovPerc:${params.minCovPerc} \\
            --params minReads:${params.minReads} \\
            --params iqtree_bootstraps:${params.iqtree_bootstraps} \\
            --params iqtree_model:'${params.iqtree_model}' \\
            --params mtbseq_minbqual:${params.mtbseq_minbqual} \\
            --params mtbseq_mincovf:${params.mtbseq_mincovf} \\
            --params mtbseq_mincovr:${params.mtbseq_mincovr} \\
            --params mtbseq_minphred20:${params.mtbseq_minphred20} \\
            --params mtbseq_minfreq:${params.mtbseq_minfreq} \\
            --params mtbseq_unambig:${params.mtbseq_unambig} \\
            --params mtbseq_window:${params.mtbseq_window}
        """

}