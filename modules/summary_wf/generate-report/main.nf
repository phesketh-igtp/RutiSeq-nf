process GENERATE_REPORT {

    publishDir "${params.outDir}/results/${params.runID}/"

    input:
        val runID
        path(processed_clusters)
        path(analysis_summary)
        path(who_resistance)
        path(tbdb_resistance)

    output:
        file("${runID}_report.html")

    script:
        """
        quarto render ${params.script_dir}/quarto/report.qmd \\
            --output ${runID}_report.html \\
            --execute --to html \\
            --execute-env R \\
            -e "analysisSummary=${analysis_summary},processedClusters=${processed_clusters},whoResistance=${who_resistance},tbdbResistance=${tbdb_resistance},reportDir=${params.outDir}/results/${params.runID}/"
        """
}
