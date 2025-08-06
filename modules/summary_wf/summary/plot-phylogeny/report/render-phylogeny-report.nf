process RENDER_PHYLGENY_REPORT {

    input:
        tuple val(lineage),
            file(rdata_file)

    output:
        path "phylogeny-report.html"

    script:

    """
    # Render the phylogeny report using Quarto
    quarto render phylogeny-report.qmd \\
        --to html \\
        -P runID=t${params.runID} \\
        -P lineage=${lineage} \\
        -P RData=${rdata_file} \\
        --output ${lineage}-ML-phylogeny-report.html
    """

}