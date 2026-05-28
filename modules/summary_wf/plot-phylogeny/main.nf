process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_phylogeny_env

    publishDir "${params.outDir}/results/${params.runID}/phylogeny/", 
                mode: 'copy',
                overwrite: true

    input:
        tuple val(lineage), 
            path(contree, stageAs: "snp.contree"), 
            path(alignments, stageAs: "snp.aln.fasta")

        path(processed_clusters, stageAs: "clusters.tsv")
        path(unprocessed_clusters, stageAs: "unprocesses_clusters.tsv")

    output:
    // Main outputs
        path("${lineage}.ML-phylogeny.Rdata")
        path("${lineage}_ML-phylogeny.html")
    
    // Timetree channel output
        tuple val(lineage),  
            path(alignments),
            path(processed_clusters),
            path("unprocesses_clusters.tsv"), emit: timetree_ch

    script:
        """
        # Make output directory
        mkdir -p ${params.outDir}/results/${params.runID}/phylogeny/

        # How many genomes are clustered in this lineage
        clustered_genomes=\$(grep '${lineage}' ${unprocessed_clusters} | cut -f4 | sort -u | grep -v '^\$' | wc -l)

        if [[ \$clustered_genomes -gt 0 ]]; then

            echo "This lineage contains clusters - plotting phylogeny with cluster heatmap"
            # Copy quarto script
            cp ${params.scriptDir}/quarto/phylogeny-report.qmd \\
                phylogeny-report.qmd
            
            # Render script
            quarto render phylogeny-report.qmd \\
                -P runID=${params.runID} \\
                -P lineage=${lineage} \\
                --output ${lineage}_ML-phylogeny.html #2>/dev/null
            cp ${lineage}_ML-phylogeny.html ${params.outDir}/results/${params.runID}/phylogeny/

        else
            echo "This lineage contains no clusters - plotting phylogeny without cluster heatmap"
        fi
        """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 2.0.0
@description:
    Plot the main phylogeny for the analysis, using the ML tree and cluster heatmap if clusters are present.
    If no clusters are present, plot the ML tree without the cluster heatmap.
@changelog:
    v1.0.0-2025-04-01: Initial version.
    v1.1.0-2025-10-01: Remove static tree PDF and replaced with quarto report generation.
    v1.1.1-2025-10-28: Corrected the Quarto report, that was not producing the MonoPhyl output.
    v2.0.0-2026-05-18: Removed the seperate Rscript for creating the ggtree objects
*/
