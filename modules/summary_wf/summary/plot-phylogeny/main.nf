process PLOT_MAIN_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Plot the main phylogeny for the analysis, using the ML tree and cluster heatmap if clusters are present.
        If no clusters are present, plot the ML tree without the cluster heatmap.
    @changelog:
        2025-04-01: Initial version.
        2025-10-01: Remove static tree PDF and replaced with quarto report generation.
        2025-10-28: Corrected the Quarto report, that was not producing the MonoPhyl output. 
*/

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
        tuple val(lineage),  
            path("snp.aln.fasta"),
            path("clusters.tsv")
            path("unprocesses_clusters.tsv"), emit: timetree_ch
        
        path("${lineage}.contree.Rdata"), optional: true
        path("${lineage}_phylogeny.html"), optional: true

    script:
        """
        # How many genomes are clustered in this lineage

        clustered_genomes=\$(grep '${lineage}' ${unprocessed_clusters} | grep -v 'singleton' | wc -l)

        if [[ \$clustered_genomes -gt 0 ]]; then

            echo "This lineage contains clusters - plotting phylogeny with cluster heatmap"

            # plot phylogeny with clsuter heatmap
            Rscript ${params.scriptDir}/R/plot_ML-phylogeny.R \\
                    --lineageID ${lineage} \\
                    --rlibrary ${params.scriptDir}/R/

            # Copy the phylogeny report to the output directory
            cp ${params.scriptDir}/quarto/phylogeny-report.qmd phylogeny-report.qmd

            # Render the phylogeny report using Quarto
            DENO_V8_FLAGS="--max-old-space-size=8192" \\
                quarto render phylogeny-report.qmd \\
                -P runID=${params.runID} \\
                -P lineage=${lineage} \\
                -P RData=${lineage}.contree.Rdata \\
                --output ${lineage}_phylogeny.html

            cp ${lineage}.contree.Rdata ${params.outDir}/results/${params.runID}/phylogeny/
            cp ${lineage}_ML-phylogeny.html ${params.outDir}/results/${params.runID}/phylogeny/ 2>/dev/null
            echo -e "quarto render phylogeny-report.qmd -P runID=${params.runID} -P lineage=${lineage} -P RData=${lineage}.contree.Rdata --output ${lineage}_phylogeny.html" > ${params.outDir}/results/${params.runID}/phylogeny/${lineage}.quarto.sh

            #quarto render phylogeny-report.qmd \\
            #    -P runID=${params.runID} \\
            #    -P lineage=${lineage} \\
            #    -P RData=${lineage}.contree.Rdata \\
            #    --output ${lineage}_ML-phylogeny.html
        
        else
            echo "This lineage contains no clusters - plotting phylogeny without cluster heatmap"
        fi
        """
}