process PLOT_MAIN_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Plot the main phylogeny for the analysis, using the ML tree and cluster heatmap if clusters are present.
        If no clusters are present, plot the ML tree without the cluster heatmap.
*/

    tag "${lineage}"

    conda params.r_phylogeny_env

    publishDir "${params.outDir}/results/${runID}/phylogeny/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage), 
            path(contree, stageAs: "snp.contree"), 
            path(alignments, stageAs: "snp.aln.fasta")
        path(processed_clusters, stageAs: "clusters.tsv")
        path(unprocessed_clusters, stageAs: "unprocesses_clusters.tsv")

    output:
        path("${lineage}.contree.Rdata", optional:true)
        path("${lineage}_phylogeny-report.html", optional:true)

    script:
        """
        # How many genomes are clustered in this lineage

        clustered_genomes=\$(grep '${lineage}' ${unprocessed_clusters} | grep -v 'singleton' | wc -l)

        if [[ \$clustered_genomes -gt 0 ]]; then

            echo "This lineage contains clusters - plotting phylogeny with cluster heatmap"

            # plot phylogeny with clsuter heatmap
            Rscript ${params.r_script_dir}/plot_ML-phylogeny.R \\
                    --lineageID ${lineage} \\
                    --rlibrary ${params.r_script_dir}

            # Copy the phylogeny report to the output directory
            cp ${params.script_dir}/quarto/phylogeny-report.qmd phylogeny-report.qmd

            # Render the phylogeny report using Quarto
            quarto render phylogeny-report.qmd \\
                -P runID=${params.runID} \\
                -P lineage=${lineage} \\
                -P RData=${lineage}.contree.Rdata \\
                --output ${lineage}_ML-phylogeny.html

        else 

            echo "This lineage contains no clusters - plotting phylogeny without cluster heatmap"

            Rscript ${params.r_script_dir}/plot_ML-phylogeny.no-clusters.R \\
                    --lineageID ${lineage}

        fi
        """
        
}