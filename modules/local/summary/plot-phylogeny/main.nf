process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_stats_env

    publishDir "${params.outdir}/results/${runID}/phylogeny/", mode: 'copy'

    input:
        val(runID)
        tuple val(lineage), 
            path(contree, stageAs: "snp.contree"), 
            path(alignments, stageAs: "snp.aln.fasta")
        path(processed_clusters, stageAs: "clusters.tsv")
        path(unprocessed_clusters, stageAs: "unprocesses_clusters.tsv")

    output:
        path("${lineage}_ML.contree.pdf", optional:true)
        path("${lineage}.contree.Rdata", optional:true)

    script:
        """
        # How many genomes are clustered in this lineage
        clustered_genomes=\$(grep '${lineage}' ${unprocessed_clusters} | grep -v 'ungrouped' | wc -l)

        if [[ \$clustered_genomes -gt 0 ]]; then

            echo "This lineage contains clusters - plotting phylogeny with cluster heatmap"

            # plot phylogeny with clsuter heatmap
            Rscript ${params.r_script_dir}/plot_ML-phylogeny.R \\
                    --lineageID ${lineage} \\
                    --rlibrary ${params.r_script_dir}

        else 

            echo "This lineage contains no clusters - plotting phylogeny without cluster heatmap"

            Rscript ${params.r_script_dir}/plot_ML-phylogeny.no-clusters.R \\
                    --lineageID ${lineage}

        fi
        """
        
}