process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_stats_env

    publishDir "${params.outdir}/results/phylogeny/", mode: 'copy'

    input:
        tuple val(lineage), 
            path(contree, stageAs: "snp.contree"), 
            path(alignments, stageAs: "snp.aln.fasta")
        path(processed_clusters, stageAs: "clusters.tsv")
        path(unprocessed_clusters, stageAs: "unprocesses_clusters.tsv")

    output:
        path("${lineage}_ML.contree.pdf")
        path("${lineage}.contree.Rdata")

    script:

        """

        # How many genomes are clustered in this lineage
            clustered_genomes=\$(grep '${lineage}' ${unprocessed_clusters} | grep -v 'ungrouped' | wc -l)

        if [[ \$clustered_genomes -gt 0 ]]; then
            # plot phylogeny with clsuter heatmap
            Rscript ${params.r_script_dir}/plot_ML-phylogeny.R \\
                    --contree ${contree} \\
                    --clusters ${processed_clusters} \\
                    --lineageID ${lineage} \\
                    --rlibrary ${params.r_script_dir}

        else
            # just plot phylogeny
            Rscript ${params.r_script_dir}/plot_ML-phylogeny.no-clusters.R \\
                    --contree ${contree} \\
                    --lineageID ${lineage} \\ 
                    --rlibrary ${params.r_script_dir}

        fi

        """

}