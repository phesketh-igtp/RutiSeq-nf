process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/Phylogeny", mode: 'copy'

    input:
        path(pairwise_clusters)
        tuple val(lineage), path(contree),
                            path(alignments)

    output:
        path("${lineage}_ML.contree.pdf")

    script:
    
        """
        Rscript ${params.r_script_dir}/plot_ML-phylogeny.R \\
                --contree ${contree} \\
                --clusters ${pairwise_clusters} \\
                --lineageID ${lineage} \\
                --rlibrary ${params.r_script_dir}
        """

}