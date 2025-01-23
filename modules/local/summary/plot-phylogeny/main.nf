process PLOT_MAIN_PHYLOGENY {

    tag "${lineage}"

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/Phylogeny/", mode: 'copy'

    input:
        tuple val(lineage), path(contree), path(alignments)
        path(pairwise_clusters)

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