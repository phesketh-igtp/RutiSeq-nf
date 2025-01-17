process PLOT_PHYLOGENY {

    publishDir "${params.outdir}/bbdd/results/phylogeny", mode: 'copy'

    input:
        path(pairwise_clusters)
        tuple val(lineage), path(contree),
                            path(timetree),
                            path(alignments)
        path(metadata)

    output:
        path("phylogeny/*"), emit: phylogeny_dir

    script:
        """
        
        Rscript ${params.r_script_dir}/plot_ML-phylogeny.R  \\
                --contree ${contree} \\
                --timetree ${timetree} \\
                --alignment ${alignments} \\
                --clusters ${pairwise_clusters}

        Rscript ${params.r_script_dir}/plot_ML-phylogeny.R  \\
                --contree ${contree} \\
                --timetree ${timetree} \\
                --alignment ${alignments} \\
                --clusters ${pairwise_clusters}

        """

}