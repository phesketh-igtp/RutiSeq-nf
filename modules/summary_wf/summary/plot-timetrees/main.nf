process PLOT_TIMETREES {

    tag "${lineage}"

    conda params.r_phylogeny_env

    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/", mode: 'copy', overwrite: true

    input:
        tuple val(lineage), 
                path(timetree),
                path(ancestral_fasta)
        path(pairwise_clusters)

    output:
        path("${lineage}_TimeTree.contree.pdf"), optional: true
        path("ancestors/*"),                     optional: true
        path("nexus.TT.tuple.csv"),              optional: true, emit: timetree_tuple

    script:

    def snp_fasta_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
    def snp_tab_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab"

        """
        mkdir -p ancestors/
        
        Rscript ${params.scriptDir}/R/plot_TimeTree-phylogeny.R \\
                1>>.command.out \\
                2>>.command.err || true

        # Isolate the variant positions for each cluster
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f3 | sort | uniq \\
            | grep -v 'singleton'  > unique.clusters.list

        # Create the nexus tuple file
        for clusterID in `cat unique.clusters.list`; do

            echo "${lineage},\${clusterID},${snp_fasta_path},${snp_tab_path},${params.outDir}/db/comparison/mtbseq/${lineage}/ancestors/\${clusterID}.ancestor.positions" \\
                >> nexus.TT.tuple.csv

        done

        touch nexus.TT.tuple.csv
        """

}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0
@description:
    This process generates a time tree for each cluster using the
    R package 'ggtree' and the 'ggplot2' library. The time tree is
    generated from the phylogenetic tree and the ancestral sequences
    for each cluster. The time tree is then plotted using the 'ggtree'
    package and saved as a PDF file.
*/