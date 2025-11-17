process PLOT_TIMETREES {

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

    tag "${lineage}"

    conda params.r_phylogeny_env

    publishDir "${params.outDir}/results/${params.runID}/phylogeny/", mode: 'copy', overwrite: true

    input:
        tuple val(lineage), 
                path(timetree),
                path(ancestral_fasta)
        path(pairwise_clusters)

    output:
        path("${lineage}_TimeTree.contree.pdf"), optional: true
        path("ancestors/*"),                     optional: true
        path("nexus.TT.tuple.csv"),              optional: true, emit: timetree_tuple
        path("${lineage}.time-tree.Rdata"),      optional: true

    script:
    
        """
        mkdir -p ancestors/
        
        Rscript ${params.r_script_dir}/plot_TimeTree-phylogeny.R \\
                1>>.command.out \\
                2>>.command.err || true

        # Isolate the variant positions for each cluster
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f3 | sort | uniq \\
            | grep -v 'singleton'  > unique.clusters.list

        # Create the nexus tuple file
        for clusterID in `cat unique.clusters.list`; do
            echo "${lineage},\${clusterID},${params.outDir}/db/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outDir}/db/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab,${params.outDir}/results/phylogeny/ancestors/\${clusterID}.ancestor.positions" >> nexus.TT.tuple.csv
        done

        touch nexus.TT.tuple.csv
        """

}