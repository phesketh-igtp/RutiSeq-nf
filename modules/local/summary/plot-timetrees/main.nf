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

    publishDir "${params.outDir}/results/${runID}/phylogeny/", mode: 'copy', overwrite: true

    input:
        val(runID)
        tuple val(lineage), 
                path(timetree),
                path(ancestral_fasta)
        path(pairwise_clusters)

    output:
        path("${lineage}_TimeTree.contree.pdf")
        path("ancestors/*")
        path("nexus.TT.tuple.csv"),         emit: timetree_tuple
        path("${lineage}.time-tree.Rdata")

    script:
    
        """
        mkdir -p ancestors/
        
        Rscript ${params.r_script_dir}/plot_TimeTree-phylogeny.R \\
                --lineageID ${lineage} \\
                --rlibrary ${params.r_script_dir} \\
                1>>.command.out \\
                2>>.command.err || true

        # Isolate the variant positions for each cluster
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f7 \\
            | sort \\
            | uniq \\
            | sed '1!{/^SampleID/d;}' \\
            | sed '1!{/^nX-/d;}' | sed '1!{/^NA-/d;}' > unique.clusters.list

        # Get list of genomes that have less than 5 clusters
        cut -f7 processed_clusters.tsv | awk '{count[\$1]++} END {for (word in count) if (count[word] > 4) print word}' > frequent_values.txt
        
        for clusterID in `cat unique.clusters.list`; do
            echo "${lineage},\${clusterID},${params.outDir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outDir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab,${params.outDir}/results/phylogeny/ancestors/\${clusterID}.ancestor.positions" >> nexus.TT.tuple.csv
        done

        touch nexus.TT.tuple.csv

        grep -f frequent_values.txt nexus.TT.tuple.csv > nexus.TT.tuple.csv.tmp
        mv nexus.TT.tuple.csv.tmp nexus.TT.tuple.csv
        """

}