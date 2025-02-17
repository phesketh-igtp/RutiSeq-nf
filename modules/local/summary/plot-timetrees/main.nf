process PLOT_TIMETREES {

    tag "${lineage}"

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/phylogeny/", mode: 'copy'

    input:
        tuple val(lineage), 
                path(timetree),
                path(ancestral_fasta)
        path(pairwise_clusters)

    output:
        path("${lineage}_TimeTree.contree.pdf")
        path("ancestors/*.fasta")
        path("nexus.TT.tuple.csv"),         emit: timetree_tuple
        path("${lineage}..time-tree.RData")

    script:
    
        """
        mkdir ancestors/
        
        Rscript ${params.r_script_dir}/plot_TimeTree-phylogeny.R \\
                --timetree ${timetree} \\
                --clusters ${pairwise_clusters} \\
                --fasta ${ancestral_fasta} \\
                --lineageID ${lineage} \\
                --rlibrary ${params.r_script_dir}

        # Isolate the variant positions for each cluster
        grep "${lineage}" ${pairwise_clusters} \\
            | cut -f7 \\
            | sort \\
            | uniq \\
            | sed '1!{/^SampleID/d;}' \\
            | sed '1!{/^nX-/d;}' > unique.clusters.list

        for clusterID in `cat unique.clusters.list`; do
            echo "${lineage},\${clusterID},${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta,${params.outdir}/bbdd/mtbseq/pairwise/${lineage}/Amend/${lineage}_joint_cf*_cr*_fr*_ph*_samples*_amended_u${params.mtbseq_unambig}_phylo.tab,${params.outdir}/bbdd/results/phylogeny/ancestors/\${clusterID}.ancestor.positions" >> nexus.TT.tuple.csv
        done

        touch nexus.TT.tuple.csv
        """

}