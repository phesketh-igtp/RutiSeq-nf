process SNIPPY_LINEAGE_CORE_PHYLOGENY {

    tag "${lineage}"

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/snippy/${lineage}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(lineage),
            path(full_fasta),
            path(core_aln)

    output:
        path("${lineage}.core_ml*")

    script:

    """
    # Perform phylogeny
    iqtree \\
        -s ${core_aln}  \\
        -fconst "\$(snp-sites -C ${full_fasta})" \\
        -m ${params.iqtree_model} \\
        -ntmax ${task.cpus} -T AUTO \\
        -B ${params.iqtree_bootstraps} \\
        --prefix ${lineage}.core_ml
    """
}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process performs the concatenated variable region phylogeny analysis
        for concatenated SNP alignments outputs from the snippy-core pipeline. 
        Bacterial whole genome phylogenetics often involves generating a consensus genome 
            from your data + an appropriate reference genome, then combining these consensus 
            genomes for different isolates, extracting the variant positions and using them 
            to build a phylogenetic tree.
        One problem with this is that the whole genome alignment usually consists primarily 
            of invariant (monomorphic) sites. If you don’t take this into account when 
            generating your phylogeny, the branch lengths of your tree will not be correct, 
            and the topology of parts of the tree with short branches will potentially be affected.
        The most common solution for this is to use ascertainment bias correction. However, 
            there has been some rumblings for a while now that this is not an appropriate method. 
            It should be self-evident that if we throw away information, one should not expect 
            it to be magically recovered. I was recently privy to an interesting discussion of 
            this issue, and thought it would be good to share the outcome of the discussion.
    @changelog
        v1.0.0-2025-11-04: Initial version
        v1.0.1-2025-11-27: Added the concatenation of outgroup sequences to the alignment before phylogeny.
                        Updated to use variant and invariant site alignments from Snippy-core,
                            using -fconst in IQ-Tree for better phylogeny accuracy 
                            (just a linear scaling of the branch lengths.)
                        
*/