process SNIPPY_PHYLOGENY {

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/snippy/", 
        mode: 'copy',
        overwrite: true

    input:
        path(core_aln_full)

    output:
        path("core_ml*")
        
        tuple path("core.masked.full.filt98.variants.aln"), 
            path("core.masked.full.invariants"),
            path("core_ml.contree"),
            emit: snippy_core_dated_phylogeny

    script:

    """
    # Extract variant sites for phylogeny
        coresnpfilter -e -c 0.98 ${core_aln_full} > core.masked.full.filt98.variants.aln
        coresnpfilter -e -c 0.98 ${core_aln_full} --table core.masked.full.filt98.variants.tab
        coresnpfilter -C ${core_aln_full} > core.masked.full.invariants

    # Perform phylogeny
        iqtree \\
            -s          core.masked.full.filt98.variants.aln  \\
            -fconst     core.masked.full.invariants \\
            -m          ${params.iqtree_model} \\
            -ntmax      ${task.cpus} -T AUTO \\
            -B          ${params.iqtree_bootstraps} \\
            --prefix    core_ml
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