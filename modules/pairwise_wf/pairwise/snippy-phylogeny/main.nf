process SNIPPY_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process performs the concatenated variable region phylogeny analysis
        for concatenated SNP alignments outputs from the MTBSeq pipeline. 
        It uses the MAFFT and IQ-Tree software to perform the phylogeny analysis.
    @changelog
        v1.0.0-2025-11-04: Initial version
        v1.0.1-2025-11-27: Updated to use variant and invariant site alignments from Snippy-core,
                            using -fconst in IQ-Tree for better phylogeny accuracy 
                            (just a linear scaling of the branch lengths.)
*/

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        tuple path(aln_variant_sites),
            path(aln_invariant_sites)

    output:
        path("core_ml*")

    script:

    """
    #--------------------------------------------------------------------------#
    ## Perform the alignments and phylogeny
    #--------------------------------------------------------------------------#

    # Perform phylogeny
        iqtree \\
        -s ${aln_variant_sites}  \\
        -fconst ${aln_invariant_sites} \\
        -m ${params.iqtree_model} \\
        -T AUTO \\
        -ntmax ${params.cpus} \\
        -B ${params.iqtree_bootstraps} \\
        --prefix core_ml
    """

}