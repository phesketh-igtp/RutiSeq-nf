process SNIPPY_PHYLOGENY {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        This process performs the concatenated variable region phylogeny analysis
        for concatenated SNP alignments outputs from the MTBSeq pipeline. 
        It uses the MAFFT and IQ-Tree software to perform the phylogeny analysis.
*/

    conda params.phylogeny_env

    publishDir "${params.outDir}/db/comparison/snippy/", mode: 'copy'

    input:
        path(fasta)

    output:
        path("core_ml*")

    script:

    """
    #--------------------------------------------------------------------------#
    ## Perform the alignments and phylogeny
    #--------------------------------------------------------------------------#

    # Perform alignment of sequences 
        mafft --auto --thread ${params.cpus} \\
            ${fasta} \\
            > alignment.fasta

    # Perform phylogeny
        iqtree -s alignment.fasta \\
        -m ${params.iqtree_model} \\
        -T AUTO \\
        -ntmax ${params.cpus} \\
        -B ${params.iqtree_bootstraps} \\
        --prefix core_ml
    """

}