process SNIPPY_LINEAGE_CORE_GUBBINS {

    tag "${lineage}"

    conda params.gubbins_env
    
    publishDir "${params.outDir}/db/comparison/snippy/${lineage}/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple val(lineage), 
            path(fasta)

    output:
    // Gubbins outputs
        path("${lineage}-core.nr.masked.gubbins.aln")
        path("${lineage}-core.nr.masked.gubbins.branch_base_reconstruction.embl")
        path("${lineage}-core.nr.masked.gubbins.filtered_polymorphic_sites.fasta")
        path("${lineage}-core.nr.masked.gubbins.filtered_polymorphic_sites.phylip")
        path("${lineage}-core.nr.masked.gubbins.final_tree.tre")
        path("${lineage}-core.nr.masked.gubbins.log")
        path("${lineage}-core.nr.masked.gubbins.node_labelled.final_tree.tre")
        path("${lineage}-core.nr.masked.gubbins.per_branch_statistics.csv")
        path("${lineage}-core.nr.masked.gubbins.recombination_predictions.embl")
        path("${lineage}-core.nr.masked.gubbins.recombination_predictions.gff")
        path("${lineage}-core.nr.masked.gubbins.summary_of_snp_distribution.vcf")
    // Reference
        path("${lineage}-core.ref.fa")

    script:
    """
    bedtools maskfasta \
        -fi ${fasta} \\
        -bed ${params.snippy_masking} \\
        -fo ${lineage}-core.nr.masked.aln

    # Gubbins
    run_gubbins.py ${lineage}-core.nr.masked.aln \\
        --prefix ${lineage}-core.nr.masked.gubbins \\
        --min-window-size ${params.mtbseq_window} \\
        --threads ${task.cpus} \\
        --extensive-search
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-05-06
@version: 1.0.1
@function:
    This process performs the core SNP analysis for multiple samples
        using the Snippy-core tool from the Snippy pipeline.
@details:
@references: 
    https://bitsandbugs.org/2019/11/06/two-easy-ways-to-run-iq-tree-with-the-correct-number-of-constant-sites/
@changelog
    v1.0.0-2026-05-06: Initial version
    v1.0.1-2026-05-08: Added masking with the ${param.snippy_masking}
*/