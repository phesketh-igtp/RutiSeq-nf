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

    script:
    """
    # Gubbins
    run_gubbins.py ${fasta} \\
        --prefix ${lineage}-core.nr.masked.gubbins \\
        --min-window-size ${params.mtbseq_window} \\
        --threads ${task.cpus} \\
        --extensive-search
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-05-22
@version: 1.0.0
@function:
    Run gubbins with reference free alignment from snippy-core
@details:
@changelog
    v1.0.0-2026-05-22: Initial version
*/