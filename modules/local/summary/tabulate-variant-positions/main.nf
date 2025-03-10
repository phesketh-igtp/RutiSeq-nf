process TABULATE_VARIANT_SITES{

    tag "cluster: ${clusterID}"

    conda params.r_stats_env

    publishDir "${params.outdir}/results/${runID}/snps/variants-tab/", mode: 'copy'

    input:
    val(runID)
        tuple val(lineage), 
                val(clusterID), 
                path(snp_alignments),
                path(genomic_possitions),
                path(snp_tab)

    output:
        path("${clusterID}.variant-positions.csv"),        emit: tabular_vars
        path("${clusterID}.variant-positions.counts.csv"), emit: tabular_var_counts


    script:

        """
        Rscript ${params.r_script_dir}/tabulate_variant_sites.R \\
            --cluster   ${clusterID} \\
            --fasta     ${snp_alignments} \\
            --positions ${genomic_possitions} \\
            --H37Rv     ${params.mtbseq_gene_annotations}
        """

}