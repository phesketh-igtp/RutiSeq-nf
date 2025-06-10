process TABULATE_VARIANT_SITES {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Tabulate variant sites for each cluster. This process takes the output of the
        generate_nexus process and generates a table of variant positions for each cluster.
*/

    tag "cluster: ${clusterID}"

    conda params.r_stats_env

    publishDir "${params.outDir}/results/${runID}/snps/variants-tab/", mode: 'copy'

    input:
    val(runID)
        tuple val(lineage), 
                val(clusterID), 
                path(snp_alignments),
                path(genomic_possitions),
                path(snp_tab)

    output:
        path("${clusterID}.variant-positions.csv"),        emit: tabular_vars, optional: true
        path("${clusterID}.variant-positions.counts.csv"), emit: tabular_var_counts, optional: true


    script:

        """
        Rscript ${params.r_script_dir}/tabulate_variant_sites.R \\
            --cluster   ${clusterID} \\
            --fasta     ${snp_alignments} \\
            --positions ${genomic_possitions} \\
            --H37Rv     ${params.mtbseq_gene_annotations} \\
                1>>.command.out \\
                2>>.command.err || true # prevents stopping the workflow is the
        """

}