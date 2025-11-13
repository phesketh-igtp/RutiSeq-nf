process GENERATE_SUMMARY_REPORT {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        Generates a summary report of the analysis results, including
        a summary of the clusters, phylogeny, and variant sites.
        The report is generated in XLSX format and includes tables.
*/

    conda params.r_stats_env

    publishDir "${params.outDir}/results/${runID}/", mode: 'copy'

    input:
        val runID
        path pairwise_clusters_processed
        path sequencing_summary
        path who_resistance
        path tbdb_resistance
        path(sylph_results, stageAs: 'sylph_results.tsv')
        tuple path(html_report), path(warnings, stageAs: 'warnings.out')
        

    output:
        path("${runID}_RutiSeq-results.xlsx")

    script:

        """
        # Correct the who_resistance file
        ## copy the WHO resistance file to the current working directory
            cp ${params.r_script_dir}/dict/tbprofiler_drt_rules.csv .
        ## correct the WHO resistance file
            Rscript ${params.r_script_dir}/rules_based_extended_DRT_assignment.R

        Rscript ${params.r_script_dir}/generate_summary_report.R \\
                --rlibrary ${params.r_script_dir} \\
                --output ${runID}_RutiSeq-results.xlsx
        """
}