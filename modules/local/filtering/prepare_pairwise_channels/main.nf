process PREPARE_PAIRWISE_CHANNELS {

    conda params.r_stats_env

    publishDir "${params.outdir}/results/", mode: 'copy'
    
    input:
        val runID
        path pairwise_analysis_list
        val sampleID_list

    output:
        path "final.lineage_samples_tuple.csv", emit: lineage_sample_tuple
        path "final.skipped-lineages_tuple.csv", emit: skipped_lineages_tuple

    script:

        """
        # Get the list of sampleIDs from this analysis run
            
            echo '${sampleID_list.join("\n")}' | sort | uniq > run_sample_ids.txt

            echo '${params.lineage_pairwise_sub.join('\n')}' > selected_sub-lineage_split.list
            echo '${params.lineage_pairwise_main.join('\n')}' > selected_main-lineage_split.list

        # Run the script to generate pairwise analysis tuples
            Rscript ${params.r_script_dir}/create_pairwise_analaysi_tuple.R \\
                1>>.command.out \\
                2>>.command.err || true # maybe?

        """

}