process PREPARE_PAIRWISE_CHANNELS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description:
        In this module creates the pairwise analysis tuples from the lineage_samples_paths.csv
        and the lineage_pairwise_sub and lineage_pairwise_main lists.
        The output is a tuple of the form (lineage, sampleID) for each sampleID in the analysis.
*/

    conda params.r_stats_env

    publishDir "${params.outDir}/results/", mode: 'copy'
    
    input:
        val(runID)
        path(pairwise_analysis_list)
        val(sampleID_list)

    output:
        path "final.lineage_samples_tuple.csv", emit: lineage_sample_tuple
        path "final.skipped-lineages_tuple.csv", emit: skipped_lineages_tuple

    script:

        """

        if [[ ${params.pairwise_lv == "sub"} ]]; then
            echo "Pairwise analysis at sub-lineage level"

                # Get the list of sampleIDs from this analysis run
                    
                    echo '${sampleID_list.join("\n")}' | sort | uniq > run_sample_ids.txt

                    echo '${params.lineage_pairwise_sub.join('\n')}' > selected_sub-lineage_split.list
                    echo '${params.lineage_pairwise_main.join('\n')}' > selected_main-lineage_split.list

                # Run the script to generate pairwise analysis tuples
                    Rscript ${params.r_script_dir}/create_pairwise_analysis_tuple_sub.R \\
                        1>>.command.out \\
                        2>>.command.err || true # i think this helps (?)

                # remove headers
                sed '/^lineage,SampleID/d' final.lineage_samples_tuple.csv | sort > tmp.final.lineage_samples_tuple.csv
                    mv tmp.final.lineage_samples_tuple.csv final.lineage_samples_tuple.csv

                sed '/^lineage,SampleID/d' final.skipped-lineages_tuple.csv | sort > tmp.final.skipped-lineages_tuple.csv
                    mv tmp.final.skipped-lineages_tuple.csv final.skipped-lineages_tuple.csv

        elif [[ ${params.pairwise_lv == "main"} ]]; then
            echo "Pairwise analysis at main-lineage level"

                # Get the list of sampleIDs from this analysis run
                    
                    echo '${sampleID_list.join("\n")}' | sort | uniq > run_sample_ids.txt

                    echo '${params.lineage_pairwise_main.join('\n')}' > selected_main-lineage_split.list

                # Run the script to generate pairwise analysis tuples
                    Rscript ${params.r_script_dir}/create_pairwise_analysis_tuple_main.R \\
                        1>>.command.out \\
                        2>>.command.err || true # i think this helps (?)

                # remove headers
                sed '/^lineage,SampleID/d' final.lineage_samples_tuple.csv | sort > tmp.final.lineage_samples_tuple.csv
                    mv tmp.final.lineage_samples_tuple.csv final.lineage_samples_tuple.csv

                sed '/^lineage,SampleID/d' final.skipped-lineages_tuple.csv | sort > tmp.final.skipped-lineages_tuple.csv
                    mv tmp.final.skipped-lineages_tuple.csv final.skipped-lineages_tuple.csv

        elif [[ ${params.pairwise_lv == "none"} ]]; then

                # Get the list of sampleIDs from this analysis run and append with a 
                ## 'no-split' denotion for the lineage, to indicate that no split is performed
                
                    echo '${sampleID_list.join("\n")}' | sort | uniq > run_sample_ids.txt
                    sed -i 's/^/no-split,/g' run_sample_ids.txt > final.lineage_samples_tuple.csv

        else
            echo "Invalid pairwise level specified: ${params.pairwise_lv}"
            exit 1
        fi


        """

}