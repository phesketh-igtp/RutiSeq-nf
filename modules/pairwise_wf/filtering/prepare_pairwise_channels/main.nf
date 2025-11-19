process PREPARE_PAIRWISE_CHANNELS {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 2.0.0
    @description:
        In this module creates the pairwise analysis tuples from the lineage_samples_paths.csv
        and the lineage_pairwise_sub and lineage_pairwise_main lists.
        The output is a tuple of the form (lineage, sampleID) for each sampleID in the analysis.
        There are three options for the pairwise analysis (specified by the params.pairwise_split):
            - sub: pairwise analysis at sub-lineage level
            - main: pairwise analysis at main-lineage level
            - none: pairwise analysis of all samples without lineage split
    @changelog:
        v2.0.0-2025-04-01: Updated to use the new lineage_pairwise_sub and lineage_pairwise_main lists
        v1.0.1-2024-11-01: Added error handling for invalid pairwise level
*/

    conda params.r_stats_env

    publishDir "${params.outDir}/comparison/src/", mode: 'copy'
    
    input:
        path(pairwise_analysis_list)
        val(sampleID_list)

    output:
        path "final.lineage_samples_tuple.csv", emit: lineage_sample_tuple
        path "final.skipped-lineages_tuple.csv", emit: skipped_lineages_tuple

    script:

        """
        # Prapare the lists of sampleIDs and (sub)lineages for possible splits

            echo '${sampleID_list.join("\n")}' | sort | uniq > run_sample_ids.txt
            echo '${params.lineage_pairwise_sub.join('\n')}' | sort | uniq > selected_sub-lineage_split.list
            echo '${params.lineage_pairwise_main.join('\n')}' | sort | uniq > selected_main-lineage_split.list

        if [[ ${params.pairwise_split} == "sub" ]]; then

            echo "Pairwise analysis at sub-lineage level lineage split"

                # Run the script to generate pairwise analysis tuples
                    Rscript ${params.r_scriptDir}/create_pairwise_analysis_tuple_sub.R \\
                        1>>.command.out \\
                        2>>.command.err || true # i think this helps (?)

                # remove headers
                sed '/^lineage,SampleID/d' final.lineage_samples_tuple.csv | \\
                    sed '/^sub_lineage,SampleID/d' | sort > tmp.final.lineage_samples_tuple.csv
                    mv tmp.final.lineage_samples_tuple.csv final.lineage_samples_tuple.csv

                sed '/^lineage,SampleID/d' final.skipped-lineages_tuple.csv | \\
                    sed '/^sub_lineage,SampleID/d' | sort > tmp.final.skipped-lineages_tuple.csv
                    mv tmp.final.skipped-lineages_tuple.csv final.skipped-lineages_tuple.csv

        elif [[ ${params.pairwise_split} == "main" ]]; then

            echo "Pairwise analysis at main-lineage level lineage split"

                # Run the script to generate pairwise analysis tuples
                    Rscript ${params.r_scriptDir}/create_pairwise_analysis_tuple_main.R \\
                        1>>.command.out \\
                        2>>.command.err || true # i think this helps (?)

                # remove headers
                sed '/^lineage,SampleID/d' final.lineage_samples_tuple.csv | \\
                    sed '/^sub_lineage,SampleID/d' | sort > tmp.final.lineage_samples_tuple.csv
                    mv tmp.final.lineage_samples_tuple.csv final.lineage_samples_tuple.csv

                sed '/^lineage,SampleID/d' final.skipped-lineages_tuple.csv | \\
                    sed '/^sub_lineage,SampleID/d' | sort > tmp.final.skipped-lineages_tuple.csv
                    mv tmp.final.skipped-lineages_tuple.csv final.skipped-lineages_tuple.csv

        elif [[ ${params.pairwise_split} == "none" ]]; then

            echo "Pairwise analysis of all samples without lineage split"

                # Get the list of sampleIDs from this analysis run and append with a 
                ## 'no-split' denotion for the lineage, to indicate that no split is performed
                
                    cut -f1 ${pairwise_analysis_list} > run_sample_ids.txt
                    sed 's/^/All,/g' run_sample_ids.txt > final.lineage_samples_tuple.csv
                    
                    touch final.skipped-lineages_tuple.csv

        else

            echo "
            Invalid pairwise level specified: ${params.pairwise_split}. 
            Choose 'sub', 'main', or 'none'., the re-run the workflow with the correct parameter and '-resume-'
            "
            exit 1

        fi


        """

}