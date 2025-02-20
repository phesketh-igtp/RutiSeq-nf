process COMPILE_SEQUENCING_STATS2 {

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/", mode: 'copy'

    input:
        val(runID)
        path(pairwise_analysis_list)
        val(sampleID_list)

    output:
        path("final.lineage_samples_tuple.csv"),    emit: lineage_sample_tuple
        path("skipped-lineages_tuple.csv"),         emit: skipped_lineages_tuple

    script:

    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Get the list of sampleIDs from this analysis run

        echo '${sampleID_list.join("\n")}' > run_sample_ids.txt

    # extract the lineages from the params.config file

        echo '${params.lineage_pairwise_sub.join('\n')}' > selected_sub-lineage_split.list
        echo '${params.lineage_pairwise_main.join('\n')}' > selected_main-lineage_split.list

    # Get the main lineages

        while read -r filt_lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv
            grep -E "\${filt_lineage}" pairwise_analysis.list.csv | \\
                                                while IFS=';' read -r sampleID main_lineage sub_lineage; do

                    
                if [[ "\${filt_lineage}" == "\${main_lineage}"* ]]; then

                    # Append the result to the output file
                    echo "\${filt_lineage},\${sampleID}" >> lineage_samples_tuple.main.csv   

                fi
            done
        done < selected_main-lineage_split.list

    # get the sub-lineages

        while read -r filt_lineage; do
            # Use grep to find matching lines from pairwise_analysis.list.csv
            grep -E "\${filt_lineage}" pairwise_analysis.list.csv | \\
                                                while IFS=';' read -r sampleID main_lineage sub_lineage; do

                # Append the result to the output file
                echo "\${filt_lineage},\${sampleID}" >> lineage_samples_tuple.sub.csv   

            done
        done < selected_sub-lineage_split.list

    #··································································································#

    # Merge the two files ensuing no duplicates my removing an main-lineage classifications that are 
    ## already represented by a sub-lineage classifications (i.e. lineage 4)

        cut -d ',' -f2 lineage_samples_tuple.sub.csv > tmp.lineage_samples_tuple.sub.sampleID
        grep -v -f tmp.lineage_samples_tuple.sub.sampleID lineage_samples_tuple.main.csv > lineage_samples_tuple.main-final.csv

        cat lineage_samples_tuple.main-final.csv lineage_samples_tuple.sub.csv > lineage_samples_tuple.1.csv

    #··································································································#

    # Only create a tuple for analysis of the lineages that are represented by the genomes in this analysis

        grep -f run_sample_ids.txt lineage_samples_tuple.1.csv | cut -d ',' -f1 | sort | uniq > lineages.to.keep

        grep -f lineages.to.keep lineage_samples_tuple.1.csv > final.lineage_samples_tuple.csv

        grep -v -f lineages.to.keep lineage_samples_tuple.1.csv > skipped-lineages_tuple.csv

    """
}    