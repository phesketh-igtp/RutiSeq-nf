process COMPILE_SEQUENCING_STATS {

/*
    In this module theWhoverall sequencing statistics for the BBDD are
        calculated
        TODO: Fix the Rscript compile-sequencing-statistics.R as it is generating
            and incorrect pairwise_analysis.list.csv
*/

    conda params.r_stats_env

    publishDir "${params.outdir}/bbdd/results/", mode: 'copy'

    input:
        val(runID)
        path(tbdb_results)
        path(who_results)
        path(mtbseq_compiled_strains)
        path(mtbseq_compiled_map_stats)
        path(lineage_fractions)
        val(sampleID_list)

    output:
        path("${runID}.sequencing_summary.csv")

        path("sequencing_summary.csv"),      emit: analysis_summary
        path("who_resistance_summary.csv"),  emit: who_resistance
        path("tbdb_resistance_summary.csv"), emit: tbdb_resistance

        path("pairwise_analysis.list.csv"),  emit: pairwise_analysis_list

    script:

    def additional_args = task.ext.compile_sequencing_stats ?: ''

    """
    # Create the lienage fraction strings
        Rscript ${params.r_script_dir}/tbprofiler_lineage_fractions.R \\
                        --tbprofiler    ${tbdb_results} \\
                        --lineages      lineages.fractions.txt

    # Generate summary statistics and create the sampleID,lineage df for
    ## creating into a channel TODO: need to fix this script in generating the output for tuplec creation
        # the production of the pairwise_analysis.list.csv doest work
        Rscript ${params.r_script_dir}/compile-sequencing-statistics.R \\
                    --minimum_coverage ${params.mtbseq_min_depth} \\
                    --runID ${runID} \\
                    --dictionary_path ${params.r_script_dir}

    # Convert the list of sample IDs to a format suitable for grep
        echo '${sampleID_list.join("\n")}' > run_sample_ids.txt

    # Seperate out the genomes from this run into their own results file
        grep -f run_sample_ids.txt ${runID}.sequencing_summary.csv > tmp.${runID}.sequencing_summary.csv
        mv tmp.${runID}.sequencing_summary.csv ${runID}.sequencing_summary.csv

    # Create the file to go to the tuple seperation
    Rscript -e 'library(tidyverse)
                df <- read_delim("Mapping_and_Variant_Statistics.tab", delim = "\t", col_names = FALSE) |> 
                    distinct() |> filter(X5 >= ${params.mtbseq_min_reads} & X16 >= ${params.mtbseq_min_cov} & X19 >= ${params.mtbseq_min_depth}) |>
                    select(X4) |> distinct()
                write.csv(df, "min.qual.genomes", quote = FALSE, row.names = FALSE)
                '
        sed 's/\t/,/g' tbdb-tbprofiler.txt | cut -d ',' -f1,2,3 > tmp.pairwise_analysis.list.csv
            # remove any ';' which is used in the mixed lienages
        grep -f min.qual.genomes tmp.pairwise_analysis.list.csv | grep -v ';' > pairwise_analysis.list.csv
        touch pairwise_analysis.list.csv

    """
}
