
process SYLPH_CLASSIFICATION_INSPECTION {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0
    @description: 
        Use sylph to classify reads from a sample.
        This process takes the sample ID and the Sylph database as input,
        and outputs a merged Sylph sequence abundance file.
    @changelog
        v1.0.0-2025-04-01: Initial version
*/

    conda params.readQC_env

    publishDir "${params.outDir}/bbdd/read-qc/", mode: 'copy'

    input:
        path(sylph_results)
        path(failed_samples_report)

    output:
        path("${params.runID}.failed-samples.report"), emit: read_qc_res

    script:

    """
    Rscript --vanilla -e '
        library("sylph")
        library("dplyr")
        library("readr")

        # Load the Sylph results
        sylph_results <- read_sylph("${sylph_results}")

        # Identify reads that have mixed taxonomy
        mixed_taxonomy <- sylph_results %>%
            filter(!is.na(taxonomy) & taxonomy != "Unclassified") %>%
            group_by(sample_id) %>%
            summarise(mixed = n_distinct(taxonomy) > 1)

        # Write the results to a file
        write_csv(mixed_taxonomy, "${params.runID}.failed-samples.report")
    """


}