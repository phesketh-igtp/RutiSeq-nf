process ASSESS_SAMPLES {

    conda params.r_stats_env

    storeDir "${params.outDir}/db/comparison/src/${params.runID}/"

    input:
        path(pairwise_analysis_list)
        val(sampleID_list)

    output:
        path "final.lineage_samples_tuple.csv", emit: lineage_sample_tuple
        path "final.skipped-lineages_tuple.csv", emit: skipped_lineages_tuple

    script:
        def sub_lineages = params.lineage_pairwise_sub.join('\\n')
        def main_lineages = params.lineage_pairwise_main.join('\\n')

    """
    #!/usr/bin/env Rscript --vanilla

    # Load libraries
    packages <- c("dplyr", "stringr", "readr")

    # Identify missing packages
    missing_pkgs <- packages[!packages %in% installed.packages()[, "Package"]]

    # Install missing packages
    if (length(missing_pkgs) > 0) {
        install.packages(missing_pkgs, dependencies = TRUE)
    }

    # Load all packages
    invisible(lapply(packages, library, character.only = TRUE))

    set.seed(1234)
    options(warn=-1)

    ########################################################################################################
    # SETUP AND PARAMETER HANDLING
    ########################################################################################################

    # Get parameters from Nextflow
    pairwise_split <- "${params.pairwise_split}"
    sampleID_string <- "${sampleID_list}"
    pairwise_analysis_file <- "${pairwise_analysis_list}"

    cat("=== ASSESS_SAMPLES R Script ===\\n")
    cat("Pairwise split parameter:", pairwise_split, "\\n")
    cat("Input file:", pairwise_analysis_file, "\\n")

    ########################################################################################################
    # PREPARE INPUT FILES
    ########################################################################################################

    # Create sample IDs file from comma-separated string
    sample_ids <- unlist(strsplit(sampleID_string, ","))
    sample_ids <- sort(unique(trimws(sample_ids)))
    writeLines(sample_ids, "run_sample_ids.txt")

    # Create lineage lists
    sub_lineages_list <- unlist(strsplit("${sub_lineages}", "\\n"))
    sub_lineages_list <- sort(unique(trimws(sub_lineages_list)))
    writeLines(sub_lineages_list, "selected_sub-lineage_split.list")

    main_lineages_list <- unlist(strsplit("${main_lineages}", "\\n"))
    main_lineages_list <- sort(unique(trimws(main_lineages_list)))
    writeLines(main_lineages_list, "selected_main-lineage_split.list")

    cat("Number of sample IDs:", length(sample_ids), "\\n")
    cat("Number of sub-lineages:", length(sub_lineages_list), "\\n")
    cat("Number of main lineages:", length(main_lineages_list), "\\n")

    ########################################################################################################
    # IMPORT AND PREPARE DATA
    ########################################################################################################

    # Import dataframes
    main_lineages_df <- data.frame(selected_main_lineage = main_lineages_list)
    sub_lineages_df <- data.frame(selected_sub_lineage = sub_lineages_list)
    run_ids_df <- data.frame(SampleID = sample_ids)

    # Read the pairwise analysis file
    meta <- readr::read_delim(pairwise_analysis_file,
                            col_names = TRUE, 
                            show_col_types = FALSE) |>
        distinct()
    colnames(meta) <- c("SampleID", "main_lineage", "sub_lineage")

    cat("Original meta data rows:", nrow(meta), "\\n")

    # Filter out any sample that contains 'CN-' and clean data
    meta <- meta |>
        filter(!is.na(main_lineage) & !str_detect(main_lineage, ";")) |>
        filter(!str_detect(SampleID, "CN-"))

    cat("Filtered meta data rows:", nrow(meta), "\\n")

    ########################################################################################################
    # ANALYSIS TYPE SPECIFIC PROCESSING
    ########################################################################################################

    if (pairwise_split == "sub") {

        cat("Processing sub-lineage analysis\\n")

        # Filter out the lineages at sub_lineage level
        filtered_meta <- meta |>
            mutate(
                # Check if sub_lineage belongs to any selected_sub_lineage
                matched_sub = sapply(sub_lineage, function(sub) {
                    if (is.na(sub) || sub == "") return(NA)
                    match <- sub_lineages_df\$selected_sub_lineage[
                        str_detect(sub, paste0("^", sub_lineages_df\$selected_sub_lineage, "\\\\b"))
                    ]
                    if (length(match) > 0) match[1] else NA
                }),
                # Determine final lineage
                final_lineage = case_when(
                    !is.na(matched_sub) ~ matched_sub,  # If sub_lineage matches, use it
                    main_lineage %in% main_lineages_df\$selected_main_lineage ~ main_lineage,  # Otherwise, use main_lineage if valid
                    TRUE ~ sub_lineage  # Default to sub_lineage if no match
                )
            ) |>
            select(final_lineage, SampleID) |>
            rename(lineage = final_lineage)

    } else if (pairwise_split == "main") {

        cat("Processing main lineage analysis\\n")

        # Select only the main lineages
        filtered_meta <- meta |>
            select(lineage = main_lineage, SampleID)

    } else if (pairwise_split == "none") {

        cat("Processing all samples without lineage split\\n")

        # For 'none' type, assign all samples to a single group
        filtered_meta <- meta |>
            mutate(lineage = "All") |>
            select(lineage, SampleID)

    } else {

        stop("Invalid pairwise level specified: ", pairwise_split, 
            ". Choose 'sub', 'main', or 'none', then re-run the workflow with the correct parameter and '-resume'")

    }

    ########################################################################################################
    # COMMON PROCESSING AND FILTERING
    ########################################################################################################

    cat("Processed filtered_meta rows:", nrow(filtered_meta), "\\n")

    # Find lineages with too few samples (< 5)
    skipped_lineages_to_few <- filtered_meta |>
        group_by(lineage) |>
        count() |>
        filter(n < 5)

    cat("Lineages with fewer than 5 samples (will be skipped):\\n")
    if (nrow(skipped_lineages_to_few) > 0) {
        print(skipped_lineages_to_few)
    } else {
        cat("None\\n")
    }

    # Get unique lineages from run_ids that pass the minimum sample threshold
    filtered_lineages_pass <- filtered_meta |>
        filter(SampleID %in% run_ids_df\$SampleID) |>
        filter(!lineage %in% skipped_lineages_to_few\$lineage) |>
        count(lineage) |>
        select(lineage) |>
        distinct()

    cat("Lineages that will be processed:\\n")
    if (nrow(filtered_lineages_pass) > 0) {
        print(filtered_lineages_pass)
    } else {
        cat("None\\n")
    }

    # Partition filtered_meta based on filtered_lineages_pass
    filtered_meta_forward <- filtered_meta |>
        filter(lineage %in% filtered_lineages_pass\$lineage)  |>
        filter(SampleID != "sample") |>
        filter(SampleID != "SampleID") |>
        filter(lineage != "main_lineage")  |>
        filter(lineage != "lineage")

    filtered_meta_skip <- filtered_meta |>
        filter(!lineage %in% filtered_lineages_pass\$lineage) |>
        filter(SampleID != "sample") |>
        filter(SampleID != "SampleID") |>
        filter(lineage != "sub_lineage")

    cat("Samples to be processed:", nrow(filtered_meta_forward), "\\n")
    cat("Samples to be skipped:", nrow(filtered_meta_skip), "\\n")

    ########################################################################################################
    # OUTPUT FILES
    ########################################################################################################

    # Export the csv files WITHOUT headers using write.table instead of write.csv
    write.table(filtered_meta_forward, "final.lineage_samples_tuple.csv",
                sep = ",",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)

    write.table(filtered_meta_skip, "final.skipped-lineages_tuple.csv",
                sep = ",",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE)

    cat("Analysis completed successfully!\\n")
    cat("================================\\n")
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 3.0.0
@description:
    In this module creates the pairwise analysis tuples from the lineage_samples_paths.csv
    and the lineage_pairwise_sub and lineage_pairwise_main lists.
    The output is a tuple of the form (lineage, sampleID) for each sampleID in the analysis.
    There are three options for the pairwise analysis (specified by the params.pairwise_split):
        - sub: pairwise analysis at sub-lineage level
        - main: pairwise analysis at main-lineage level
        - none: pairwise analysis of all samples without lineage split
    Now implemented entirely in R for better data handling and consistency.
@changelog:
    v3.0.0-2026-04-13: Complete rewrite in R, eliminating bash script components
    v2.0.0-2025-04-01: Updated to use the new lineage_pairwise_sub and lineage_pairwise_main lists
    v1.0.1-2024-11-01: Added error handling for invalid pairwise level
*/