process CONCATENATE_CLUSTERS {

    conda params.r_stats_env

    publishDir "${params.outDir}/db/results/main/", mode: 'copy'

    input:
        path(clusters)
        path(sequencing_summary, stageAs: "sequencing_summary.csv")

    output:
        path("unprocessed_clusters.tsv"),   emit: pairwise_clusters_unprocessed
        path("processed_clusters.tsv"),     emit: pairwise_clusters_processed

    script:

        """
        #!/usr/bin/env Rscript --vanilla

        # Enable comprehensive error reporting
        options(error = function() {
            cat("ERROR occurred at:", date(), "\\n")
            cat("Traceback:\\n")
            traceback()
            #quit(status = 1)
        })

        cat("=== CONCATENATE_CLUSTERS DEBUG START ===\\n")
        cat("R version:", R.version.string, "\\n")
        cat("Working directory:", getwd(), "\\n")
        cat("Files in working directory:\\n")
        print(list.files(recursive = TRUE))

        #--------------------------------------------------------------------------#
        # Load packages
        #--------------------------------------------------------------------------#

        # Set a default CRAN mirror (avoids prompt)
        options(repos = c(CRAN = "https://cloud.r-project.org"))

        # List of packages
        packages <- c("tidyverse", "dplyr")

        # Check installed packages
        cat("Installed packages:\\n")
        installed_pkgs <- installed.packages()[, "Package"]
        cat("tidyverse installed:", "tidyverse" %in% installed_pkgs, "\\n")
        cat("dplyr installed:", "dplyr" %in% installed_pkgs, "\\n")

        # Identify missing packages
        missing_pkgs <- packages[!packages %in% installed_pkgs]
        cat("Missing packages:", paste(missing_pkgs, collapse = ", "), "\\n")

        # Install missing packages (non-interactive)
        if (length(missing_pkgs) > 0) {
            cat("Installing missing packages...\\n")
            tryCatch({
                install.packages(missing_pkgs, dependencies = TRUE)
                cat("Package installation completed\\n")
            }, error = function(e) {
                cat("ERROR installing packages:", e\$message, "\\n")
                #quit(status = 1)
            })
        }

        # Load all packages
        cat("Loading packages...\\n")
        for (pkg in packages) {
            tryCatch({
                library(pkg, character.only = TRUE)
                cat("Successfully loaded:", pkg, "\\n")
            }, error = function(e) {
                cat("ERROR loading package", pkg, ":", e\$message, "\\n")
                #quit(status = 1)
            })
        }

        #--------------------------------------------------------------------------#
        # Check input files and directories
        #--------------------------------------------------------------------------#
        cat("\\n=== CHECKING INPUTS ===\\n")

        # Check sequencing summary file
        if (!file.exists("sequencing_summary.csv")) {
            cat("ERROR: sequencing_summary.csv not found\\n")
            #quit(status = 1)
        }
        cat("sequencing_summary.csv found\\n")

        # Define base directory
        base_dir <- "${params.outDir}/db/comparison/mtbseq/"
        cat("Base directory:", base_dir, "\\n")
        
        if (!dir.exists(base_dir)) {
            cat("ERROR: Base directory does not exist:", base_dir, "\\n")
            cat("Available directories:\\n")
            print(list.dirs("${params.outDir}", recursive = TRUE))
            #quit(status = 1)
        }
        cat("Base directory exists\\n")

        # List all cluster and singleton files
        cat("\\n=== FINDING FILES ===\\n")
        cluster_files <- list.files(path = file.path(base_dir), 
                                    pattern = "_d.*\\\\.processed\\\\.clusters\\\\.tsv\$", 
                                    recursive = TRUE, 
                                    full.names = TRUE)
        
        singleton_files <- list.files(path = file.path(base_dir), 
                                    pattern = "_d.*\\\\.singletons\\\\.tsv\$", 
                                    recursive = TRUE, 
                                    full.names = TRUE)

        cat("Found", length(cluster_files), "cluster files\\n")
        cat("Found", length(singleton_files), "singleton files\\n")

        # Combine all files into a single vector
        all_files <- c(cluster_files, singleton_files)
        cat("Total files to process:", length(all_files), "\\n")

        if (length(all_files) == 0) {
            cat("ERROR: No cluster or singleton files found\\n")
            cat("Files in base directory:\\n")
            print(list.files(base_dir, recursive = TRUE))
            #quit(status = 1)
        }

        cat("Files to process:\\n")
        print(all_files)

        #--------------------------------------------------------------------------#
        # Import all the dataframes and merge
        #--------------------------------------------------------------------------#
        cat("\\n=== READING FILES ===\\n")

        # Read all files (no header)
        df_list <- list()
        for (i in seq_along(all_files)) {
            f <- all_files[i]
            cat("Reading file", i, "of", length(all_files), ":", basename(f), "\\n")
            tryCatch({
                df <- read.table(f, sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
                cat("  - Rows:", nrow(df), "Cols:", ncol(df), "\\n")
                df_list[[i]] <- df
            }, error = function(e) {
                cat("ERROR reading file", f, ":", e\$message, "\\n")
                #quit(status = 1)
            })
        }

        cat("Merging dataframes...\\n")
        # Merge all into one data frame
        merged_df <- do.call(rbind, df_list) |>
            distinct() |>
            tidyr::separate_wider_delim(
                V4,
                delim = "-",
                names = c("num", "dist", "lin"),
                too_few = "align_start"
            ) |>
            dplyr::mutate(
                num  = stringr::str_pad(num,  width = 3, pad = "0"),
                dist = stringr::str_pad(dist, width = 2, pad = "0")
            ) |> 
            mutate(V4 = paste0(num, "-", 
                                dist, "-", 
                                lin,
                                sep = "")) |>
            select(V1, V2, V3, V4)

        cat("Merged dataframe: Rows:", nrow(merged_df), "Cols:", ncol(merged_df), "\\n")

        #--------------------------------------------------------------------------#
        # Wrangle the files
        #--------------------------------------------------------------------------#
        cat("\\n=== PROCESSING DATA ===\\n")

        raw_clust <- merged_df
        colnames(raw_clust) <- c("lineage", "dSNP", "genomes", "int_clusterID")

        raw_clust <- raw_clust |>
            mutate(dSNP = gsub("dist_", "t=", dSNP))

        cat("Data after initial processing: Rows:", nrow(raw_clust), "\\n")

        # Process the sampleIDs
        cat("Reading sequencing summary...\\n")
        tryCatch({
            sample_id_df <- read.delim("sequencing_summary.csv", header = TRUE, sep = ",") |>
                select(Sample) |>
                mutate(genomes = Sample) |>
                separate_wider_delim(genomes, delim = "_",
                            names = c("id", "library"),
                            too_few = "align_end",
                            too_many = "drop") |>
                distinct() |>
                select(SampleID = Sample, genomes = id)
            cat("Sample ID dataframe: Rows:", nrow(sample_id_df), "\\n")
        }, error = function(e) {
            cat("ERROR processing sequencing summary:", e\$message, "\\n")
            #quit(status = 1)
        })

        # Merge and pivot the data
        cat("Merging with sample IDs...\\n")
        raw_clust <- left_join(raw_clust, sample_id_df) |>
            distinct() |>
            select(SampleID, lineage, dSNP, int_clusterID) |>
            pivot_wider(names_from = dSNP, values_from = int_clusterID)

        cat("Data after pivot: Rows:", nrow(raw_clust), "Cols:", ncol(raw_clust), "\\n")

        # Extract and sort t= column names numerically
        t_cols <- names(raw_clust)[grepl("^t=", names(raw_clust))]
        t_cols_sorted <- t_cols[order(as.numeric(gsub("t=", "", t_cols)))]

        cat("Found t= columns:", paste(t_cols_sorted, collapse = ", "), "\\n")

        # Reorder the full data frame
        raw_clust <- raw_clust |>
            select(SampleID, lineage, all_of(t_cols_sorted))

        # Create the merged clusterID
        raw_clust <- raw_clust |>
            unite("merged_clusterID", all_of(t_cols_sorted), sep = "/", remove = FALSE) |>
            select(SampleID, lineage, all_of(t_cols_sorted), merged_clusterID)

        cat("Final processed data: Rows:", nrow(raw_clust), "Cols:", ncol(raw_clust), "\\n")

        #--------------------------------------------------------------------------#
        # Export the dataframe
        #--------------------------------------------------------------------------#
        cat("\\n=== WRITING OUTPUT ===\\n")

        tryCatch({
            write.table(raw_clust, "processed_clusters.tsv",
                        sep = "\\t", 
                        row.names = FALSE, 
                        quote = FALSE)
            cat("Successfully wrote processed_clusters.tsv\\n")
        }, error = function(e) {
            cat("ERROR writing processed_clusters.tsv:", e\$message, "\\n")
            #quit(status = 1)
        })

        tryCatch({
            write.table(merged_df, file = "unprocessed_clusters.tsv", 
                        sep = "\\t", 
                        row.names = FALSE, 
                        col.names = FALSE, 
                        quote = FALSE)
            cat("Successfully wrote unprocessed_clusters.tsv\\n")
        }, error = function(e) {
            cat("ERROR writing unprocessed_clusters.tsv:", e\$message, "\\n")
            #quit(status = 1)
        })
        """
}

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 2.0.0
    @description:
        This process concatenates all the cluster files generated by the pairwise analysis
        into a single file. The output is a tab-separated file with the following columns:
            - lineage
            - distance
            - genomes
            - group
        The output file is named "unprocessed_clusters.tsv".
    @changelog
        v1.0.0-2025-04-01: Initial version
        v1.0.1-2025-04-04: Added comprehensive error handling and debugging information
        v2.0.0-2025-12-03: Improved package installation and loading with error handling, 
                            added detailed input checks and logging
                        Moved from running a seperate Rscript to embedding the Rcode directly
                            into the nextflow module.
*/
