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

    publishDir "${params.outDir}/results/${params.runID}/", 
        mode: 'copy', 
        overwrite: true

    input:
        path pairwise_clusters_processed

        path sequencing_summary

        path who_resistance
        path tbdb_resistance

        tuple path(sylph_sequence_abundance, stageAs: 'sylph_results.tsv'),
            path(sylph_relative_abundance, stageAs: 'sylph_relative_abundance.tsv'),
            path(sylph_coverage, stageAs: 'sylph_coverage.tsv')

        tuple path(html_report), 
            path(warnings, stageAs: 'warnings.out')
        

    output:
        path("${params.runID}_RutiSeq-results.xlsx")
        path(html_report)

    script:

        """
        #!/usr/bin/env Rscript

        # load libraries
        packages <- c(
            "dplyr", "tidyverse", "openxlsx")

        # Identify missing packages
        missing_pkgs <- packages[!packages %in% installed.packages()[, "Package"]]

        # Install missing packages
        if (length(missing_pkgs) > 0) {
            install.packages(missing_pkgs, dependencies = TRUE)
        }

        # Load all packages
        invisible(lapply(packages, library, character.only = TRUE))

        set.seed(1234)

        #··············································································#
        # Functions
        #··············································································#
        
        dictionary_rename <- function(df, dict_path) { 
            # import dictionary
                dict <- read.csv(dict_path)
            # create vector name
                dict_names <- dict |> 
                    select(new.name,old.name) |>
                    deframe()
            # Create vector for cols to keep (using 'new.name' since the cols will be renamed)
                cols.to.keep <- dict |> 
                        filter(final == "Y") |>
                            select(new.name) |>
                                deframe()
            # Rename cols
                df <- df |> 
                    rename(all_of(dict_names)) |> # Rename cols using dict_names
                        select(all_of(cols.to.keep)) # Keep only cols.to.keep
        }

        #··············································································#
        # Import all dataframes
        #··············································································#

        summary <- read.delim("sequencing_summary.csv",  header = TRUE,
                sep = ",", 
                check.names = FALSE
                ) |> 
            distinct()

        who_res <- read.delim("who_resistance_summary.csv",  header = TRUE,
                sep = ",",
                check.names = FALSE
                ) |> 
            distinct()

        tbdb_res <- read.delim("tbdb_resistance_summary.csv", header = TRUE,
                sep = ",", 
                check.names = FALSE
                ) |> 
            distinct()

        clusters <- read.delim("processed_clusters.tsv", header = TRUE,
                sep = "\t",
                check.names = FALSE
                ) |> 
            distinct()

        sylph_out <- read.delim("sylph_results.tsv", 
                header = TRUE,
                sep = "\t", 
                check.names = FALSE
                )

        #··············································································#
        # Create the final results summary of the dataframe
        #··············································································#

        who_res.min <- who_res |> 
            select(Sample=sample, `DRtype (WHO)` = `DR type`,
            INH, EMB, PZA, MFX, LFX, BDQ, DLM, Pa, LZD, STM, 
            AMK, KAN, CAP, CFZ, ETO, PAC, CYR)
        tbdb_res.min <- tbdb_res |> 
            select(Sample=sample, `DRtype (TBDB)` = `DR type`)

        resistance_summary <- left_join(who_res.min, tbdb_res.min, by = "Sample")

        summary <- summary |> 
            mutate(OriginalID = sub("_.*", "", Sample))

        summary.tmp <- left_join(summary, resistance_summary, by = "Sample") |> 
            distinct()

        summary_xlsx <- dictionary_rename(
            df = summary.tmp,
            dict_path = "${params.scriptDir}/R/dict/xlsx_sheet1.dict.csv"
            )

        #··············································································#
        # Assemble the XLSX spreadsheet
        #··············································································#

        # Create a new workbook
        wb <- createWorkbook()

        # Add each sheet to the workbook
        addWorksheet(wb, "sequencing_summary")
            writeData(wb, 
                "sequencing_summary", 
                summary_xlsx, 
                rowNames = FALSE
                )

        addWorksheet(wb, "resistance_table_who")
            writeData(wb, "resistance_table_who", 
                who_res,
                rowNames = FALSE
                )

        addWorksheet(wb, "resistance_table_tbdb")
            writeData(wb,
                "resistance_table_tbdb",
                tbdb_res, 
                rowNames = FALSE
                )

        addWorksheet(wb, "transmission_clusters")
            writeData(wb,
                "transmission_clusters",
                clusters,
                rowNames = FALSE)

        addWorksheet(wb, "sylph_classification")
            writeData(wb,
                "sylph_classification",
                sylph_out,
                rowNames = FALSE)

        # Save the workbook
        saveWorkbook(wb, 
                file = "${params.runID}_RutiSeq-results.xlsx", 
                overwrite = TRUE
            )
        """
}