process COMPILE_SEQUENCING_STATS {

    conda params.r_stats_env

    publishDir "${params.outDir}/db/comparison/summary/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple path(tbdb_results, stageAs: "tbdb-tbprofiler.txt"), 
            path(who_results, stageAs: "who-tbprofiler.txt"),
            path(lineage_fractions, stageAs: "lineages.fractions.txt")

        path(mtbseq_strains, stageAs: "Strain_Classification.tab")
        path(mtbseq_stats, stageAs: "Mapping_and_Variant_Statistics.tab")

    output:
        path("${params.runID}.sequencing_summary.csv")

        path("sequencing_summary.csv"),      emit: analysis_summary
        path("who_resistance_summary.csv"),  emit: who_resistance
        path("tbdb_resistance_summary.csv"), emit: tbdb_resistance

        path("pairwise_analysis.list.csv"),  emit: pairwise_analysis_list

    script:

    """
    #!/usr/bin/env Rscript

    # Set a default CRAN repo to avoid mirror prompts
    options(repos = c(CRAN = "https://cloud.r-project.org"))

    # load libraries
    packages <- c("tidyverse", "dplyr", "tidyr")

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
    # Params
    #··············································································#

    dictionary  <- "${params.scriptDir}/R/dict/mtbseq_classification.dict.csv"
    runID       <- "${params.runID}"
    
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
    # Create Lineage fractions file
    #··············································································#

    # Read the input files
    tbprof <- read.delim("tbdb-tbprofiler.txt", header = TRUE) |>
        mutate(across(everything(), ~ str_trim(as.character(.)))) |>
        select(SampleID=sample,main_lineage,sub_lineage) |> # designate 'clonal' or 'mixed
        mutate(infection_type = ifelse(grepl(';', main_lineage) | grepl(';', sub_lineage), 'mixed', 'clonal'))

    tbprof.infec <- tbprof |> select(SampleID,infection_type)

    lineage.frac <- read.delim("lineages.fractions.txt", header = TRUE) |>
        filter(Fraction != 'Fraction')

    # Split 'clonal' and 'clonal'
    # clonal: add fraction of lineage designation
    mixed <- tbprof |> filter(infection_type == 'mixed')
        if (nrow(mixed) == 0) { 
            # Check if 'mixed' is empty
            mixed.frac.5 <- mixed

            } else { # Proceed with processing if 'mixed' is not empty

            mixed.frac <- left_join(mixed, lineage.frac) |>
                select(SampleID, main_lineage, sub_lineage, infection_type, Lineage, Fraction)
            mixed.frac <- mixed.frac |> mutate(across(Fraction, as.numeric))
            mixed.frac <- mixed.frac |> mutate(Perc = 100 * Fraction)

            mixed.frac.1 <- mixed.frac |>
                mutate(Lineage_p = paste0(Lineage, " (", Perc, "%)")) |>
                mutate(Lineage_p = gsub("lineage", "L", Lineage_p))

            mixed.frac.2 <- mixed.frac.1 |>
                select(SampleID, sub_lineage, Lineage_p) |>
                mutate(sub_lineage = gsub("lineage", "L", sub_lineage)) |>
                mutate(Lineage = Lineage_p) |>
                separate_wider_delim(cols = Lineage, " ", names = c("Lineage", NA))

            mixed.frac.2.keep <- mixed.frac.2 |>
                select(SampleID, sub_lineage) |>
                distinct(SampleID, sub_lineage, .keep_all = TRUE) |>
                separate_rows(sub_lineage, sep = ";") |>
                select(SampleID, Lineage = sub_lineage)

            mixed.frac.3 <- left_join(mixed.frac.2.keep, mixed.frac.2) |>
                                group_by(SampleID) |>
                                mutate(Lineage_pL = paste(unique(Lineage_p), collapse = "; ")) |>
                                slice(1) |>  # Keeps only the first row per SampleID after mutating
                                ungroup() |>
                                select(SampleID, Lineage_frac = Lineage_pL)

            mixed.frac.4 <- mixed.frac.1 |> 
                                group_by(SampleID) |>  # Group by SampleID
                                arrange(desc(Perc), .by_group = TRUE) |>   # Sort by Perc in descending order within each group
                                slice(1) |> # Retain only the first row of each group
                                filter(Perc >= 90) |> # Filter rows where Perc >= 90
                                mutate(Mixed_90perc = paste0(Lineage, " (", Perc, "%)")) |>  # Create new column
                                ungroup() # Ungroup after manipulation (optional but recommended)

            mixed.frac.5 <- left_join(mixed.frac.3, mixed.frac.4) |>
                            select(SampleID, Lineage_frac, Mixed_90perc)
    }

    # Clonal: add fraction of lineage designation
    clonal <- tbprof |> filter(infection_type == "clonal")

    if (nrow(clonal) == 0) { # Check if 'clonal' is empty
            clonal.frac.3 <- clonal
        } else { 
        # Proceed with processing if 'clonal' is not empty
        clonal.frac <- left_join(clonal, lineage.frac) |>
            select(SampleID, main_lineage, sub_lineage, 
                    infection_type, Lineage, Fraction)
        clonal.frac <- clonal.frac |> mutate(across(Fraction, as.numeric))
        clonal.frac <- clonal.frac |> mutate(Perc = 100 * Fraction)

        clonal.frac.1 <- clonal.frac |>
                                    mutate(Lineage_p = paste0(Lineage, " (", Perc, "%)")) |>
                                    mutate(Lineage_p = gsub("lineage", "L", Lineage_p))

        clonal.frac.2 <- clonal.frac.1 |>
                                    select(SampleID, sub_lineage, Lineage_p) |>
                                    mutate(sub_lineage = gsub("lineage", "L", sub_lineage)) |>
                                    mutate(Lineage = Lineage_p) |>
                                    separate_wider_delim(cols = Lineage, " ", names = c("Lineage", NA))

        clonal.frac.2.keep <- clonal.frac.2 |>
            select(SampleID, sub_lineage) |>
            distinct(SampleID, sub_lineage, .keep_all = TRUE) |>
            separate_rows(sub_lineage, sep = ";") |>
            select(SampleID, Lineage = sub_lineage)

        clonal.frac.3 <- left_join(clonal.frac.2.keep, clonal.frac.2) |>
            group_by(SampleID) |>
            mutate(Lineage_pL = paste(unique(Lineage_p), collapse = "; ")) |>
            slice(1) |>  # Keeps only the first row per SampleID after mutating
            ungroup() |>
            select(SampleID, Lineage_frac = Lineage_pL) |>
            distinct() |>
            mutate(Mixed_90perc = NA)
    }

    # Combine clonal and clonal results
    tbdb_lin_fract <- rbind(clonal.frac.3, 
        mixed.frac.5) |>
        arrange(desc(SampleID))

    tbdb_lin_fract_final <- left_join(tbdb_lin_fract, tbprof.infec) |>
        select(SampleID,infection_type,
        Lineage_frac,Mixed_90perc)

    # Output the result
    write.table(tbdb_lin_fract_final, 
        file = "tbprofiler.lineages.fractions.txt",
        sep = ";", 
        row.names = FALSE, 
        col.names = TRUE, 
        quote = FALSE)

    #··············································································#
    # Sequencing Summary
    #··············································································#

    # Import dataframes
    mtbseq_stats    <- read.delim("Mapping_and_Variant_Statistics.tab", header = TRUE)
    mtbseq_class    <- read.delim("Strain_Classification.tab", header = TRUE)
    tbprofiler_tbdb <- read.delim("tbdb-tbprofiler.txt", header = TRUE)
    tbprofiler_who  <- read.delim("who-tbprofiler.txt", header = TRUE)
    lineage_frac    <- read.delim("tbprofiler.lineages.fractions.txt", sep = ";", header = TRUE) |> 
                                    select(SampleID, Lineage_frac, Mixed_90perc)

    # Merge the dataframes
    merge1 <- left_join(mtbseq_statistics.df, mtbseq_classification.df)
    merge2 <- left_join(merge1, tbprofiler_tbdb, by = c("FullID" = "sample"))
    merge3 <- left_join(merge2, lineage_frac, by = c("FullID" = "SampleID"))
    full.df <- merge3 # rename the dataframe

    # Filter out the genomes using a minimum coverage 
    full.df.final <- full.df |>
        mutate(infection_type = if_else(
            grepl(";", main_lineage) | grepl(";", sub_lineage),
                "Mixed","Clonal")
        )

    # Creating the outputs for later assembly into final results
    sequencing_summary.df <- dictionary_rename(
        df = full.df.final,
        dict_path = "${params.scriptDir}/R/dict/sequencing-summary.dict.csv")

    resistance_profiles_TBDB.df <- dictionary_rename(
        df = tbprofiler_tbdb,
        dict_path = "${params.scriptDir}/R/dict/resistance_profiles_TBDB.csv")

    resistance_profiles_WHO.df <- dictionary_rename(
        df = tbprofiler_who,
        dict_path = "${params.scriptDir}/R/dict/resistance_profiles_WHO.csv")

    # Create the list of genomes for pairwise analysis
    pairwise_analysis.df <- full.df.final |>
            select(SampleID=FullID,
                    main_lineage,
                    sub_lineage) |>
            filter(main_lineage != "NA" & !str_detect(main_lineage, ";") & !str_detect(sub_lineage, ";"))

    # export all the ouputs (broken!)
    write.csv(sequencing_summary.df, 
        "sequencing_summary.csv",
        quote = TRUE, 
        row.names = FALSE
        )

    write.csv(sequencing_summary.df,
        "${params.runID}.sequencing_summary.csv",
        quote = TRUE, 
        row.names = FALSE
        )

    write.csv(resistance_profiles_TBDB.df,
        "tbdb_resistance_summary.csv",
        quote = TRUE, 
        row.names = FALSE
        )

    write.csv(resistance_profiles_WHO.df,
        "who_resistance_summary.csv",
        quote = TRUE, 
        row.names = FALSE
        )

    write.csv2(pairwise_analysis.df,
        "pairwise_analysis.list.csv",
        quote = TRUE, 
        row.names = FALSE
        )

    #··············································································#
    # Create split tuple
    #··············································································#

    # Import the samplesheet
    samplelist <- read.csv("${params.samplesheet}", 
            header = TRUE) |> 
        select(sampleID) |>
        deframe()

    # Use the MTBseq mapping statistics for identifying good genomes
    minQual_genomes_list <- read.delim("Mapping_and_Variant_Statistics.tab", header = FALSE) |>
        distinct() |> 
        filter(V5 >= ${params.filt_min_reads} & 
                V16 >= ${params.filt_min_cov} & 
                V19 >= ${params.filt_min_depth}) |>
        select(V4) |> 
        distinct() |> deframe()

    # get list of clonal genomes
    clonal_genomes_list <- tbdb_lin_fract_final |> 
        filter(infection_type == "clonal") |>
        select(SampleID) |> deframe()
    
    # Get a list of the genomes with their lineages
    minQual_genomes <- tbprofiler_tbdb |>
        filter(sample %in% minQual_genomes_list) |>
        select(sample, main_lineage, sub_lineage) |>
        filter(sample %in% clonal_genomes_list)
    
    # Export the pairwise analysis list
    write.csv(minQual_genomes, 
        "pairwise_analysis.list.csv", 
        quote = FALSE, 
        row.names = FALSE,
        col.names = FALSE
        )
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.1
@description:
    In this module the sequencing statistics for the db are
    calculated with Rscripts. The output is a summary of the
    sequencing statistics, the tbdb and who resistance  summaries
    and a list of the genomes which pass the minimum quality
    requirements for the pairwise analysis. This is then used to create a 
    filtered list of genomes for the pairwise analysis, which is converated into
    a tuple/channel for downstream processing.
@changelog:
    v1.0.0_2024-04-01: Inital version of module
    v1.0.1-2026-05-13: Removed renaming the MTBSeq concatenated files as they now include headers
*/
