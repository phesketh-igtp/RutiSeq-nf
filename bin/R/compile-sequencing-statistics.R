# load libraries
    library(argparse)
    library(tidyverse)

# Create a parser for the script
    parser <- ArgumentParser(description = "Script to process MTBseq and TBProfiler data")

# Define arguments
    parser$add_argument("--mtbseq_statistics", required=TRUE, help="Path to MTBseq compiled map statistics file")
    parser$add_argument("--mtbseq_classification", required=TRUE, help="Path to MTBseq compiled strains file")
    parser$add_argument("--tbprofiler_tbdb", required=TRUE, help="Path to TBProfiler TBDB results file")
    parser$add_argument("--tbprofiler_who", required=TRUE, help="Path to TBProfiler TBDB results file")
    parser$add_argument("--minimum_coverage", required=TRUE, type="integer", help="Minimum coverage threshold")
    parser$add_argument("--dictionary_path", default=NULL, help="Path to R dictioanry for renaming files")
    ## parser$add_argument("--additional_args", default=NULL, help="Optional additional arguments from config file")
    parser$add_argument("--output", required=TRUE, help="Output file name for sequencing summary CSV")

# Parse the arguments
    args <- parser$parse_args()

# Assign arguments to variables
    mtbseq_statistics <- args$mtbseq_statistics
    mtbseq_classification <- args$mtbseq_classification
    tbprofiler_tbdb <- args$tbprofiler_tbdb
    output_file <- args$output
    minimum_coverage <- args$minimum_coverage
    dictionary_path <- args$dictionary_path
## additional_args <- args$additional_args

# Debugging: Print the arguments (optional)
    cat("Arguments received:\n")
    cat("MTBseq Statistics: ", mtbseq_statistics, "\n")
    cat("MTBseq Classification: ", mtbseq_classification, "\n")
    cat("TBProfiler TBDB: ", tbprofiler_tbdb, "\n")
    cat("Output File: ", output_file, "\n")
    cat("Minimum Coverage: ", minimum_coverage, "\n")

# Import dataframes
    mtbseq_stats        <- read.delim(mtbseq_statistics, header = FALSE); cat("Loading MTBseq statistics...\n")
    mtbseq_class        <- read.delim(mtbseq_classification, header = FALSE); cat("Loading MTBseq classifications...\n")
    tbprofiler_tbdb.df  <- read.delim(tbprofiler_tbdb, header = TRUE); cat("Loading TB-profiler TBDB results...\n")
    tbprofiler_who      <- read.delim(tbprofiler_who, header = TRUE); cat("Loading TB-profiler WHO results...\n")

# Add appropriate headers to dataframes using dictionaries

    # import local functions
        source(paste0(dictionary_path,"/functions/dictionary-rename.R"))

    mtbseq_statistics.df <- dictionary_rename(df = mtbseq_stats,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_statistics.dict.csv"))

    mtbseq_classification.df  <- dictionary_rename(df = mtbseq_class,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_classification.dict.csv"))

    merge1 <- left_join(mtbseq_statistics.df, mtbseq_classification.df)
    merge2 <- left_join(merge1, tbprofiler_tbdb.df, by = c("FullID" = "sample"))
    full.df <- merge2 # rename the dataframe

# Filter out the genomes using a minimum coverage 
    full.df.final <- full.df |>
                mutate(infection_type = if_else(
                    grepl(";", main_lineage) | grepl(";", sub_lineage),
                        "Mixed","Clonal")
                    )

# Creating the outputs for later assembly into final results
    sequencing_summary.df <- dictionary_rename(df = full.df.final,
                                dict_path = paste0(dictionary_path,
                                "/dict/sequencing-summary.dict.csv"))

    resistance_profiles_TBDB.df <- dictionary_rename(df = tbprofiler_tbdb.df,
                                dict_path = paste0(dictionary_path,
                                "/dict/resistance_profiles_TBDB.csv"))

# Create the list of genomes for pairwise analysis

    pairwise_analysis.df <- full.df.final |>
                filter(`(Unambiguous) Coverage median` >= minimum_coverage) |> 
                filter(infection_type == "Clonal") |> #only clonal genomes
                filter(sub_lineage != "NA") |> # no unclassified genomes
                select(SampleID=FullID,sub_lineage)

    pairwise_analysis.df

# export all the ouputs (broken!)
write.csv2(sequencing_summary.df,
            "sequencing_summary.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(resistance_profiles_TBDB.df,
            "resistance_profiles_TBDB.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(pairwise_analysis.df,
            "pairwise_analysis.list.csv",
            quote=FALSE, row.names=FALSE)
