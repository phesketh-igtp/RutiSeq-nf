# load libraries
library(argparse)
library(tidyverse)

# Create a parser for the script
    parser <- ArgumentParser(description = "Script to process MTBseq and TBProfiler data")

# Define arguments
parser$add_argument("--minimum_coverage", required=TRUE, type="integer", help="Minimum coverage threshold")
parser$add_argument("--dictionary_path", default=NULL, help="Path to R dictioanry for renaming files")
parser$add_argument("--runID", required=TRUE, help="RunID")

## parser$add_argument("--additional_args", default=NULL, help="Optional additional arguments from config file")

# Parse the arguments
args <- parser$parse_args()

# Assign arguments to variables
mtbseq_statistics     <- args$mtbseq_statistics
mtbseq_classification <- args$mtbseq_classification
tbprofiler_tbdb       <- args$tbprofiler_tbdb
tbprofiler_who        <- args$tbprofiler_who
minimum_coverage      <- args$minimum_coverage
dictionary_path       <- args$dictionary_path
lineage_fractions     <- args$lineage_fractions
runID                 <- args$runID
## additional_args <- args$additional_args

# Import dataframes
mtbseq_stats       <- read.delim("Mapping_and_Variant_Statistics.tab", header = FALSE)
mtbseq_class       <- read.delim("Strain_Classification.tab", header = FALSE)
tbprofiler_tbdb.df <- read.delim("tbdb-tbprofiler.txt", header = TRUE)
tbprofiler_who     <- read.delim("who-tbprofiler.txt", header = TRUE)
lineage_frac       <- read.delim("tbprofiler.lineages.fractions.txt", header = TRUE) |> 
                                select(SampleID, Lineage_frac, Mixed_90perc)

# Add appropriate headers to dataframes using dictionaries

# import local functions ## dictionary_path = "/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R/"
source(paste0(dictionary_path,"/functions/dictionary-rename.R"))

mtbseq_statistics.df     <- dictionary_rename(df = mtbseq_stats,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_statistics.dict.csv"))

mtbseq_classification.df <- dictionary_rename(df = mtbseq_class,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_classification.dict.csv"))

merge1 <- left_join(mtbseq_statistics.df, mtbseq_classification.df)
merge2 <- left_join(merge1, tbprofiler_tbdb.df, by = c("FullID" = "sample"))
merge3 <- left_join(merge2, lineage_frac, by = c("FullID" = "SampleID"))
full.df <- merge3 # rename the dataframe

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

resistance_profiles_WHO.df <- dictionary_rename(df = tbprofiler_who,
                        dict_path = paste0(dictionary_path,
                        "/dict/resistance_profiles_WHO.csv"))

# Create the list of genomes for pairwise analysis

pairwise_analysis.df <- full.df.final |>
                filter(`Unambiguous Coverage median` >= minimum_coverage) |>
                filter(infection_type == "Clonal") |> #only clonal genomes
                filter(sub_lineage != "NA") |> # no unclassified genomes
                select(SampleID=FullID,sub_lineage)

# export all the ouputs (broken!)
write.csv2(sequencing_summary.df,
            "sequencing_summary.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(sequencing_summary.df,
            paste0(runID, ".sequencing_summary.csv", sep = ""),
            quote=FALSE, row.names=FALSE)

write.csv2(resistance_profiles_TBDB.df,
            "tbdb_resistance_summary.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(resistance_profiles_WHO.df,
            "who_resistance_summary.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(pairwise_analysis.df,
            "pairwise_analysis.list.csv",
            quote=FALSE, row.names=FALSE)