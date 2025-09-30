#!/usr/bin/env R

# load libraries
packages <- c(
    "dplyr", "tidyverse", "argparse", "openxlsx")

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
#··············································································#

# Initialize the argument parser
parser <- ArgumentParser(description = "Create summary XLSX file for results")

# Define command-line options
parser$add_argument("--output", required=TRUE,help="Path of output file")
parser$add_argument("--rlibrary", required=TRUE,help="Path to directory containing R scripts and functions")


# Parse command-line arguments
# Parse the arguments
args <- parser$parse_args() 

# Import the function for creating the palette
#source("/home/phesketh/Documents/GitHub/TBSEQ.cat-nf/bin/Rfunctions/dictionary-rename.R")
#source("/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R/functions/dictionary-rename.R")
#rlibrary="/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R"
rlibrary <- args$rlibrary
source(paste(rlibrary, "/functions/dictionary-rename.R", sep=""))

OUT <- args$output

#··············································································#
#··············································································#

# Import all dataframes
summary <- read.delim("sequencing_summary.csv",  header = TRUE,
                      sep = ",",  check.names = FALSE)

who_res <- read.delim("who_resistance_summary.csv",  header = TRUE,
                      sep = ",",  check.names = FALSE) |>
  select(Sample = sample, `DR type`)

tbdb_res <- read.delim("tbdb_resistance_summary.csv", header = TRUE,
                       sep = ",",  check.names = FALSE)

clusters <- read.delim("processed_clusters.tsv", header = TRUE,
                       sep = "\t", check.names = FALSE)

who_corr <- read.delim("tbprofiler-DR-corrections.csv", header = TRUE,
                       sep = ";", check.names = FALSE) |>
  select(Sample = sample,
         corr_DRType = DRT,
         corr_DRTypeExt = DRT_ext)

#··············································································#
#··············································································#

# Create the final results summary of the dataframe

who_res_min <- left_join(who_res, who_corr)

summary.tmp <- left_join(summary, who_res_min)

summary_xlsx <- dictionary_rename(df = summary.tmp,
                                  dict_path = paste(rlibrary, "/dict/xlsx_sheet1.dict.csv",
                                  sep="")
                                  )

#··············································································#
#··············································································#

# Assemble the XLSX spreadsheet

# Create a new workbook
wb <- createWorkbook()

# Add each sheet to the workbook
addWorksheet(wb, "sequencing_summary")
writeData(wb, "sequencing_summary", summary_xlsx, rowNames=FALSE)

addWorksheet(wb, "resistance_table_who")
writeData(wb, "resistance_table_who", who_res, rowNames=FALSE)

addWorksheet(wb, "resistance_table_tbdb")
writeData(wb, "resistance_table_tbdb", tbdb_res, rowNames=FALSE)

addWorksheet(wb, "transmission_clusters")
writeData(wb, "transmission_clusters", clusters, rowNames=FALSE)

# Save the workbook
saveWorkbook(wb, file = OUT, overwrite = TRUE)

#··············································································#
#··············································································#
