library(dplyr,      quietly = TRUE)
library(tidyverse,  quietly = TRUE)
library(argparse,   quietly = TRUE)
library(openxlsx,   quietly = TRUE)

#··············································································#
#··············································································#

# Initialize the argument parser
parser <- ArgumentParser(description = "Create summary XLSX file for results")

# Define command-line options
parser$add_argument("--summary",     required=TRUE,help="Path to the summary file")
parser$add_argument("--who_res", required=TRUE,help="TB-Profiler WHO results")
parser$add_argument("--tbdb_res",  required=TRUE,help="TB-Profiler TBDB results")
parser$add_argument("--clusters",    required=TRUE,help="Path to processed cluster file")
parser$add_argument("--output", required=TRUE,help="Path of output file")
parser$add_argument("--rlibrary", required=TRUE,help="Path to directory containing R scripts and functions")

# Parse command-line arguments
# Parse the arguments
args <- parser$parse_args()

# Import the function for creating the palette
#source("/home/phesketh/Documents/GitHub/TBSEQ.cat-nf/bin/R/functions/dictionary-rename.R")
#rlibrary <- "/home/phesketh/Documents/GitHub/TBSEQ.cat-nf/bin/R"
rlibrary <- args$rlibrary
source(paste(args$rlibrary, "/functions/dictionary-rename.R", sep=""))

#··············································································#
#··············································································#

# Import all dataframes

#summary   <- read.delim("sequencing_summary.csv", header=T, sep=";",check.names = FALSE);who_res   <- read.delim("who_resistance_summary.csv", header=T, sep=";",check.names = FALSE);tbdb_res  <- read.delim("tbdb_resistance_summary.csv", header=T, sep=";",check.names = FALSE);clusters  <- read.delim("processed_clusters.tsv", header=T, sep="\t",check.names = FALSE);output="test"
summary   <- read.delim(args$summary,  header=T, sep=";",  check.names = FALSE)
who_res   <- read.delim(args$who_res,  header=T, sep=";",  check.names = FALSE)
tbdb_res  <- read.delim(args$tbdb_res, header=T, sep=";",  check.names = FALSE)
clusters  <- read.delim(args$clusters, header=T, sep="\t", check.names = FALSE)

#··············································································#
#··············································································#

# Create the final results summary of the dataframe
who_res.min <- who_res |> select(Sample=sample,`DR type`)
head(summary);head(who_res)
summary.tmp <- left_join(summary, who_res.min)
summary_xlsx <- dictionary_rename(df = summary.tmp,
                                dict_path = paste(rlibrary,
                                        "/dict/xlsx_sheet1.dict.csv",
                                        sep=""))

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
saveWorkbook(wb, file = output, overwrite = TRUE)
