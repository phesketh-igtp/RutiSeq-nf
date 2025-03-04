library(dplyr,      quietly = TRUE)
library(tidyverse,  quietly = TRUE)
library(argparse,   quietly = TRUE)
library(openxlsx,   quietly = TRUE)

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
#summary   <- read.delim("sequencing_summary.csv", header=T, sep=";",check.names = FALSE);who_res   <- read.delim("who_resistance_summary.csv", header=T, sep=";",check.names = FALSE);tbdb_res  <- read.delim("tbdb_resistance_summary.csv", header=T, sep=";",check.names = FALSE);clusters  <- read.delim("processed_clusters.tsv", header=T, sep="\t",check.names = FALSE);output="test"
summary   <- read.delim("sequencing_summary.csv",  header=T, sep=";",  check.names = FALSE)
who_res   <- read.delim("who_resistance_summary.csv",  header=T, sep=";",  check.names = FALSE)
tbdb_res  <- read.delim("tbdb_resistance_summary.csv", header=T, sep=";",  check.names = FALSE)
clusters  <- read.delim("processed_clusters.tsv", header=T, sep="\t", check.names = FALSE)

#··············································································#
#··············································································#

# Create the final results summary of the dataframe
who_res.min <- who_res |> select(Sample=sample,`DR type`)
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
saveWorkbook(wb, file = OUT, overwrite = TRUE)

#··············································································#
#··············································································#
