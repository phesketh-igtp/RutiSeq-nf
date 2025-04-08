#' negative-control-compile.R
#' @author: Poppy J Hesketh Best
#' @date: 2025-04-07
#' @version: 1.0.0
#' @description:
#'  Take all the outputs from the negative control workflow part of 
#' RutiSeq-nd pipeline, and compiles them into a single excell workbook
#' @input:
#'  - negative-controls.stats.csv
#'  - negative-controls.k2.report.csv
#'  - tbprofile.results.txt
#'  - Strain_Classification.tab
#'  - Mapping_and_Variant_Statistics.tab
#' @output: negative-control-report.xlsx
#' @dependencies: base-R, tidyverse.
#' @changelog:
#'      v1.0.0-2025-04-07: Initial version of the script

#··············································································#
#··············································································#

# Load the libraries

library(tidyverse, quietly = TRUE)
library(openxlsx, quietly = TRUE)

## Import the dataframes from the negative control analysis

read_stats.df <- read.delim("negative-controls.stats.csv", sep = ',')
reads_class.df <- read.delim("negative-controls.k2.report.csv", sep = ',')
tbprofiler.df <- read.delim("tbprofiler.txt", sep = '\t')
mtbseq_class.df <- read.delim("Strain_Classification.tab",
                              header = FALSE, sep = '\t')
mtbseq_stats.df <- read.delim("Mapping_and_Variant_Statistics.tab",
                              header = FALSE, sep = '\t')

## Create xlsx workbook
wb <- createWorkbook()
addWorksheet(wb, "Summary")
addWorksheet(wb, "TB-Profiler")
addWorksheet(wb, "MTBSeq-Statistics")
addWorksheet(wb, "MTBSeq-Classification")

#··············································································#
#··············································································#

## Add headers to "MTBSeq-Statistics"

# create the headers
mtbseq_stats.headers <- c(
  "Date",
  "SampleID",
  "LibraryID",
  "FullID",
  "Total Reads",
  "Mapped Reads",
  "% Mapped Reads",
  "Genome Size",
  "Genome GC",
  "(Any) Total Bases",
  "% (Any) Total Bases",
  "(Any) GC-Content",
  "(Any) Coverage mean",
  "(Any) Coverage median",
  "(Unambiguous) Total Bases",
  "% (Unambiguous) Total Bases",
  "(Unambiguous) GC-Content",
  "(Unambiguous) Coverage mean",
  "(Unambiguous) Coverage median",
  "SNPs",
  "Deletions",
  "Insertions",
  "Uncovered",
  "Substitutions (Including Stop Codons)"
)
mtbseq_class.headers <- c(
  "Date",
  "SampleID",
  "LibraryID",
  "FullID",
  "Homolka species",
  "Homolka lineage",
  "Homolka group",
  "Quality",
  "Coll lineage (branch)",
  "Coll lineage_name (branch)",
  "Coll quality (branch)",
  "Coll lineage (easy)",
  "Coll lineage_name (easy)",
  "Coll quality (easy)",
  "Beijing lineage (easy)",
  "Beijing quality (easy)"
)

# add the headers to the dataframes
colnames(mtbseq_stats.df) <- mtbseq_stats.headers
colnames(mtbseq_class.df) <- mtbseq_class.headers

#··············································································#
#··············································································#

## Create the long string of the top 5 taxa for each negative control with their
## corresponding percentages

# Mutate to create a column with the tax and the percentage



#··············································································#
#··············································································#

## ## Add the dataframes to the workbook
writeData(wb, "TB-Profiler", tbprofiler.df)
writeData(wb, "MTBSeq-Statistics", mtbseq_stats.final)
writeData(wb, "MTBSeq-Classification", mtbseq_class.final)
