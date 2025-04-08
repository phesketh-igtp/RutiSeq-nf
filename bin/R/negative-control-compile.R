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
rm(mtbseq_stats.headers);rm(mtbseq_class.headers)

#··············································································#
#··············································································#

## Create the long string of the top 5 taxa for each negative control with their
## corresponding percentages

# Mutate to create a column with the tax and the percentage
reads_class.string <- reads_class.df |>
  filter(num_reads > 10, taxonomy != "root") |>
  mutate(classification = word(taxonomy, -1, sep = ";")) |>
  filter(!classification %in% c("Bacteria", "Virus", 
                                "Eukaryota", "cellular organisms")) |>
    group_by(sampleID, classification) |>
  summarise(percentage = sum(percentage), .groups = "drop") |>
  group_by(sampleID) |>
  arrange(desc(percentage), .by_group = TRUE) |>
  slice_head(n = 5) |>
  mutate(Classification_perc = paste(classification,
                                     " (",percentage,"%)", sep ="")) |>
  select(sampleID,Classification_perc) |> ungroup() |>
  group_by(sampleID) |>
  summarise(Classification_perc_str = str_c(Classification_perc, collapse = "; "))

#··············································································#
#··············································································#

## Summary the read stats for PE-reads
head(read_stats.df)

read_stats.summary <- read_stats.df |>
  group_by(sampleID) |>
  summarise(
    num_seqs = sum(num_seqs),
    across(
      .cols = where(is.numeric) & !matches("num_seqs"),
      .fns = mean,
      .names = "{.col}_avg"
    )
  )

read_stats.final <- left_join(reads_class.string, read_stats.summary)


read_stats.final.colnames <- c("sampleID", "Top 5 classifications (%)",
                                  "PE reads", "Sum length", 
                                  "Min length", "Average length", "Max length", 
                                  "Q1", "Q2", "Q3", "Sum gap", "N50", "N50 number", 
                                  "Q20", "Q30", "Average quality", "GC%", "Sum number ave")

colnames(read_stats.final) <- read_stats.final.colnames

read_stats.final <- read_stats.final |> 
  mutate(
    Flags = case_when(
      `PE reads` >= 500000 ~ "Warning: >500,000 reads",
      `PE reads` >= 100000 ~ "Caution: Between 100,000 reads",
      TRUE ~ NA_character_
    )
  )

read_stats.final <- read_stats.final |> select(
  `sampleID`, `Flags`, `PE reads`, `Sum length`, 
  `Min length`, `Average length`, `Max length`,
  `Q1`, `Q2`, `Q3`, `Sum gap`, `N50`,
  `Q20`, `Q30`, `Average quality`, `GC%`, `Top 5 classifications (%)`)


#··············································································#
#··············································································#

## Rename TB-Profiler outouts

tbprofiler.colnames <- c(
  "sampleID", "Main lineage", "Sub-lineage", "Spoligotype", "DF type (TBDB)", 
  "Target median depth", "% read mapped", "No. reads mapped", "No. DR variants", 
  "No. other variants", "RIF", "ISO", "ETH", "PYRAZ", "MOXI", "LEVO", 
  "BED", "DEL", "PRE", "LINEZ", "STREP", "AMIK", "KAN", "CAP", 
  "CLOF", "ETHION", "PAS", "CYCLSER"
)

colnames(tbprofiler.df) <- tbprofiler.colnames
tbprofiler.df <- tbprofiler.df |> mutate(sampleID = gsub("tbdb-", "", sampleID))

#··············································································#
#··············································································#

## Create xlsx workbook
wb <- createWorkbook()
addWorksheet(wb, "Read statistics")
addWorksheet(wb, "Read taxonomy")
addWorksheet(wb, "TB-Profiler")
addWorksheet(wb, "MTBSeq-Statistics")
addWorksheet(wb, "MTBSeq-Classification")

## ## Add the dataframes to the workbook
writeData(wb, "Read statistics", read_stats.final)
writeData(wb, "Read taxonomy", reads_class.df)
writeData(wb, "TB-Profiler", tbprofiler.df)
writeData(wb, "MTBSeq-Statistics", mtbseq_stats.df)
writeData(wb, "MTBSeq-Classification", mtbseq_class.df)

# Save the workbook to a file
saveWorkbook(wb, "negative-controls.xlsx", overwrite = TRUE)