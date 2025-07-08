#!/usr/bin/env R

library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(readr, quietly = TRUE)

# Import dataframes
main_lineages <- readr::read_delim("selected_main-lineage_split.list",
                                   col_names = FALSE, delim = ",")

colnames(main_lineages) <- "selected_main_lineage"

run_ids <- readr::read_delim("run_sample_ids.txt",
                             col_names = FALSE, delim = ",")

colnames(run_ids) <- "SampleID"

meta <- readr::read_delim("pairwise_analysis.list.csv",
                          col_names = FALSE) |>
  distinct()

colnames(meta) <- c("SampleID", "main_lineage", "sub_lineage")

# Filter out any sample that contains 'CN-'
meta <- meta |>
  filter(!is.na(main_lineage) & !str_detect(main_lineage, ";")) |>
  filter(!str_detect(SampleID, "CN-"))

# Select only the main lineages that are in the selected_main_lineage list
filtered_meta <- meta |>
  select(lineage = main_lineage, SampleID)

skipped_lineages_to_few <- filtered_meta |>
  group_by(lineage) |>
  count() |>
  filter(n < 5)

# Get unique lineages from run_ids
filtered_lineages_pass <- filtered_meta |>
  filter(SampleID %in% run_ids$SampleID) |>
  filter(!lineage %in% skipped_lineages_to_few$lineage) |>
  count(lineage) |>
  select(lineage) |>
  distinct()  # Use distinct() instead of unique() for dplyr consistency

# Partition filtered_meta based on filtered_lineages_pass
filtered_meta_forward <- filtered_meta %>%
  filter(lineage %in% filtered_lineages_pass$lineage)

filtered_meta_skip <- filtered_meta %>%
  filter(!lineage %in% filtered_lineages_pass$lineage)

# Export the csv
write.csv(filtered_meta_forward, "final.lineage_samples_tuple.csv",
          quote = FALSE,
          row.names = FALSE
          )

write.csv(filtered_meta_skip, "final.skipped-lineages_tuple.csv",
          quote = FALSE,
          row.names = FALSE
          )
