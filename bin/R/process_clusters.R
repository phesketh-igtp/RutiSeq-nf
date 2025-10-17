#!/usr/bin/env R

# load libraries
packages <- c("tidyverse", "dplyr")

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

# Import the dataframes
raw_clust <- read.delim("unprocessed_clusters.tsv",
                        header = FALSE, sep = "\t") |>
  distinct()

colnames(raw_clust) <- c("lineage", "dSNP", "genomes", "int_clusterID")

raw_clust <- raw_clust |>
  mutate(dSNP = gsub("dist_", "t=", dSNP))

# Process the sampleIDs
sample_id_df <- read.delim("sequencing_summary.csv",
                           header = TRUE, sep = ",") |>
  select(Sample) |>
  mutate(genomes = Sample) |>
  separate_wider_delim(genomes, delim = "_",
                       names = c("genomes", "library"),
                       too_few = "align_end") |>
  select(SampleID = Sample, genomes)

# Merge and pivot the data
raw_clust <- left_join(raw_clust, sample_id_df) |>
  distinct() |>
  select(SampleID, lineage, dSNP, int_clusterID) |>
  pivot_wider(names_from = dSNP, values_from = int_clusterID)

# Extract and sort t= column names numerically
t_cols <- names(raw_clust)[grepl("^t=", names(raw_clust))]
t_cols_sorted <- t_cols[order(as.numeric(gsub("t=", "", t_cols)))]

# Reorder the full data frame
raw_clust <- raw_clust %>%
  select(SampleID, lineage, all_of(t_cols_sorted))

# Create the merged clusterID
raw_clust <- raw_clust %>%
  unite("merged_clusterID", all_of(t_cols_sorted), sep = "/", remove = FALSE) |>
  select(SampleID, lineage, all_of(t_cols_sorted), merged_clusterID)

#··············································································#
#··············································································#

write.table(raw_clust, "processed_clusters.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)