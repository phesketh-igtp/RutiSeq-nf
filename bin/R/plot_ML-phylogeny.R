
#··············································································#
#··············································································#

#!/usr/bin/env R

# load libraries
packages <- c(
    "ape", "ggtree", "ggrepel", "tidyverse", "patchwork",
    "randomcoloR", "tidytree", "argparse")

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
parser <- ArgumentParser(description = "Plot Maximum-likelihood phylogeny")

# Define command-line options
parser$add_argument("--lineageID",  required = TRUE,
                    help = "Name of the lineage - short (i.e. L4.1, or L4.3)")

parser$add_argument("--rlibrary", required = TRUE,
                    help = "Path to directory containing R scripts/functions")

# Parse command-line arguments
args <- parser$parse_args()

# Import the function for creating the palette
source(paste(args$rlibrary,
             "/functions/isolate_cluster_ancestal_sequence.R", sep = ""))
source(paste(args$rlibrary, "/functions/create-ggtree-palette.R", sep = ""))

#··············································································#
#··············································································#

## Test commands
#lineageID <- "lineage4.8"
#tree <- read.tree("lineage4.8_ML.contree")
#rlibrary="/home/phesketh/Documents/GitHub/TBSEQ.cat-nf/bin/R"
#rlibrary="/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R"
#source(paste(rlibrary, "/functions/isolate_cluster_ancestal_sequence.R", sep=""))
#source(paste(rlibrary, "/functions/create-ggtree-palette.R", sep=""))

#··············································································#
#··············································································#

lineageID <- args$lineageID

lineage <- gsub("lineage", "L", lineageID)

## Import the CSV file for clusters
clusters <- read.delim("clusters.tsv", header = TRUE)

### Read trees
tree <- read.tree("snp.contree")  # Assuming the tree file is in Newick format
tree_rooted <- root(tree, "MTB_anc", resolve.root=TRUE, edgelabel=TRUE)

# Filter rows based on the pattern in any of the specified columns
filtered_clusters <- clusters |> filter(lineage == lineageID)

#··············································································#
#··············································································#

# Create tree specific metadata from clusters

# parsing tree tip names
tree.tips <- tree$tip.label

# Get the custer information only for tips in the tree
tree.clusters <- clusters %>% filter(Sample %in% tree.tips) |> distinct()

# For some reason ggtree wont read the first collumn and so you need to use
## the second column to call the tip names
tree.clusters <- tree.clusters |>
  mutate(Tip_lable = Sample)

# Create the additional rows as a data frame
additional_rows <- data.frame(Sample = c("H37Rv", "MTB_anc"),
  Tip_lable = c("H37Rv", "MTB_anc")
  )
# Find missing columns in additional_rows compared to tree.clusters
missing_cols <- setdiff(names(tree.clusters), names(additional_rows))
# Add missing columns with NA (as character) to additional_rows
additional_rows_full <- additional_rows
for (col in missing_cols) {
  additional_rows_full[[col]] <- NA
}
# Reorder columns to match tree.clusters
additional_rows_full <- additional_rows_full[, names(tree.clusters)]

# Combine
tree.clusters <- rbind(tree.clusters, additional_rows_full) |> 
  distinct()

# Create a matrix of the clusters
tree.clusters.df <- tree.clusters |>
  select(!lineage) |>
  tibble::column_to_rownames(var = "Sample")  # Set Sample as row names

# Replace unclustered values with NA so that there is no legend for those values and improves
# readability of the trees
tree.clusters.df <- tree.clusters.df |>
  mutate(across(everything(), ~replace(., . == "singleton", NA))) |>
  mutate(across(everything(), ~replace(., . == "singleton/singleton/singleton", NA))
  )

tree.clusters <- tree.clusters |>
  mutate(across(everything(), ~replace(., . == "singleton", NA))) |>
  mutate(across(everything(), ~replace(., . == "singleton/singleton/singleton", NA))
  )

#··············································································#
#··············································································#

# Create tree specific metadata palette
color_palette <- create_tree_palette(
        input=tree.clusters.df, lin = lineage)
color_palette <- trimws(color_palette)
color_palette

#··············································································#
#··············································································#

outgroups <- c("MTB_anc", "H37Rv")

#··············································································#
#··············································································#

#·············· Export RData for using in later plots ··············#

save.image(file = paste0(lineageID, ".contree.Rdata"))
