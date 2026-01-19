#!/usr/bin/env R

# load libraries
packages <- c(
        "ape", "ggtree", "ggrepel", "patchwork", "randomcoloR",
        "tidytree", "argparse", "Biostrings", "tibble", "dplyr", "stringr")

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

process_cluster <- function(clusterID,df) {

  library(ggtree)
  library(ape)

  # Filter the selection based on clusterID
  selection <- df |>
    filter(cluster == clusterID) |>
    select(SampleID) |>
    deframe()

  # Visualize and find the most recent common ancestor (MRCA)
  p <- ggtree(tree) + geom_tiplab() + geom_nodelab()
  mrca_node <- MRCA(p, selection)
  viewClade(p, MRCA(p, selection))
  ancestor_nodeID <- p$data |> filter(node == mrca_node) |> pull(label)
  ancestor_nodeID

  # Filter the fast sequence to isoalte the ancestral sequence
  filtered_sequence <- fasta_sequences[ancestor_nodeID]
  filtered_sequence

  # Produce output fast file name
  output_fasta <- paste("ancestors/", clusterID, ".ancestor.fasta", sep = "")
  writeXStringSet(filtered_sequence, filepath = output_fasta, format = "fasta")

  # Convert the filtered_sequence to a character string (e.g., "ATGCCG...")
  seq_str <- as.character(filtered_sequence)

  # Create a data frame with position (1 to n) and the corresponding base pair
  seq_df <- data.frame(
    Node = rep(ancestor_nodeID,
               # Repeating the name for each base
               nchar(seq_str)),
    Cluster = rep(clusterID, nchar(seq_str)),
    # Position is just 1 to length of the sequence
    Position = 1:nchar(seq_str),
    # Split the sequence into individual base pairs
    Allel = unlist(strsplit(seq_str, split = "")
    )
  )

  # Display the reshaped data frame
  head(seq_df)
  output_tsv <- paste("ancestors/", clusterID, ".ancestor.positions", sep = "")
  write.table(seq_df, output_tsv, sep = "\t",
              quote = FALSE, row.names = FALSE)

}

#··············································································#
#··············································································#

## Import the CSV file for clusters
clusters <- read.delim("processed_clusters.tsv", header = TRUE)

### Read trees
tree <- read.nexus("timetree.nexus")  # Assuming the tree file is in Newick format

# Get tree tips
tree.tips <- tree$tip.label

# Filter rows based on the pattern in any of the specified columns
filtered_clusters <- clusters |> filter(SampleID %in% tree.tips) |>
        filter(merged_clusterID != "NA/NA/NA")

#··············································································#
#··············································································#

# Create tree specific metadata from clusters
# For some reason ggtree wont read the first collumn and so you need to use
## the second column to call the tip names
tree.clusters <- filtered_clusters
rownames(tree.clusters) <- tree.clusters$SampleID
dist_col <- tree.clusters |> select(contains("t.")) |> colnames()
dist_col <- dist_col[order(as.numeric(sub("t\\.", "", dist_col)))]

# Create a matrix of the clusters
tree.clusters.df <- tree.clusters |>
        select(dist_col)

#··············································································#
#··············································································#

fasta_sequences <- readDNAStringSet("ancestral_sequences.fasta")

dist_col_lowest <- dist_col[1]

recent.clusters <- tree.clusters |> 
        select(SampleID, all_of(dist_col_lowest)) |>
        mutate(across(all_of(dist_col_lowest), as.character)) |>
        filter(!is.na(.data[[dist_col_lowest]])) |>
        filter(.data[[dist_col_lowest]] != "singleton") |>
        select(SampleID, cluster=all_of(dist_col_lowest)) |>
        distinct()

recent.clusters.deframed <- recent.clusters |>
        select(cluster) |>
        distinct() |> deframe()

# Apply the function to each clusterID in the dataframe
for (unique_cluster in recent.clusters.deframed) {
        process_cluster(
                clusterID = unique_cluster,
                df    = recent.clusters
                )
        }

#··············································································#
#··············································································#
