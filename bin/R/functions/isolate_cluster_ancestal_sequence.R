process_cluster <- function(clusterID,df) {

  library(ggtree)
  library(ape)

  # Filter the selection based on clusterID
  selection <- df |>
    filter(cluster == clusterID) |>
    select(SampleID) |>
    deframe()

  # Visualize and find the most recent common ancestor (MRCA)
  tree_rooted <- root(tree, "MTB_anc",
                      resolve.root = TRUE,
                      edgelabel = TRUE)
  p <- ggtree(tree_rooted) + geom_tiplab() + geom_nodelab()
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