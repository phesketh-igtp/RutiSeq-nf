process_cluster <- function(clusterID,df) {
  
  # Filter the selection based on clusterID
  selection <- df |>
    filter(cluster == clusterID) |>
    select(SampleID) |>
    deframe()
  
  # Visualize and find the most recent common ancestor (MRCA)
  p <- ggtree(tree) + geom_tiplab() #+ geom_nodelab()
  mrca_node <- MRCA(p, selection)
  viewClade(p, MRCA(p, selection))
  
  # Subset the tree to the cluster MRCA
  mrca_tree <- tree_subset(tree, node=mrca_node, levels_back=0)
  
  # Get the nodeID for the very first node in the subsetted tree
  ancestor_nodeID <- mrca_tree$node.label[1]
  
  # Filter the fast sequence to isoalte the ancestral sequence
  filtered_sequence <- fasta_sequences[ancestor_nodeID]
  
  # Produce output fast file name
  output_fasta <- paste("ancestors/",clusterID,".ancestor.fasta", sep="")
  writeXStringSet(filtered_sequence, filepath = output_fasta)
  
  
  # Convert the filtered_sequence to a character string (e.g., "ATGCCG...")
  seq_str <- as.character(filtered_sequence)
  
  # Create a data frame with position (1 to n) and the corresponding base pair
  seq_df <- data.frame(
    Node = rep(ancestor_nodeID, nchar(seq_str)),  # Repeating the name for each base
    Cluster = rep(clusterID,nchar(seq_str)),
    Position = 1:nchar(seq_str),                  # Position is just 1 to length of the sequence
    Allel = unlist(strsplit(seq_str, split = ""))  # Split the sequence into individual base pairs
  )
  
  # Display the reshaped data frame
  head(seq_df)
  output_tsv <- paste("ancestors/",clusterID,".ancestor.positions", sep="")
  write.table(seq_df, output_tsv, 
              quote = FALSE, row.names = FALSE)
  
}