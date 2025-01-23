fasta_to_dataframe <- function(fasta) {
# Function to create a data frame from a FASTA file
    
# Read sequences from the FASTA file
sequences <- read.fasta(file = fasta, as.string = TRUE, forceDNAtolower = FALSE)
    
# Get sequence names and sequences
seq_names <- names(sequences)
sequence_list <- lapply(sequences, function(seq) unlist(strsplit(seq, "")))
    
# Determine the maximum sequence length
max_length <- max(sapply(sequence_list, length))
    
# Align sequences by padding shorter sequences with NA
aligned_sequences <- lapply(sequence_list, function(seq) {
    c(seq, rep(NA, max_length - length(seq)))
})
    
# Combine aligned sequences into a data frame
result <- as.data.frame(do.call(cbind, aligned_sequences))
    
# Assign sequence names as column names
colnames(result) <- seq_names
    
# Return the data frame
return(result)

}