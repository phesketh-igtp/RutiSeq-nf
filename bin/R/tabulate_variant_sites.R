# load libraries
packages <- c("argparse", "tidyverse", "dplyr", "tidyr", "seqinr")

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

# Create a parser for the script
parser <- ArgumentParser(description = "Script to process MTBseq and TBProfiler data")

# Define arguments
parser$add_argument("--cluster",    required=TRUE, help="")
parser$add_argument("--fasta",      required=TRUE, help="")
parser$add_argument("--positions",  required=TRUE, help="")
parser$add_argument("--H37Rv",required=TRUE, help="")

# Parse the arguments
args <- parser$parse_args()

#··············································································#
#··············································································#

# Source the necessaryfunctions
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

#··············································································#
#··············································································#

clusterID       <- args$cluster
fasta_df        <- fasta_to_dataframe(paste(args$fasta))
positions_raw   <- read.delim(args$positions, header = FALSE) |> select(positions = V1)
H37Rv_genes     <- read.delim(args$H37Rv, header = TRUE, check.names = FALSE)

#fasta_df <- fasta_to_dataframe(paste("n1-L4.3_refseq.fasta"));positions_raw <- read.delim("n1-L4.3_genomic_positions", header = FALSE) |> select(positions = V1); clusterID <- "n1-L4.3"
#H37Rv_genes     <- read.delim("H37Rv.genes.txt", header = TRUE, check.names = FALSE)

#··············································································#
#··············································································#

# Join the possitions with the allels
positions.wide <- cbind(fasta_df, positions_raw)

# Reshape the df from wide to long and create the counts data
positions_df <- positions.wide |> 
                pivot_longer(!positions, 
                names_to = "Sample", values_to = "allel") |>
                mutate(cluster = clusterID) |> ungroup() |>
                mutate(
                        genome_type = case_when(
                        Sample == "H37Rv" ~ "H37Rv",
                        Sample == "MTB_anc" ~ "MTB_anc",
                        TRUE ~ "samples"  # Default case for all other samples
                        )
                    )
                            
positions_counts <- positions_df |>
            group_by(positions, allel, cluster,genome_type) |>
            summarize(count = n(), .groups = 'drop') |>
            select(positions,allel,cluster,genome_type,freq=count)

positions_list <- positions_counts |>
            select(positions) |>
            distinct()

# Annotate the genomic possitions (bidirectionallity)
annotated_positions <- positions_list |>
    rowwise() |>
    mutate(
    match = list(H37Rv_genes |>
            filter(
            (start <= positions & stop >= positions) |  # Normal range
            (start >= positions & stop <= positions)    # Reversed range
                )
            )
        ) |>
    unnest(match, keep_empty = TRUE)

# merge annotations with df for final df to export
positions_counts.final  <- left_join(positions_counts, annotated_positions)
positions_df.final      <- left_join(positions_df, annotated_positions)

#··············································································#
#··············································································#

# Export the dataframes as csv files

write.table(positions_df.final, quote=FALSE, sep = ";",
            file = paste0(clusterID,".variant-positions.csv",sep=""),
            col.names = TRUE, row.names = FALSE)

write.table(positions_counts.final, quote=FALSE, sep = ";",
            file = paste0(clusterID,".variant-positions.counts.csv",sep=""),
            col.names = TRUE, row.names = FALSE)

#··············································································#
#··············································································#