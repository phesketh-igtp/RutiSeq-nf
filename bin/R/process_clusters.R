library(tidyverse,  quietly = TRUE, verbose = FALSE)
library(dplyr,      quietly = TRUE, verbose = FALSE)
library(argparse,   quietly = TRUE, verbose = FALSE)

#··············································································#
#··············································································#

# Initialize the argument parser
parser <- ArgumentParser(description = "Process cluster and summary files")

# Define command-line options
parser$add_argument("--clusters", required = TRUE, help = "Path to the clusters file")
parser$add_argument("--summary", required = TRUE, help = "Path to the summary file CSV")

# Parse the arguments
args <- parser$parse_args()

#··············································································#
#··············································································#

# Import the dataframes
raw_clust <- read.delim(args$clusters, header=TRUE, sep = "\t")
raw_clust <- raw_clust |> 
                    pivot_wider(names_from = distance, values_from = group)
colnames(raw_clust) <- c("lineage","genomes","dSNP_5","dSNP_10","dSNP_15")

sampleIDs <- read.delim(args$summary, header=TRUE, sep = ";") |>
                select(Sample) |>
                mutate(genomes=Sample) |>
                separate_wider_delim(genomes, delim = "_", names = c("genomes", "library")) |>
                select(SampleID=Sample,genomes)

#··············································································#
#··············································································#

# Simplify the names of the groups
raw_clust$dSNP_5<-gsub("group_","n",as.character(raw_clust$dSNP_5))
raw_clust$dSNP_5<-gsub("ungrouped","nX",as.character(raw_clust$dSNP_5))
raw_clust$dSNP_10<-gsub("group_","i",as.character(raw_clust$dSNP_10))
raw_clust$dSNP_10<-gsub("ungrouped","iX",as.character(raw_clust$dSNP_10))
raw_clust$dSNP_15<-gsub("group_","v",as.character(raw_clust$dSNP_15))
raw_clust$dSNP_15<-gsub("ungrouped","vX",as.character(raw_clust$dSNP_15))

# Create grouping IDs that contain the Lineage information
transmission <- raw_clust |> mutate(SNP_nd5.id10.vd15 = paste0(dSNP_5,".",dSNP_10,".",dSNP_15, "-", lineage))
transmission <- transmission |> mutate(SNP_d5_L = paste0(dSNP_5,"-", lineage))
transmission <- transmission |> mutate(SNP_d10_L = paste0(dSNP_10,"-", lineage))
transmission <- transmission |> mutate(SNP_d15_L = paste0(dSNP_15,"-", lineage))

transmission$SNP_nd5.id10.vd15<-gsub("lineage","L",as.character(transmission$SNP_nd5.id10.vd15))
transmission$SNP_d5_L<-gsub("lineage","L",as.character(transmission$SNP_d5_L))
transmission$SNP_d10_L<-gsub("lineage","L",as.character(transmission$SNP_d10_L))
transmission$SNP_d15_L<-gsub("lineage","L",as.character(transmission$SNP_d15_L))

transmission <- transmission |> distinct()

# Get the real sample ID instead of the truncated one
transmission.master.clusters <- left_join(transmission, sampleIDs) |>
            select(SampleID,lineage,dSNP_5,dSNP_10,
            dSNP_15,SNP_nd5.id10.vd15,SNP_d5_L,SNP_d10_L,SNP_d15_L) |>
            distinct()

write.table(transmission.master.clusters, "processed_clusters.tsv", 
                        sep = "\t", row.names = FALSE, quote = FALSE)