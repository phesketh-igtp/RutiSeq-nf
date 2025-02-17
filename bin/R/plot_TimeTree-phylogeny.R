library(ape)
library(ggtree)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(randomcoloR)
library(tidytree)
library(argparse)
library(Biostrings)

#··············································································#
#··············································································#
# Initialize the argument parser
parser <- ArgumentParser(description = "Plot Maximum-likelihood TimeTree phylogeny")

# Define command-line options
parser$add_argument("--timetree", required=TRUE, help="Path to the tree file (contre format)")
parser$add_argument("--clusters", required=TRUE, help="Path to the cluster file")
parser$add_argument("--lineageID",required=TRUE, help="Name of the lineage - short (i.e. L4.1, or L4.3)")
parser$add_argument("--fasta",    required=TRUE, help="Path to file containing the ancestral node outputs FASTA file")
parser$add_argument("--rlibrary", required=TRUE, help="Path to directory containing R scripts and functions")

# Parse command-line arguments
args <- parser$parse_args()

# Import the function for creating the palette
source(paste(args$rlibrary, "/functions/isolate_cluster_ancestal_sequence.R", sep=""))
source(paste(args$rlibrary, "/functions/create-ggtree-palette.R", sep=""))

#··············································································#
#··············································································#

## Test commands
#lineageID <- "lineage4.8"
#clusters <- read.delim("processed_clusters.tsv", header=T) |> select(SampleID,SNP_d5_L,SNP_d10_L,SNP_d15_L,SNP_nd5.id10.vd15)
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
clusters <- read.delim(args$cluster, header = TRUE) |> 
        select(SampleID,SNP_d5_L,SNP_d10_L,SNP_d15_L,SNP_nd5.id10.vd15)

### Read trees
tree <- read.nexus(args$timetree)  # Assuming the tree file is in Newick format
tree_rooted <- root(tree, "MTB_anc", resolve.root=TRUE, edgelabel=TRUE)

# Filter clusters to the correct analysis group lineage
pattern <- paste0("-", lineage, "$")
# Filter rows based on the pattern in any of the specified columns
filtered_clusters <- clusters[
        grepl(pattern, clusters$SNP_d5_L) |
        grepl(pattern, clusters$SNP_d10_L) |
        grepl(pattern, clusters$SNP_d15_L) |
        grepl(pattern, clusters$SNP_nd5.id10.vd15),
        ]; clusters <- filtered_clusters; rm(filtered_clusters)

#··············································································#
#··············································································#

# Create tree specific metadata from clusters

# parsing tree tip names
tree.tips <- tree$tip.label

# Get the custer information only for tips in the tree
tree.clusters <- clusters %>% filter(SampleID %in% tree.tips) |> distinct()

# For some reason ggtree wont read the first collumn and so you need to use
## the second column to call the tip names
tree.clusters <- tree.clusters |> 
        mutate(Tip_lable=SampleID) |>
        select(SampleID,Tip_lable,
                SNP_d5_L,SNP_d10_L,SNP_d15_L,SNP_nd5.id10.vd15
                )

# Create the additional rows as a data frame
additional_rows <- data.frame(SampleID = c("H37Rv", "MTB_anc"),Tip_lable = c("H37Rv", "MTB_anc"),
        SNP_d5_L = NA,SNP_d10_L = NA,SNP_d15_L = NA,SNP_nd5.id10.vd15 = NA)
# Add the new rows to the tree.clusters data frame
tree.clusters <- rbind(tree.clusters, additional_rows)

# Create a matrix of the clusters
tree.clusters.df <- tree.clusters |>
        select(SampleID, d5 = SNP_d5_L, d10 = SNP_d10_L, 
                d15 = SNP_d15_L, d5.10.15 = SNP_nd5.id10.vd15) |>
        tibble::column_to_rownames(var = "SampleID")  # Set SampleID as row names

# Replace unclustered values with NA so that there is no legend for those values and improves
# readability of the trees
pattern <- paste0("nX-", lineage)
        tree.clusters.df$d5[tree.clusters.df$d5 == pattern] <- NA
pattern <- paste0("iX-", lineage)
        tree.clusters.df$d10[tree.clusters.df$d10 == pattern] <- NA
pattern <- paste0("vX-", lineage)
        tree.clusters.df$d15[tree.clusters.df$d15 == pattern] <- NA
pattern <- paste0("nX.iX.vX-", lineage)
        tree.clusters.df$d5.10.15[tree.clusters.df$d5.10.15 == pattern] <- NA

#··············································································#
#··············································································#

# Create tree specific metadata palette
color_palette <- create_tree_palette(
        input=tree.clusters.df, lin = lineage)
color_palette <- trimws(color_palette)        
color_palette

#··············································································#
#··············································································#

# Final modifications to the plotting dataframes
tree.clusters.df.tmp <- tree.clusters.df |> select(d5,d10,d15)

# create a dummy rf for the outgroups
outgroups <- data.frame(sampleID = c("MTB_anc", "H37Rv"),
                        d5 = c(NA, NA),
                        d10 = c(NA, NA),
                        d15 = c(NA, NA)) %>% 
        column_to_rownames(var = "sampleID")

tree.clusters.df <- rbind(tree.clusters.df.tmp, outgroups)
rm(tree.clusters.df.tmp); rm(outgroups)
outgroups <- c("MTB_anc","H37Rv")

#··············································································#
#··············································································#

# Timetree corrected phylogeny

tree.p.c <- ggtree(tree_rooted, 
                linewidth=0.1, 
                layout = "circular",
                branch.length = 'none'
                ) %<+% 
        tree.clusters +
        geom_tiplab(aes(
                label=Tip_lable,
                color = label %in% outgroups),
                size = 2, 
                align = TRUE
                ) +
        scale_color_manual(
                values=c("#000000", "#FF0000")
        )

p3 <- gheatmap(tree.p.c, 
                tree.clusters.df,
                offset = 10, 
                width = 0.3,
                colnames_angle = 45, 
                colnames_offset_y = -0.8,
                font.size = 5, 
                color = "#3a3a3a",
                ) +
        scale_fill_manual(
                values = color_palette,
                na.value = "white", 
                name = "Unique\nclusterID\n(with lineage)") +  # Default color for unmatched values
        ggtitle(paste0("IQ-Tree ML TimeTree | Branch-length ignored | ", lineage)) +
        theme(legend.position = "none")
p3


#··············································································#
#··············································································#

fasta_sequences <- readDNAStringSet(args$fasta)

d5.tree.clusters <- clusters |> 
                filter(SampleID %in% tree.tips) |> 
                filter(!grepl("nX-", SNP_d5_L)) |>
                distinct() |>
                select(SampleID,cluster=SNP_d5_L)

d5.tree.clusters.groups <- d5.tree.clusters |> 
                select(cluster) |> 
                distinct()

d5.tree.clusters.groups.deframed <- deframe(d5.tree.clusters.groups)

# Apply the function to each clusterID in the dataframe
for (unique_cluster in d5.tree.clusters.groups.deframed) {
        process_cluster(
                clusterID = unique_cluster,
                df    = d5.tree.clusters
                )
        }

#··············································································#
#··············································································#

#·············· Export trees ··············#

pdf(file = paste0(lineageID,"_TimeTree.contree.pdf"))
p3
dev.off()

#·············· Export RData for using in later plots ··············#

save.image(file = paste0(lineageID,".time-tree.RData"))
