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
parser$add_argument("--lineageID",  required=TRUE,help="Name of the lineage - short (i.e. L4.1, or L4.3)")

# Parse command-line arguments
args <- parser$parse_args()

#··············································································#
#··············································································#

lineageID <- args$lineageID
lineage <- gsub("lineage", "L", lineageID)

## Import the CSV file for clusters
clusters <- read.delim("clusters.tsv", header = TRUE) |> 
        select(SampleID,SNP_d5_L,SNP_d10_L,SNP_d15_L,SNP_nd5.id10.vd15)

### Read trees
tree <- read.tree("snp.contree")  # Assuming the tree file is in Newick format
tree_rooted <- root(tree, "MTB_anc", resolve.root=TRUE, edgelabel=TRUE)

# Filter clusters to the correct analysis group lineage
pattern <- paste0("-", lineage, "$")

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

# Rectangular trees

# Plot base trees - rectangular and circular
tree.p.r <- ggtree(tree_rooted, linewidth=0.5,
                layout = "rectangular"
                ) %<+% 
            tree.clusters +
            geom_tiplab(aes(label=Tip_lable,
                color = label %in% outgroups),
                size = 2.5, align = TRUE
                ) +
            scale_color_manual(values=c("#000000", "#FF0000")
            )

# Circular (dendrogram) trees
tree.p.c <- ggtree(tree_rooted,
                linewidth=0.1,
                layout = "circular",
                branch.length = 'none'
                ) %<+%
            tree.clusters +
            geom_tiplab(aes(label=Tip_lable,
                color = label %in% outgroups),
                size = 2.5, align = TRUE
                ) +
            scale_color_manual(values=c("#000000", "#FF0000")
            )

#··············································································#
#··············································································#

#·············· Export trees ··············#

pdf(file = paste0(lineageID,"_ML.contree.pdf"))
tree.p.r
tree.p.c
dev.off()

#·············· Export RData for using in later plots ··············#

save.image(file = paste0(lineageID,".contree.Rdata"))