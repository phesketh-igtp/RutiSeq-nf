# load libraries

library(optparse)
library(tidyverse)

#··········································································#
#··········· DEFINE INPUT AND OUTPUT FLAGS FOR THE SCRIPT ·················#
#··········································································#

# Define the command-line options
option_list <- list(
    make_option(c("--summary"), type="character", help="Path to analysis summary file", metavar="FILE"),
    make_option(c("--who_res"), type="character", help="Path to WHO resistance file", metavar="FILE"),
    make_option(c("--tbdb_res"), type="character", help="Path to TBDB resistance file", metavar="FILE"),
    make_option(c("--clusters"), type="character", help="Path to pairwise clusters file", metavar="FILE"),
    make_option(c("--matrices"), type="character", help="Path to pairwise matrix file", metavar="FILE"),
    make_option(c("--output"), type="character", help="Output file name", metavar="FILE")
)

# Parse the options
parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

# Validate required arguments
if (is.null(args$summary) || is.null(args$who_res) || is.null(args$tbdb_res) ||
    is.null(args$clusters) || is.null(args$matrices) || is.null(args$output)) {
    print_help(parser)
    stop("All arguments are required.")
}

# Assign arguments to variables
analysis_summary <- read_csv(args$summary, header = TRUE)
who_resistance <- read_csv(args$who_res, header = TRUE)
tbdb_resistance <- read_csv(args$tbdb_res, header = TRUE)
pairwise_clusters <- read_csv(args$clusters, header = TRUE)
pairwise_matrix <- read_csv(args$matrices, header = TRUE)
output_file <- args$output

# Import the dataframes

#··········································································#
#································ MAIN ····································#
#··········································································#

# Create the
