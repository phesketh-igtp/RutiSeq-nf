# load libraries

    library(optparse)
    library(tidyverse)

#··········································································································#
#··························· DEFINE INPUT AND OUTPUT FLAGS FOR THE SCRIPT ·································#
#··········································································································#

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
    analysis_summary <- args$summary
    who_resistance <- args$who_res
    tbdb_resistance <- args$tbdb_res
    pairwise_clusters <- args$clusters
    pairwise_matrix <- args$matrices
    output_file <- args$output

#··········································································································#
#·················································· MAIN ··················································#
#··········································································································#

# Create the
