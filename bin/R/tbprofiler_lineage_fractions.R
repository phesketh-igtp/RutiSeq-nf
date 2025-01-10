cat("Running: R/tbprofiler_lineage_fractions.R\n")

# load libraries
library(argparse); library(tidyverse)
library(dplyr); library(tidyr)

# Create a parser for the script
    parser <- ArgumentParser(description = "Script to process MTBseq and TBProfiler data")

# Define arguments
parser$add_argument("--tbprofiler", required=TRUE, help="Path to TBprofile compiled results")
parser$add_argument("--lineages", required=TRUE, help="Path to lineages tables with fractions")

## parser$add_argument("--additional_args", default=NULL, help="Optional additional arguments from config file")

# Parse the arguments
args <- parser$parse_args()

# Assign arguments to variables
tbprof_file     <- args$tbprofiler
lineage_frac_file <- args$lineages
OUT <- "tbprofiler.lineages.fractions.txt"

# Read the input files
tbprof <- read.delim(tbprof_file, header = TRUE) %>%
    mutate(across(everything(), ~ str_trim(as.character(.)))) %>%
    select(SampleID=sample,main_lineage,sub_lineage) %>% # designate 'clonal' or 'mixed
    mutate(infection_type = ifelse(grepl(';', main_lineage) | grepl(';', sub_lineage), 'mixed', 'clonal'))

tbprof.infec <- tbprof %>% select(SampleID,infection_type)

lineage.frac <- read.delim(lineage_frac_file, , header = TRUE) %>%
                filter(Fraction != 'Fraction')

# Split 'clonal' and 'clonal'
# clonal: add fraction of lineage designation
mixed <- tbprof %>% filter(infection_type == 'mixed')

if (nrow(mixed) == 0) { # Check if 'mixed' is empty

    mixed.frac.5 <- mixed

} else { # Proceed with processing if 'mixed' is not empty

    mixed.frac <- left_join(mixed, lineage.frac) %>%
        select(SampleID, main_lineage, sub_lineage, infection_type, Lineage, Fraction)
    mixed.frac$Fraction <- as.numeric(mixed.frac$Fraction)
    mixed.frac <- mixed.frac %>% mutate(Perc = 100 * Fraction)

    mixed.frac.1 <- mixed.frac %>%
        mutate(Lineage_p = paste0(Lineage, " (", Perc, "%)")) %>%
        mutate(Lineage_p = gsub("lineage", "L", Lineage_p))

    mixed.frac.2 <- mixed.frac.1 %>%
        select(SampleID, sub_lineage, Lineage_p) %>%
        mutate(sub_lineage = gsub("lineage", "L", sub_lineage)) %>%
        mutate(Lineage = Lineage_p) %>%
        separate_wider_delim(cols = Lineage, " ", names = c("Lineage", NA))

    mixed.frac.2.keep <- mixed.frac.2 %>%
        select(SampleID, sub_lineage) %>%
        distinct(SampleID, sub_lineage, .keep_all = TRUE) %>%
        separate_rows(sub_lineage, sep = ";") %>%
        select(SampleID, Lineage = sub_lineage)

    mixed.frac.3 <- left_join(mixed.frac.2.keep, mixed.frac.2) %>%
                        group_by(SampleID) %>%
                        mutate(Lineage_pL = paste(unique(Lineage_p), collapse = "; ")) %>%
                        slice(1) %>%  # Keeps only the first row per SampleID after mutating
                        ungroup() %>%
                        select(SampleID, Lineage_frac = Lineage_pL)

    mixed.frac.4 <- mixed.frac.1 %>% 
                        group_by(SampleID) %>%  # Group by SampleID
                        arrange(desc(Perc), .by_group = TRUE) %>%   # Sort by Perc in descending order within each group
                        slice(1) %>% # Retain only the first row of each group
                        filter(Perc >= 90) %>% # Filter rows where Perc >= 90
                        mutate(Mixed_90perc = paste0(Lineage, " (", Perc, "%)")) %>%  # Create new column
                        ungroup() # Ungroup after manipulation (optional but recommended)

    mixed.frac.5 <- left_join(mixed.frac.3, mixed.frac.4) %>%
                    select(SampleID, Lineage_frac, Mixed_90perc)

}

# Clonal: add fraction of lineage designation
clonal <- tbprof %>% filter(infection_type == "clonal")

if (nrow(clonal) == 0) { # Check if 'clonal' is empty

    clonal.frac.3 <- clonal

} else { # Proceed with processing if 'clonal' is not empty

    clonal.frac <- left_join(clonal, lineage.frac) %>%
                        select(SampleID, main_lineage, sub_lineage, 
                                infection_type, Lineage, Fraction)
    clonal.frac$Fraction <- as.numeric(clonal.frac$Fraction)
    clonal.frac <- clonal.frac %>% mutate(Perc = 100 * Fraction)

    clonal.frac.1 <- clonal.frac %>%
                        mutate(Lineage_p = paste0(Lineage, " (", Perc, "%)")) %>%
                        mutate(Lineage_p = gsub("lineage", "L", Lineage_p))

    clonal.frac.2 <- clonal.frac.1 %>%
                        select(SampleID, sub_lineage, Lineage_p) %>%
                        mutate(sub_lineage = gsub("lineage", "L", sub_lineage)) %>%
                        mutate(Lineage = Lineage_p) %>%
                        separate_wider_delim(cols = Lineage, " ", names = c("Lineage", NA))

    clonal.frac.2.keep <- clonal.frac.2 %>%
                        select(SampleID, sub_lineage) %>%
                        distinct(SampleID, sub_lineage, .keep_all = TRUE) %>%
                        separate_rows(sub_lineage, sep = ";") %>%
                        select(SampleID, Lineage = sub_lineage)

    clonal.frac.3 <- left_join(clonal.frac.2.keep, clonal.frac.2) %>%
                        group_by(SampleID) %>%
                        mutate(Lineage_pL = paste(unique(Lineage_p), collapse = "; ")) %>%
                        slice(1) %>%  # Keeps only the first row per SampleID after mutating
                        ungroup() %>%
                        select(SampleID, Lineage_frac = Lineage_pL) %>%
                        distinct() %>%
                        mutate(Mixed_90perc = NA)

}

# Combine clonal and clonal results
TBProfiler.lineages <- rbind(clonal.frac.3, mixed.frac.5) %>%
                        arrange(desc(SampleID))

TBProfiler.lineages.final <- left_join(TBProfiler.lineages, tbprof.infec) %>%
                                select(SampleID,infection_type,
                                        Lineage_frac,Mixed_90perc)

# Output the result
write.table(TBProfiler.lineages.final, file = OUT,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("Finished: R/tbprofiler_lineage_fractions.R\n")