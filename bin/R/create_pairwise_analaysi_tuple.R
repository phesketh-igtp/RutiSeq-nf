#!/usr/bin/env R

library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(readr, quietly = TRUE)

# Import dataframes
main_lineages <- readr::read_delim("selected_main-lineage_split.list", col_names = FALSE, delim = ",")
colnames(main_lineages) <- "selected_main_lineage"

sub_lineages <- readr::read_delim("selected_sub-lineage_split.list", col_names = FALSE, delim = ",")
colnames(sub_lineages) <- "selected_sub_lineage"

run_ids <- readr::read_delim("run_sample_ids.txt", col_names = FALSE, delim = ",")
colnames(run_ids) <- "SampleID"

meta <- readr::read_delim("pairwise_analysis.list.csv", col_names = TRUE, delim = ",") |> 
    distinct() |>
    filter(!is.na(main_lineage) & !str_detect(main_lineage, ";")) |>
    filter(!str_detect(sample, "CN-"))  # This line filters out any sample that contains 'CN-'
colnames(meta) <- c("SampleID", "main_lineage", "sub_lineage")

# Filter out the lineages at sub_lineage level
filtered_meta <- meta %>%
    mutate(
            # Check if sub_lineage belongs to any selected_sub_lineage
            matched_sub = sapply(sub_lineage, function(sub) {
            match <- sub_lineages$selected_sub_lineage[str_detect(sub, paste0("^", sub_lineages$selected_sub_lineage, "\\b"))]
            if (length(match) > 0) match[1] else NA
            }),
            
            # Determine final lineage
            final_lineage = case_when(
            !is.na(matched_sub) ~ matched_sub,  # If sub_lineage matches, use it
            main_lineage %in% main_lineages$selected_main_lineage ~ main_lineage,  # Otherwise, use main_lineage if valid
            TRUE ~ sub_lineage  # Default to sub_lineage if no match
            )
        ) |> 
    select(final_lineage, SampleID) |> 
    rename(lineage = final_lineage)

# Step 2: Get unique lineages from run_ids
filtered_lineages_forward <- filtered_meta |> 
    filter(SampleID %in% run_ids$SampleID) |> 
    count(lineage) |> 
    filter(n > 4) |> # rm lineage with less than 4 genome (min for MTBSeq)
    select(lineage) |> 
    distinct()  # Use distinct() instead of unique() for dplyr consistency


# Step 3: Partition filtered_meta based on filtered_lineages_forward
filtered_meta_forward <- filtered_meta %>%
    filter(lineage %in% filtered_lineages_forward$lineage)

filtered_meta_skip <- filtered_meta %>%
    filter(!lineage %in% filtered_lineages_forward$lineage)

# Step 4: export the csv
write.csv(filtered_meta_forward, "final.lineage_samples_tuple.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.csv(filtered_meta_skip, "final.skipped-lineages_tuple.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
