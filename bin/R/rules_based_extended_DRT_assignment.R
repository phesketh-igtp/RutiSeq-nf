#!/bin/bash/R

#' @tbprofiler_extended_drt_classified.R
#' @author Poppy J hesketh Best, Antoni Bordoy
#' @description
#' Script that takes information from `tb-profiler profiler --db who` output
#' and add a new drug-resistance type classification based on the CDC TB
#' drug resistance classifications.
#' @workflow
#' 1. Imports the WHO TB-Profiler outouts from `tb-profiler collate`
#' 2. import the tbprofiler who rules (CSV)
#' 3. Creates a list of the antibiotics to be inspected and replaces '-' with
#'        'S' to denote that no mutation was found, thus the genome is predicted
#'        to have the sensistive phenotpy
#' 4. Add the function that applies the rules
#' @describeIn -i, --input TB-profiler collated output
#' @describeIn -r, --rules Rules CSV of assignment of extended drug resistance
#' @version 1.0.0
#' @date 2025-05-13
#' @chagelog
#'    v1.0.0-2025+05-13: Intial versions (draft)

################################################################################

library(tidyverse)
set.seed(1234)

################################################################################

# Import the rules
rule_df <- readr::read_delim("tbprofiler_drt_rules.csv",
                             quote = "\"",
                             col_types = cols(.default = "c"),
                             comment = "#")

# Create a categories list
categories <- rule_df |> select(DRT = DRT_type, extended_DRT = DRT)

# Import the tb-profiler collated results
df <- read.delim("who_resistance_summary.csv",
                 sep = ",") |>
  select(sample,
         RIF, INH, EMB, PZA,
         MFX, LFX, BDQ, DLM,
         Pa, LZD, STM, AMK, KAN,
         CAP, CFZ, ETO, PAC,
         CYR)

# Columns you want to update
cols_to_update <- c(
  "RIF", "INH", "EMB", "PZA",
  "MFX", "LFX", "BDQ", "DLM",
  "Pa", "LZD", "STM", "AMK",
  "KAN", "CAP", "CFZ", "ETO",
  "PAC", "CYR"
)

# Replace "-" with "S" in those columns
df[cols_to_update] <- lapply(df[cols_to_update],
                             function(x) {
                                          ifelse(is.na(x), NA,
                                                 ifelse(x == "-", "S", "R"))})

rm(cols_to_update)

# Function to parse each row
parse_rule <- function(row) {

  rule_name <- row[["DRT"]]
  required <- list()
  optional_raw <- list()

  for (drug in names(row)[-1]) {
    val <- row[[drug]]
    if (is.na(val) || val == "") next
    if (val %in% c('== "S"', '!= "S"')) {
      required[[drug]] <- val
    } else if (grepl("!=", val) && grepl("==", val)) {
      optional_raw[[drug]] <- TRUE
    }
  }

  # Group optional drugs by category if needed
  optional_group <- list()

  if (length(optional_raw) > 0) {
    # Custom logic: group known types like FQs and injectables
    FQs <- intersect(names(optional_raw), c("MFX", "LFX"))
    INJ <- intersect(names(optional_raw), c("AMK", "KAN", "CAP"))
    OTHER <- setdiff(names(optional_raw), c(FQs, INJ))
    
    if (length(FQs) > 0) optional_group <- append(optional_group, list(FQs))
    if (length(INJ) > 0) optional_group <- append(optional_group, list(INJ))
    if (length(OTHER) > 0) optional_group <- append(optional_group, list(OTHER))
  }
  
  list(
    name = rule_name,
    required = required,
    optional_groups = optional_group,
    description = rule_name
  )
}

# Apply parsing across all rules
rules_list <- apply(rule_df, 1, parse_rule)

################################################################################

apply_rules <- function(samples_df, rules_list) {
  results <- vector("character", nrow(samples_df))

  for (i in seq_len(nrow(samples_df))) {
    sample <- samples_df[i, ]

    matched_rule <- NA

    for (rule in rules_list) {
      # Check required conditions
      required <- rule$required
      passed_required <- TRUE
      for (drug in names(required)) {
        if (!drug %in% colnames(samples_df)) {
          passed_required <- FALSE
          break
        }
        expr <- required[[drug]]
        value <- sample[[drug]]
        condition <- eval(parse(text = paste0('"', value, '" ', expr)))
        if (!condition) {
          passed_required <- FALSE
          break
        }
      }

      if (!passed_required) next

      # Check optional groups (at least one drug in each group must be not "S")
      optional_pass <- TRUE
      for (group in rule$optional_groups) {
        values <- sample[intersect(group, colnames(sample))]
        if (all(values == "S" | is.na(values))) {
          optional_pass <- FALSE
          break
        }
      }

      if (optional_pass) {
        matched_rule <- rule$name
        break  # first match wins
      }
    }

    results[i] <- ifelse(is.na(matched_rule), "Other", matched_rule)
  }

  return(results)
}

results <- apply_rules(df, rules_list)

df$extended_DRT <- results
colnames(df)

final_df <- left_join(df, categories, relationship = "many-to-many")

final_df <- final_df |> select(sample,DRT,DRT_ext=extended_DRT)

################################################################################

# Export
write.csv2(final_df, "tbprofiler-DR-corrections.csv",
           row.names = FALSE,
           quote = TRUE,
           na = "N/A",
           fileEncoding = "UTF-8")
