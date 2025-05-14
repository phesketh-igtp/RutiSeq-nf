classify_tbprofiler_extended_drt <- function(input_file,
                                             rules_file,
                                             output_file) {
  library(tidyverse)

  # Import rules
  rule_df <- readr::read_delim(rules_file,
                               quote = "\"",
                               col_types = cols(.default = "c"),
                            comment = "#")

  # Extract categories
  categories <- rule_df |>
      select(DRT = DRT_type, extended_DRT = DRT)

  # Import TB-Profiler collated results
  df <- read.delim(input_file) |> 
    select(sample, rifampicin, isoniazid, ethambutol, pyrazinamide,
          moxifloxacin, levofloxacin, bedaquiline, delamanid, 
          pretomanid, linezolid, streptomycin, amikacin, kanamycin, 
          capreomycin, clofazimine, ethionamide, 
          para.aminosalicylic_acid, cycloserine)
  
  # Replace "-" with "S" in relevant columns
  cols_to_update <- c("rifampicin", "isoniazid", "ethambutol", "pyrazinamide", 
                      "moxifloxacin", "levofloxacin", "bedaquiline", "delamanid", 
                      "pretomanid", "linezolid", "streptomycin", "amikacin", 
                      "kanamycin", "capreomycin", "clofazimine", "ethionamide", 
                      "para.aminosalicylic_acid", "cycloserine")
  
  df[cols_to_update] <- lapply(df[cols_to_update], function(x) {
    ifelse(is.na(x), NA, ifelse(x == "-", "S", "R"))
  })
  
  # Rule parsing function
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
    
    optional_group <- list()
    if (length(optional_raw) > 0) {
      FQs <- intersect(names(optional_raw), c("moxifloxacin", "levofloxacin"))
      INJ <- intersect(names(optional_raw), c("amikacin", "kanamycin", "capreomycin"))
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
  
  # Apply rule parsing
  rules_list <- apply(rule_df, 1, parse_rule)
  
  # Rule application function
  apply_rules <- function(samples_df, rules_list) {
    results <- vector("character", nrow(samples_df))
    
    for (i in seq_len(nrow(samples_df))) {
      sample <- samples_df[i, ]
      matched_rule <- NA
      
      for (rule in rules_list) {
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
          break
        }
      }
      
      results[i] <- ifelse(is.na(matched_rule), "Other", matched_rule)
    }
    
    return(results)
  }
  
  # Apply rules
  df$extended_DRT <- apply_rules(df, rules_list)
  
  # Join categories
  final.df <- left_join(df, categories, relationship = "many-to-many") |>
    select(sample, DRT = extended_DRT, DRT_ext = DRT)
  
  # Export result
  write.csv2(final.df, output_file, row.names = FALSE)
}
