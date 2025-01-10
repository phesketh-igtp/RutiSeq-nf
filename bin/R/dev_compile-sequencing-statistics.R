# load libraries
library(tidyverse)

# Import dataframes
mtbseq_stats       <- read.delim("Mapping_and_Variant_Statistics.tab", header = FALSE)
mtbseq_class       <- read.delim("Strain_Classification.tab", header = FALSE)
tbprofiler_tbdb.df <- read.delim("tbdb-tbprofiler.txt", header = TRUE)
tbprofiler_who     <- read.delim("who-tbprofiler.txt", header = TRUE)
minimum_coverage   <- 50
runID              <- 'test-03'


# Add appropriate headers to dataframes using dictionaries

# import local functions
source("/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R/functions/dictionary-rename.R")

dictionary_path <-  "/imppc/labs/emlab/share/GitHub/RutiSeq-nf/bin/R/"

# Rename the dfs according to dictionary
mtbseq_statistics.df     <- dictionary_rename(df = mtbseq_stats,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_statistics.dict.csv"))

mtbseq_classification.df <- dictionary_rename(df = mtbseq_class,
                            dict_path = paste0(dictionary_path,
                            "/dict/mtbseq_classification.dict.csv"))

merge1 <- left_join(mtbseq_statistics.df, mtbseq_classification.df)
merge2 <- left_join(merge1, tbprofiler_tbdb.df, by = c("FullID" = "sample"))
full.df <- merge2 # rename the dataframe

# Filter out the genomes using a minimum coverage 
full.df.final <- full.df |>
            mutate(infection_type = if_else(
                grepl(";", main_lineage) | grepl(";", sub_lineage),
                        "Mixed","Clonal")
                )

# Creating the outputs for later assembly into final results
sequencing_summary.df <- dictionary_rename(df = full.df.final,
                            dict_path = paste0(dictionary_path,
                            "/dict/sequencing-summary.dict.csv"))

resistance_profiles_TBDB.df <- dictionary_rename(df = tbprofiler_tbdb.df,
                                dict_path = paste0(dictionary_path,
                                "/dict/resistance_profiles_TBDB.csv"))

resistance_profiles_WHO.df <- dictionary_rename(df = tbprofiler_who,
                        dict_path = paste0(dictionary_path,
                        "/dict/resistance_profiles_WHO.csv"))

# Create the list of genomes for pairwise analysis
pairwise_analysis.df <- full.df.final |>
                filter(`(Unambiguous) Coverage median` >= minimum_coverage) |>
                filter(infection_type == "Clonal") |> #only clonal genomes
                filter(sub_lineage != "NA") |> # no unclassified genomes
                select(SampleID=FullID,sub_lineage)

# export all the ouputs (broken!)
write.csv2(sequencing_summary.df,
            "sequencing_summary.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(sequencing_summary.df,
            paste(runID, ".sequencing_summary.csv"),
            quote=FALSE, row.names=FALSE)

write.csv2(resistance_profiles_TBDB.df,
            "resistance_profiles_TBDB.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(resistance_profiles_WHO.df,
            "resistance_profiles_WHO.csv",
            quote=FALSE, row.names=FALSE)

write.csv2(pairwise_analysis.df,
            "pairwise_analysis.list.csv",
            quote=FALSE, row.names=FALSE)
