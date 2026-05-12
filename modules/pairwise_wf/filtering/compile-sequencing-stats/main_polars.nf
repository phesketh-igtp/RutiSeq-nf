process COMPILE_SEQUENCING_STATS {

    conda params.r_stats_env

    publishDir "${params.outDir}/db/comparison/summary/", 
        mode: 'copy',
        overwrite: true

    input:
        tuple path(tbdb_results, stageAs: "tbdb-tbprofiler.txt"), 
            path(who_results, stageAs: "who-tbprofiler.txt"),
            path(lineage_fractions, stageAs: "lineages.fractions.txt")

        path(mtbseq_strains, stageAs: "Strain_Classification.tab")
        path(mtbseq_stats, stageAs: "Mapping_and_Variant_Statistics.tab")

    output:
        path("${params.runID}.sequencing_summary.csv")

        path("sequencing_summary.csv"),      emit: analysis_summary
        path("who_resistance_summary.csv"),  emit: who_resistance
        path("tbdb_resistance_summary.csv"), emit: tbdb_resistance

        path("pairwise_analysis.list.csv"),  emit: pairwise_analysis_list

    script:

"""
#!/usr/bin/env python
import polars as pl

# ------------------------------------------------------------------
# Params
# ------------------------------------------------------------------
dictionary = "${params.scriptDir}/R/dict/mtbseq_classification.dict.csv"
runID = "${params.runID}"

# ------------------------------------------------------------------
# Helper: dictionary rename
# ------------------------------------------------------------------
def dictionary_rename(df, dict_path):
    dict_df = pl.read_csv(dict_path)

    rename_map = dict(
        zip(dict_df["old.name"], dict_df["new.name"])
    )

    cols_to_keep = (
        dict_df
        .filter(pl.col("final") == "Y")
        ["new.name"]
        .to_list()
    )

    return (
        df.rename(rename_map)
          .select(cols_to_keep)
    )

# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
tbprof = (
    pl.read_csv("tbdb-tbprofiler.txt", separator="\t")
    .with_columns(pl.all().cast(pl.Utf8).str.strip_chars())
    .select([
        pl.col("sample").alias("SampleID"),
        "main_lineage",
        "sub_lineage"
    ])
    .with_columns(
        pl.when(
            pl.col("main_lineage").str.contains(";") |
            pl.col("sub_lineage").str.contains(";")
        )
        .then("mixed")
        .otherwise("clonal")
        .alias("infection_type")
    )
)

lineage_frac = (
    pl.read_csv("lineages.fractions.txt", separator="\t")
    .filter(pl.col("Fraction") != "Fraction")
    .with_columns(pl.col("Fraction").cast(pl.Float64))
    .with_columns((pl.col("Fraction") * 100).alias("Perc"))
)

# ------------------------------------------------------------------
# Function to process lineage fractions (replaces duplicate code)
# ------------------------------------------------------------------
def process_fraction(df):
    if df.height == 0:
        return df

    return (
        df.join(lineage_frac, on="SampleID", how="left")
        .with_columns([
            (pl.col("Lineage") + " (" + pl.col("Perc").round(2).cast(pl.Utf8) + "%)")
            .str.replace("lineage", "L")
            .alias("Lineage_p")
        ])
        .group_by("SampleID")
        .agg([
            pl.col("Lineage_p").unique().str.concat("; ").alias("Lineage_frac"),
            pl.when(pl.col("Perc") >= 90)
              .then(pl.col("Lineage_p"))
              .otherwise(None)
              .drop_nulls()
              .first()
              .alias("Mixed_90perc")
        ])
    )

# ------------------------------------------------------------------
# Split datasets
# ------------------------------------------------------------------
mixed   = tbprof.filter(pl.col("infection_type") == "mixed")
clonal  = tbprof.filter(pl.col("infection_type") == "clonal")

mixed_frac  = process_fraction(mixed)
clonal_frac = process_fraction(clonal).with_columns(
    pl.lit(None).alias("Mixed_90perc")
)

tbdb_lin_fract_final = (
    pl.concat([mixed_frac, clonal_frac])
    .join(tbprof.select(["SampleID", "infection_type"]),
          on="SampleID",
          how="left")
)

# ------------------------------------------------------------------
# Write lineage fraction output
# ------------------------------------------------------------------
tbdb_lin_fract_final.write_csv(
    "tbprofiler.lineages.fractions.txt",
    separator=";"
)

# ------------------------------------------------------------------
# Sequencing summary section
# ------------------------------------------------------------------
mtbseq_stats = pl.read_csv("Mapping_and_Variant_Statistics.tab", separator="\t", has_header=False)
mtbseq_class = pl.read_csv("Strain_Classification.tab", separator="\t", has_header=False)

tbprof_tbdb = pl.read_csv("tbdb-tbprofiler.txt", separator="\t")
tbprof_who  = pl.read_csv("who-tbprofiler.txt", separator="\t")

lineage_frac_short = tbdb_lin_fract_final.select([
    "SampleID", "Lineage_frac", "Mixed_90perc"
])

# Rename using dictionaries
mtbseq_statistics_df = dictionary_rename(
    mtbseq_stats,
    "${params.scriptDir}/R/dict/mtbseq_statistics.dict.csv"
)

mtbseq_classification_df = dictionary_rename(
    mtbseq_class,
    "${params.scriptDir}/R/dict/mtbseq_classification.dict.csv"
)

# Merge all
full_df = (
    mtbseq_statistics_df
    .join(mtbseq_classification_df, on="FullID", how="left")
    .join(tbprof_tbdb, left_on="FullID", right_on="sample", how="left")
    .join(lineage_frac_short, left_on="FullID", right_on="SampleID", how="left")
)

# Add infection type
full_df = full_df.with_columns(
    pl.when(
        pl.col("main_lineage").str.contains(";") |
        pl.col("sub_lineage").str.contains(";")
    )
    .then("Mixed")
    .otherwise("Clonal")
    .alias("infection_type")
)

# ------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------
sequencing_summary_df = dictionary_rename(
    full_df,
    "${params.scriptDir}/R/dict/sequencing-summary.dict.csv"
)

sequencing_summary_df.write_csv("sequencing_summary.csv")
sequencing_summary_df.write_csv(f"{runID}.sequencing_summary.csv")

dictionary_rename(
    tbprof_tbdb,
    "${params.scriptDir}/R/dict/resistance_profiles_TBDB.csv"
).write_csv("tbdb_resistance_summary.csv")

dictionary_rename(
    tbprof_who,
    "${params.scriptDir}/R/dict/resistance_profiles_WHO.csv"
).write_csv("who_resistance_summary.csv")

# ------------------------------------------------------------------
# Pairwise analysis
# ------------------------------------------------------------------
pairwise_analysis_df = (
    full_df
    .select([
        pl.col("FullID").alias("SampleID"),
        "main_lineage",
        "sub_lineage"
    ])
    .filter(
        (pl.col("main_lineage") != "NA") &
        ~pl.col("main_lineage").str.contains(";") &
        ~pl.col("sub_lineage").str.contains(";")
    )
)

pairwise_analysis_df.write_csv("pairwise_analysis.list.csv")
"""
}


/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.0
@description:
    In this module the sequencing statistics for the db are
    calculated with Rscripts. The output is a summary of the
    sequencing statistics, the tbdb and who resistance  summaries
    and a list of the genomes which pass the minimum quality
    requirements for the pairwise analysis. This is then used to create a 
    filtered list of genomes for the pairwise analysis, which is converated into
    a tuple/channel for downstream processing.
@changelog:
    v1.0.0_2024-04-01: Inital version of module
*/
