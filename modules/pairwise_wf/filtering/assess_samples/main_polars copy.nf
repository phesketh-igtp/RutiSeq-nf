process ASSESS_SAMPLES {

    conda params.r_stats_env

    storeDir "${params.outDir}/db/comparison/src/${params.runID}/"

    input:
        path(pairwise_analysis_list)
        val(sampleID_list)

    output:
        path "final.lineage_samples_tuple.csv", emit: lineage_sample_tuple
        path "final.skipped-lineages_tuple.csv", emit: skipped_lineages_tuple

    script:
        def sub_lineages = params.lineage_pairwise_sub.join('\\n')
        def main_lineages = params.lineage_pairwise_main.join('\\n')

"""
import polars as pl

# ------------------------------------------------------------------
# Parameters (Nextflow injects these)
# ------------------------------------------------------------------
pairwise_split = "${params.pairwise_split}"
sampleID_string = "${sampleID_list}"
pairwise_analysis_file = "${pairwise_analysis_list}"

# ------------------------------------------------------------------
# Prepare inputs
# ------------------------------------------------------------------
sample_ids = sorted(set(s.strip() for s in sampleID_string.split(",")))

sub_lineages_list = sorted(set(
    s.strip() for s in "${sub_lineages}".split("\n") if s.strip()
))

main_lineages_list = sorted(set(
    s.strip() for s in "${main_lineages}".split("\n") if s.strip()
))

# Convert to DataFrames
run_ids_df = pl.DataFrame({"SampleID": sample_ids})
main_lineages_df = pl.DataFrame({"selected_main_lineage": main_lineages_list})
sub_lineages_df = pl.DataFrame({"selected_sub_lineage": sub_lineages_list})

# ------------------------------------------------------------------
# Read metadata
# ------------------------------------------------------------------
meta = (
    pl.read_csv(pairwise_analysis_file)
    .unique()
    .rename({
        "column_1": "SampleID",
        "column_2": "main_lineage",
        "column_3": "sub_lineage"
    })
    .filter(pl.col("main_lineage").is_not_null())
    .filter(~pl.col("main_lineage").str.contains(";"))
    .filter(~pl.col("SampleID").str.contains("CN-"))
)

# ------------------------------------------------------------------
# Analysis type logic
# ------------------------------------------------------------------
if pairwise_split == "sub":

    # Build regex for sub_lineage matching
    sub_patterns = "^(" + "|".join(sub_lineages_list) + r")\b"

    filtered_meta = (
        meta
        .with_columns(
            pl.when(pl.col("sub_lineage").str.contains(sub_patterns))
            .then(
                pl.col("sub_lineage")
                .str.extract(sub_patterns)
            )
            .otherwise(None)
            .alias("matched_sub")
        )
        .with_columns(
            pl.when(pl.col("matched_sub").is_not_null())
            .then(pl.col("matched_sub"))
            .when(pl.col("main_lineage").is_in(main_lineages_list))
            .then(pl.col("main_lineage"))
            .otherwise(pl.col("sub_lineage"))
            .alias("lineage")
        )
        .select(["lineage", "SampleID"])
    )

elif pairwise_split == "main":

    filtered_meta = meta.select([
        pl.col("main_lineage").alias("lineage"),
        "SampleID"
    ])

elif pairwise_split == "none":

    filtered_meta = meta.select([
        pl.lit("All").alias("lineage"),
        "SampleID"
    ])

else:
    raise ValueError(
        f"Invalid pairwise level: {pairwise_split}"
    )

# ------------------------------------------------------------------
# Filter small lineages (<5 samples)
# ------------------------------------------------------------------
counts = (
    filtered_meta
    .group_by("lineage")
    .count()
)

skipped = counts.filter(pl.col("count") < 5)

passed = (
    filtered_meta
    .filter(pl.col("SampleID").is_in(sample_ids))
    .filter(~pl.col("lineage").is_in(skipped["lineage"]))
    .select("lineage")
    .unique()
)

# ------------------------------------------------------------------
# Split forward vs skipped
# ------------------------------------------------------------------
filtered_meta_forward = (
    filtered_meta
    .filter(pl.col("lineage").is_in(passed["lineage"]))
    .filter(~pl.col("SampleID").is_in(["sample", "SampleID"]))
    .filter(~pl.col("lineage").is_in(["main_lineage", "lineage"]))
)

filtered_meta_skip = (
    filtered_meta
    .filter(~pl.col("lineage").is_in(passed["lineage"]))
    .filter(~pl.col("SampleID").is_in(["sample", "SampleID"]))
    .filter(pl.col("lineage") != "sub_lineage")
)

# ------------------------------------------------------------------
# Write outputs (no headers)
# ------------------------------------------------------------------
filtered_meta_forward.write_csv(
    "final.lineage_samples_tuple.csv",
    include_header=False
)

filtered_meta_skip.write_csv(
    "final.skipped-lineages_tuple.csv",
    include_header=False
)
"""
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 4.0.0
@description:
    In this module creates the pairwise analysis tuples from the lineage_samples_paths.csv
    and the lineage_pairwise_sub and lineage_pairwise_main lists.
    The output is a tuple of the form (lineage, sampleID) for each sampleID in the analysis.
    There are three options for the pairwise analysis (specified by the params.pairwise_split):
        - sub: pairwise analysis at sub-lineage level
        - main: pairwise analysis at main-lineage level
        - none: pairwise analysis of all samples without lineage split
    Now implemented entirely in R for better data handling and consistency.
@changelog:
    v4.0.1-2026-05-12: Convert from R to Polars (python)
    v3.0.0-2026-04-13: Complete rewrite in R, eliminating bash script components
    v2.0.0-2025-04-01: Updated to use the new lineage_pairwise_sub and lineage_pairwise_main lists
    v1.0.1-2024-11-01: Added error handling for invalid pairwise level
*/