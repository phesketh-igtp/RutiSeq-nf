process COMPILE_SEQUENCING_STATS {

    conda params.stats_env
    //conda "conda-forge::polars-lts-cpu=0.20.31 conda-forge::python=3.11"

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
    runID = "${params.runID}"

    # ------------------------------------------------------------------
    # Helper: dictionary rename//renaming maps
    # ------------------------------------------------------------------
    rename_map_final = {
        "FullID": "Sample",
        "LibraryID": "Library",
        "main_lineage": "Main lineage",
        "sub_lineage": "Sub-lineage",
        "Total Reads": "Total Reads",
        "Unambiguous Total Bases": "Unambiguous Total Bases",
        "Unambiguous Total Bases (%)": "Unambiguous Total Bases (%)",
        "Unambiguous GC-Content": "Unambiguous GC Content",
        "Unambiguous Coverage median": "Unambiguous Coverage median",
        "SNPs": "nSNPs",
        "drtype": "Drug resistance type",
        "infection_type": "Infection type",
        "lineages_frac": "Lineages (fractions)",
        "status": "Status"
    }

    # Define the fucntion for breaking down the lineages
    def lineage_level(expr: pl.Expr, n: int) -> pl.Expr:
        parts = expr.cast(pl.Utf8).str.strip_chars().str.split(".")
        return (
            pl.when(parts.list.len() >= n)
            .then(parts.list.slice(0, n).list.join("."))
            .otherwise(None)
        )

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    tbprof = (
        pl.read_csv(
            "tbdb-tbprofiler.txt",
            separator="\t")
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
        .then(pl.lit("mixed"))
        .otherwise(pl.lit("clonal"))
        .alias("infection_type")
        )
    )

    # ------------------------------------------------------------------
    # Create the lineage fractions column
    # ------------------------------------------------------------------

    lineage_frac = (
        pl.read_csv("lineages.fractions.txt", separator="\t")
        .filter(pl.col("Fraction") != "Fraction")
        .with_columns(pl.col("Fraction").cast(pl.Float64))
        .with_columns((pl.col("Fraction") * 100).alias("Perc"))
    )

    lineage_frac_collapsed = (
        lineage_frac
        .with_columns(
            pl.col("Lineage").str.replace(r"^lineage", "L").alias("Lineage")
            # If you want to be extra safe about case:
            # pl.col("Lineage").str.replace(r"(?i)^lineage", "L").alias("Lineage")
        )
        .sort(["SampleID", "Fraction", "Lineage"], descending=[False, True, False])
        .group_by("SampleID")
        .agg(
            pl.concat_str(
                [
                    pl.col("Lineage"),
                    pl.lit(" ("),
                    pl.col("Fraction").round(3).cast(pl.Utf8),
                    pl.lit(")")
                ],
                separator=""
            ).str.join(", ").alias("lineages_frac")
        )
    )

    # Join the final tbprofiler file
    tbprof_final = tbprof.join(
        lineage_frac_collapsed,
        on="SampleID"
    )

    # ------------------------------------------------------------------
    # Write lineage fraction output
    # ------------------------------------------------------------------
    tbprof_final.write_csv(
        "tbprofiler.lineages.fractions.txt",
        separator=";"
    )

    # ------------------------------------------------------------------
    # Sequencing summary section
    # ------------------------------------------------------------------
    mtbseq_stats = pl.read_csv(
        "Mapping_and_Variant_Statistics.tab", 
        separator="\t", 
        has_header=True)
    mtbseq_class = pl.read_csv(
        "Strain_Classification.tab",
        separator="\t",
        has_header=True)
    tbprof_tbdb = pl.read_csv(
        "tbdb-tbprofiler.txt",
        separator="\t")
    tbprof_who  = pl.read_csv(
        "who-tbprofiler.txt",
        separator="\t")

    # Merge all
    full_df = (
        mtbseq_stats
        .join(mtbseq_class, on="FullID", how="left")
        .join(tbprof_tbdb, left_on="FullID", right_on="sample", how="left")
        .join(tbprof_final, left_on="FullID", right_on="SampleID", how="left")
    )

    # ------------------------------------------------------------------
    # Sequencing outcome
    # ------------------------------------------------------------------

    full_df_with_status = (
        full_df
        .with_columns(
            # Build a list of failed checks (nulls are removed later)
            pl.concat_list([
                pl.when(pl.col("infection_type") == "clonal")
                .then(None)
                .otherwise(pl.lit("mixed infection")),

                pl.when(pl.col("Total Reads") >= 1_000_000)
                .then(None)
                .otherwise(pl.lit("Total Reads < 1,000,000")),

                pl.when(pl.col("Unambiguous Total Bases (%)") >= 0.95)
                .then(None)
                .otherwise(pl.lit("Unambiguous Total Bases (%) < 0.95")),

                pl.when(pl.col("Unambiguous Coverage median") >= 50)
                .then(None)
                .otherwise(pl.lit("Unambiguous Coverage median < 50")),
            ])
            .list.drop_nulls()
            .alias("fail_reasons")
        )
        .with_columns(
            # Convert fail_reasons -> status string
            pl.when(pl.col("fail_reasons").list.len() == 0)
            .then(pl.lit("PASS"))
            .otherwise(
                pl.concat_str(
                    [
                        pl.lit("FAIL :"),
                        pl.col("fail_reasons").list.join(", "),
                        pl.lit(")")
                    ],
                    separator=""
                )
            )
            .alias("status")
        )
        # If you don't want to keep fail_reasons, drop it:
        .drop("fail_reasons")
        # then select what you want in the final output:
        .select([
            pl.col("FullID"),
            pl.col("main_lineage"),
            pl.col("sub_lineage"),
            pl.col("Total Reads"),
            pl.col("infection_type"),
            pl.col("Unambiguous Total Bases (%)"),
            pl.col("Unambiguous Coverage median"),
            pl.col("status")
        ])
    )

    ##full_df_with_status.filter(pl.col("status") == "FAIL")

    # one final merge with the outputs
    # Merge all
    full_df = full_df.join(
        full_df_with_status.select(["FullID", "status"]),  # keep only what you need
        on="FullID",
        how="left"
    )

    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------

    # filter the final df for the summary page
    full_df = (
        full_df
        .rename(rename_map_final)
        .select(list(rename_map_final.values()))
    )
    sequencing_summary_df.write_csv("sequencing_summary.csv")
    sequencing_summary_df.write_csv(f"{runID}.sequencing_summary.csv")
    tbprof_tbdb.write_csv("tbdb_resistance_summary.csv")
    tbprof_who.write_csv("who_resistance_summary.csv")

    # ------------------------------------------------------------------
    # Pairwise analysis (genome that pass minimum quality)
    # ------------------------------------------------------------------
        
    # Read the samplesheet to identify the new samples
    new_sample_list = (
        pl.read_csv(f"{params.outDir}/sample-sheets/{params.runID}.csv")
        .select(
            pl.col("sampleID")
            .cast(pl.Utf8)
            .str.strip_chars()
            .alias("sampleID")
        )
        .filter(
            pl.col("sampleID").is_not_null() &
            (pl.col("sampleID") != "") &
            (pl.col("sampleID").str.to_lowercase() != "sampleid")
        )
        .unique().to_series().to_list()
    )

    # Filter to keep only samples in the list
    passed_new = (
        full_df
        .filter(
            (pl.col("Status") == "PASS") &
            (pl.col("Sample").is_in(new_sample_list))
        )
        .sort("Main lineage")
        .unique()
    )
    # Break down the lineages levels 
    lineage_lvs = (
        passed_new
        .with_columns([
            lineage_level(pl.col("Sub-lineage"), 1).alias("Main_lineage_lv1"),
            lineage_level(pl.col("Sub-lineage"), 2).alias("Sub_lineage_lv2"),
            lineage_level(pl.col("Sub-lineage"), 3).alias("Sub_lineage_lv3"),
            lineage_level(pl.col("Sub-lineage"), 4).alias("Sub_lineage_lv4"),
        ])
    )
    # Identify the lineages they belongue to
    lv1 = lineage_lvs.get_column("Main_lineage_lv1").drop_nulls().unique().to_list()
    lv2 = lineage_lvs.get_column("Sub_lineage_lv2").drop_nulls().unique().to_list()
    lv3 = lineage_lvs.get_column("Sub_lineage_lv3").drop_nulls().unique().to_list()
    lv4 = lineage_lvs.get_column("Sub_lineage_lv4").drop_nulls().unique().to_list()

    # Collect all the pass genomes that have the correct (sub-)lineages
    pairwise_analysis_genomes = (
        full_df
        .filter(
            (pl.col("Status") == "PASS") &
            (
                pl.col("Main lineage").is_in(lv1) |
                pl.col("Sub-lineage").is_in(lv2) |
                pl.col("Sub-lineage").is_in(lv3) |
                pl.col("Sub-lineage").is_in(lv4)
            )
        )
        .select(["Sample", "Main lineage", "Sub-lineage"])
        .rename({
            "Sample": "sample",
            "Main lineage": "main_lineage",
            "Sub-lineage": "sub_lineage",
        })
    )

    # Export the minimum quality genomes
    pairwise_analysis_genomes.write_csv(
        "pairwise_analysis.list.csv",
        include_header=False
    )
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2026-05-13
@version: 2.0.0
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
    v2.0.0_2026-05-13: Migrated R code to Python (polars)
*/
