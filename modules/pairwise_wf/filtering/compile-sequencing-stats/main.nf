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
    mtbseq_stats_rename_map = {
        "column_1": "Date",
        "column_2": "SampleID",
        "column_3": "LibraryID",
        "column_4": "FullID",
        "column_5": "Total Reads",
        "column_6": "Mapped Reads",
        "column_7": "Mapped Reads (%)",
        "column_8": "Genome Size",
        "column_9": "Genome GC",
        "column_10": "Any Total Bases",
        "column_11": "Any Total Bases (%)",
        "column_12": "Any GC-Content",
        "column_13": "Any Coverage mean",
        "column_14": "Any Coverage median",
        "column_15": "Unambiguous Total Bases",
        "column_16": "Unambiguous Total Bases (%)",
        "column_17": "Unambiguous GC-Content",
        "column_18": "Unambiguous Coverage mean",
        "column_19": "Unambiguous Coverage median",
        "column_20": "SNPs",
        "column_21": "Deletions",
        "column_22": "Insertions",
        "column_23": "Uncovered",
        "column_24": "Substitutions (Including Stop Codons)"
    }
    
    mtbseq_class_rename_map = {
        "column_1": "Date",
        "column_2": "SampleID",
        "column_3": "LibraryID",
        "column_4": "FullID",
        "column_5": "Homolka species",
        "column_6": "Homolka lineage",
        "column_7": "Homolka group",
        "column_8": "Quality",
        "column_9": "Coll lineage (branch)",
        "column_10": "Coll lineage_name (branch)",
        "column_11": "Coll quality (branch)",
        "column_12": "Coll lineage (easy)",
        "column_13": "Coll lineage_name (easy)",
        "column_14": "Coll quality (easy)",
        "column_15": "Beijing lineage (easy)",
        "column_16": "Beijing quality (easy)"
    }
    
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
        "Lineage_frac": "Lineages (fractions)",
        "Mixed_90perc": "Lineages (mixed fractions > 90)"
    }
    
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
            .then(pl.lit("mixed"))
            .otherwise(pl.lit("clonal"))
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
            df.join(lineage_frac, on="SampleID", how="left", coalesce=True)
            .with_columns([
                (pl.col("Lineage") + " (" + pl.col("Perc").round(2).cast(pl.Utf8) + "%)")
                .str.replace("lineage", "L")
                .alias("Lineage_p")
            ])
            .group_by("SampleID")
            .agg([
                pl.concat_str(pl.col("Lineage_p").unique(), separator="; ").alias("Lineage_frac"),
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
        clonal_frac
        .join(
            tbprof.select(["SampleID", "infection_type"]),
            on="SampleID",
            how="left"
        )
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
    
    # Rename using renaming maps
    mtbseq_stats = mtbseq_stats.rename(mtbseq_stats_rename_map)
    mtbseq_class = mtbseq_class.rename(mtbseq_class_rename_map)
    
    # Merge all
    full_df = (
        mtbseq_stats
        .join(mtbseq_class, on="FullID", how="left")
        .join(tbprof_tbdb, left_on="FullID", right_on="sample", how="left")
        .join(lineage_frac_short, left_on="FullID", right_on="SampleID", how="left")
    )
    
    # Add infection type
    full_df = full_df.with_columns(
        pl.when(
            pl.col("main_lineage").str.contains(";") |
            pl.col("sub_lineage").str.contains(";")
        )
        .then(pl.lit("Mixed"))
        .otherwise(pl.lit("Clonal"))
        .alias("infection_type")
    )
    
    # ------------------------------------------------------------------
    # Outputs
    # ------------------------------------------------------------------
    
    # filter the final df for the summary page
    sequencing_summary_df = (
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
    
    # Read the samplesheet
    samplelist = (
        pl.read_csv("${params.outDir}/sample-sheets/${params.runID}.csv")
        .select("sampleID")
        .to_series()
        .to_list()
    )
    
    minQual_genomes = (
        full_df
        # 1. keep only clonal samples
        .filter(pl.col("infection_type") == "Clonal")
        # 2. apply quality filters
        .filter(
            (pl.col("Total Reads") >= 1000000) &
            (pl.col("Unambiguous Total Bases (%)") >= 0.95) &
            (pl.col("Unambiguous Coverage median") >= 50)
        )
        # 3. select required columns
        .select([
            pl.col("FullID").alias("sample"), # or FullID depending on your naming
            pl.col("main_lineage"),
            pl.col("sub_lineage")
        ])
    )
    
    # Export the minimum quality genomes
    minQual_genomes.write_csv(
        "pairwise_analysis.list.csv",
        include_header=False
    )
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
