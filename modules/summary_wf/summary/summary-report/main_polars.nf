process GENERATE_SUMMARY_REPORT {

    conda params.r_stats_env

    publishDir "${params.outDir}/results/${params.runID}/", 
        mode: 'copy', 
        overwrite: true

    input:
        path pairwise_clusters_processed

        path sequencing_summary

        path who_resistance
        path tbdb_resistance

        tuple path(sylph_sequence_abundance, stageAs: 'sylph_results.tsv'),
            path(sylph_relative_abundance, stageAs: 'sylph_relative_abundance.tsv'),
            path(sylph_coverage, stageAs: 'sylph_coverage.tsv')

        tuple path(html_report), 
            path(warnings, stageAs: 'warnings.out')
        

    output:
        path("${params.runID}_RutiSeq-results.xlsx")
        path(html_report)

    script:

    """
    import polars as pl

    # -----------------------------
    # Parameters (replace as needed)
    # -----------------------------
    filt_min_cov = ${params.filt_min_cov}
    filt_min_depth = ${params.filt_min_depth}
    filt_min_reads = ${params.filt_min_reads}
    version = "${params.version}"
    runID = "${params.runID}"

    # -----------------------------
    # Helper: dictionary rename
    # -----------------------------

    #TODO

    # -----------------------------
    # Load data
    # -----------------------------
    summary = pl.read_csv("sequencing_summary.csv").unique()
    who_res = pl.read_csv("who_resistance_summary.csv").unique()
    tbdb_res = pl.read_csv("tbdb_resistance_summary.csv").unique()
    clusters = pl.read_csv("processed_clusters.tsv", separator="\t").unique()
    sylph_out = pl.read_csv("sylph_results.tsv", separator="\t")

    # -----------------------------
    # Resistance tables
    # -----------------------------
    who_res_min = who_res.select([
        pl.col("sample").alias("Sample"),
        pl.col("DR type").alias("DRtype (WHO)"),
        "INH","EMB","PZA","MFX","LFX","BDQ","DLM","Pa","LZD","STM",
        "AMK","KAN","CAP","CFZ","ETO","PAC","CYR"
    ])

    tbdb_res_min = tbdb_res.select([
        pl.col("sample").alias("Sample"),
        pl.col("DR type").alias("DRtype (TBDB)")
    ])

    resistance_summary = who_res_min.join(tbdb_res_min, on="Sample", how="left")

    # -----------------------------
    # Summary processing
    # -----------------------------
    summary = summary.with_columns(
        pl.col("Sample").str.replace(r"_.*", "").alias("OriginalID")
    )
    summary_tmp = summary.join(resistance_summary, on="Sample", how="left").unique()
    summary_xlsx = dictionary_rename(summary_tmp, dict_path)

    # -----------------------------
    # Status column (case_when equivalent)
    # -----------------------------
    summary_xlsx = summary_xlsx.with_columns(
        pl.when(
            (pl.col("Unambiguous Total Bases (%)") < filt_min_cov) &
            (pl.col("Unambiguous Coverage median") < filt_min_depth)
        ).then(
            f"FAILED: Percentage coverage <{filt_min_cov}, and coverage <{filt_min_depth}x"
        ).when(
            pl.col("Unambiguous Total Bases (%)") < filt_min_cov
        ).then(
            f"FAILED: Percentage coverage <{filt_min_cov}"
        ).when(
            pl.col("Unambiguous Coverage median") < filt_min_depth
        ).then(
            f"FAILED: Coverage <{filt_min_depth}X"
        ).when(
            pl.col("Total Reads") < filt_min_reads
        ).then(
            f"FAILED: <{filt_min_reads} reads"
        ).when(
            pl.col("Infection type") == "Mixed"
        ).then(
            "FAILED: mixed (sub-)lineage"
        ).otherwise("PASS").alias("Status")
    )

    summary_xlsx = summary_xlsx.with_columns([
        pl.lit(version).alias("Version Control"),
        pl.lit(runID).alias("AnalysisID")
    ])

    # -----------------------------
    # Write Excel workbook
    # -----------------------------
    # Polars doesn't natively do multi-sheet Excel writing,
    # so convert to pandas temporarily
    import pandas as pd

    with pd.ExcelWriter(f"{runID}_RutiSeq-results.xlsx", engine="openpyxl") as writer:
        summary_xlsx.to_pandas().to_excel(writer, sheet_name="sequencing_summary", index=False)
        who_res.to_pandas().to_excel(writer, sheet_name="resistance_table_who", index=False)
        tbdb_res.to_pandas().to_excel(writer, sheet_name="resistance_table_tbdb", index=False)
        clusters.to_pandas().to_excel(writer, sheet_name="transmission_clusters", index=False)
        sylph_out.to_pandas().to_excel(writer, sheet_name="sylph_classification", index=False)
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.1
@description:
    Generates a summary report of the analysis results, including
    a summary of the clusters, phylogeny, and variant sites.
    The report is generated in XLSX format and includes tables.
*/