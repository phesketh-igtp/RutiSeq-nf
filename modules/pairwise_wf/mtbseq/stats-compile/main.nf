process MTBSEQ_STATS_COMPILE {

    conda params.stats_env
    //conda "conda-forge::polars-lts-cpu=0.20.31 conda-forge::python=3.11"

    publishDir "${params.outDir}/db/comparison/mtbseq/", 
        mode: 'copy',
        overwrite: true

    input:
        val(sampleID_list)

    output:
        path("Strain_Classification.tab"),              emit: mtbseq_compiled_strains
        path("Mapping_and_Variant_Statistics.tab"),     emit: mtbseq_compiled_map_stats

    script:

    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    # ------------------------------------------------------------------
    # Collect files
    # ------------------------------------------------------------------
    stats_files = list(Path("${params.outDir}/db/samples").rglob("*/mtbseq/Statistics/*.tab"))
    class_files = list(Path("${params.outDir}/db/samples").rglob("*/mtbseq/Classification/*.tab"))

    # ------------------------------------------------------------------
    # Function to process files
    # ------------------------------------------------------------------
    def concat_clean(files):
        return pl.concat([
            pl.read_csv(f, separator="\t", has_header=False)
            for f in files
        ]).filter(
            ~pl.col("column_1").str.starts_with("Date")
        ).with_columns(
            pl.all().cast(pl.Utf8).str.replace_all("'", "")
        )

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

    # ------------------------------------------------------------------
    # Run for both datasets
    # ------------------------------------------------------------------
    stats_df = concat_clean(stats_files)
    class_df = concat_clean(class_files)

    # ------------------------------------------------------------------
    # Rename using renaming maps
    # ------------------------------------------------------------------
    stats_df = stats_df.rename(mtbseq_stats_rename_map)
    class_df = class_df.rename(mtbseq_class_rename_map)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    stats_df.write_csv(
        "Mapping_and_Variant_Statistics.tab",
        separator="\t",
        include_header=TRUE
    )

    class_df.write_csv(
        "Strain_Classification.tab",
        separator="\t",
        include_header=TRUE
    )
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 2.0.0
@description:
    This process compiles the statistics and classification files from the MTBseq analysis.
    It concatenates the files from all samples and removes the header lines.
    The output files are saved in the specified output directory.
    The output files are:
        - Mapping_and_Variant_Statistics.tab
        - Strain_Classification.tab
@changelog:
    v1.0.0-2025-04-01: Initial version of script
    v2.0.0-2026-05-12: Conversion to polars for secuirty and speed
    v2.0.1-2026-05-12: Added headers to the outputs
*/