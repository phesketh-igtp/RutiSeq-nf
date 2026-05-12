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

    # ------------------------------------------------------------------
    # Run for both datasets
    # ------------------------------------------------------------------
    stats_df = concat_clean(stats_files)
    class_df = concat_clean(class_files)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    stats_df.write_csv(
        "Mapping_and_Variant_Statistics.tab",
        separator="\t",
        include_header=False
    )

    class_df.write_csv(
        "Strain_Classification.tab",
        separator="\t",
        include_header=False
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
*/