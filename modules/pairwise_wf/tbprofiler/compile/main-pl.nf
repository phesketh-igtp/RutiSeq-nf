process TBPROFILER_COMPILE {

    conda params.tbprofiler_env

    container { 
        if (workflow.containerEngine == 'singularity') {
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cbf8de71c4b6e9b044bbbf6ef573ab58e14bf75a846c7bc84dfbe03ac0e278c1/data'
        } else { 
            'quay.io/biocontainers/tb-profiler' 
        }
    }

    publishDir "${params.outDir}/db/comparison/tbprofiler/", 
        mode: 'copy',
        overwrite: true

    input:
        val(sampleID_list)
        //val(db_compliance_check)
        //path(tbprofiler_results)

    output:
    tuple path("tbdb-tbprofiler.txt"), 
        path("who-tbprofiler.txt"),     
        path("lineages.fractions.txt"), emit: tbdb_out

    script:
    """
    #!/usr/bin/env python
    import polars as pl
    from pathlib import Path

    # ------------------------------------------------------------------
    # Collect files
    # ------------------------------------------------------------------
    tbprofiler_files = list(Path("${params.outDir}/db/samples").rglob("*/tbprofiler/tbprofiler.csv"))
    lineage_files = list(Path("${params.outDir}/db/samples").rglob("*/tbprofiler/tbdb-*.txt"))

    # ------------------------------------------------------------------
    # Function to process files
    # ------------------------------------------------------------------

    # Concatenate the python files
    def concat_clean(files, skip_header=False):
        dfs = [
            pl.read_csv(
                f,
                separator=",",
                has_header=True,  # use headers now
                skip_rows=0       # no need if has_header=True
            ).with_columns(
                pl.all().cast(pl.Utf8).str.replace_all("'", "")
            )
            for f in files
        ]

        return pl.concat(dfs, how="diagonal")

    # ------------------------------------------------------------------
    # Collate TBProfiler summaries
    # ------------------------------------------------------------------

    # collate
    tbprofiler_df = concat_clean(tbprofiler_files)

    # Split into TBDB and WHO results
    tbprofiler_df = tbprofiler_df.with_columns([
        # Create database column based on prefix
        pl.when(pl.col("sample").str.starts_with("tbdb-"))
        .then(pl.lit("TBDB"))
        .when(pl.col("sample").str.starts_with("who-"))
        .then(pl.lit("WHO"))
        .otherwise(None)
        .alias("database"),

        # Remove prefixes from sample names
        pl.col("sample")
        .str.replace("^tbdb-|^who-", "")
        .alias("sample")
    ])

    
tbprofiler_tbdb = tbprofiler_df.filter(pl.col("database") == "TBDB")
tbprofiler_who  = tbprofiler_df.filter(pl.col("database") ==
 = tbprofiler_df.filter(pl.col("database") == "TBDB")
    tbprofiler_who  = tbprofiler_df.filter(pl.col("database") == "WHO")

    # ------------------------------------------------------------------
    # Collect all the lineage information
    # ------------------------------------------------------------------

    dfs = []

    for f in lineage_files:
        # Extract ID (mimics basename without suffix)
        sample_id = f.stem.replace(".results", "")

        # Read as single-column raw text
        df = pl.read_csv(
            f,
            separator="\n",
            has_header=False,
            new_columns=["line"]
        )

        # Apply equivalent of sed filters
        df = df.with_columns(
            pl.col("line").str.replace_all("-", "").alias("line")
        ).filter(
            # remove everything before "Lineage report"
            pl.col("line").cum_count().over(
                (pl.col("line") == "Lineage report").cum_sum()
            ) > 0
        ).filter(
            # remove everything from "Resistance report" onward
            (pl.col("line") != "Resistance report")
            .cum_min()
        ).filter(
            pl.col("line").str.strip_chars() != ""
        )

        # Add sample ID column (prepend like sed)
        df = df.with_columns(
            pl.lit(sample_id).alias("sample")
        ).select(
            ["sample", "line"]
        )

        dfs.append(df)

    # Concatenate all
    lineages_df = pl.concat(dfs)

    # Optional: sort like bash
    lineages_df = lineages_df.sort(["sample", "line"])

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    tbprofiler_df.write_csv(
        "tbprofiler.txt",
        separator="\t",
        include_header=True
    )

    # Output the split files (TBDB and WHO)
    tbprofiler_tbdb.write_csv(
        "tbdb-tbprofiler.txt",
        separator="\t",
        include_header=True
    )

    tbprofiler_who.write_csv(
        "who-tbprofiler.txt",
        separator="\t",
        include_header=True
    )

    # Output the lineage files
    lineages_df.write_csv(
        "lineages.fractions.txt",
        separator="\t",
        include_header=True
    )
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-11
@version: 1.1.0
@description:
    This process compiles the TB-Profiler results from the tbdb pipeline
    into a single file. Since TB-Profiler requires the results directory to be 
    present in the current working directory, we create symbolic links to the
    results directory, the bam directory and the vcf directory. The symbolic
    links are then used to run the tb-profiler collate command. The results
    are then moved to the current working directory and renamed to remove the
    tbdb- prefix. Renaming is to prevent clashes with input files in downstream
    processes. The results are then moved to the output directory.
@changelog:
        v1.0.0-2024-12-01: Initial version added
        v1.1.0-2025-04-11: Added - a handover from the TBPROFILER db updaitng module
*/