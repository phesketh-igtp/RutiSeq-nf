process PREPROCESS_CLUSTER {

/*
    @author: Poppy J Hesketh Best
    @date: 2025-04-01
    @version: 1.0.1
    @description:
        This process runs the MTBseq TBgroups step on the joint and amend directories
        for each lineage and distance. It takes the output from the
        MTBSEQ_LINEAGE_JOINT_AMEND() process and runs the TBgroups step on the joint
        and amend directories. It also renames the output files for simplicity.
        It also wrangles the output matrix into a useful format for downstream analysis.
        Occasionally, the MTBseq TBgroups step will fail and produice empty files/no-files.
    @last_updated: 2025-04-01
    @changelog:
        v1.0.0-2025-04-01: Initial version + documnetadocumentation
        v1.0.1-2025-05-19: Removed additonal_args as it was not utilised and was creating inconsistencies
*/

    tag "${lineage}; t=${distance}"

    conda params.r_stats_env

    // TODO: container { if (workflow.containerEngine == 'singularity') { 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ce098dd570838fdcb0eb401b3afe4ebf4bc88d1038768ec18b3f970deb28c313/data'
    ///        } else { 'quay.io/biocontainers/mtbseq' }
    ///}
                
    publishDir "${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/", 
        mode: 'copy',
        overwrite: true

    input:
        // Nexus output
        tuple val(lineage), 
            val(distance), 
            path(fasta),
            path(tab),
            path(clusters, stageAs: 'clusters_input.tsv'),
            path(singletons)

    output:
        // Nexus output
        tuple val(lineage), 
            val(distance), 
            path("${params.runID}_${lineage}_snps.fasta"),
            path("${params.runID}_${lineage}_snps.tab"),
            path("${params.runID}_${lineage}_d${distance}.processed.clusters.tsv"), emit: nexus_ch

        path("${params.runID}_${lineage}_d${distance}.processed.clusters.tsv"), emit: pairwise_clusters_processed

    script:

        //def additional_args = task.ext.additional_args ?: '' // defined in the nextflow.config file

"""
import polars as pl

# ------------------------------------------------------------------
# Read input
# ------------------------------------------------------------------
df = pl.read_csv(
    "clusters_input.tsv",
    separator="\t",
    has_header=False
)

# ------------------------------------------------------------------
# Process clusters
# ------------------------------------------------------------------
df = (
    df
    .unique()
    .with_columns(
        pl.col("column_4")
        .str.split_exact("-", 2)
        .struct.rename_fields(["num", "dist", "lin"])
        .alias("parts")
    )
    .unnest("parts")
    .with_columns([
        pl.col("num").str.zfill(3),
        pl.col("dist").str.zfill(2)
    ])
    .with_columns(
        (pl.col("num") + "-" +
         pl.col("dist") + "-" +
         pl.col("lin"))
        .alias("column_4")
    )
    .select(["column_1", "column_2", "column_3", "column_4"])
)

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------
df.write_csv(
    "${params.runID}_${lineage}_d${distance}.processed.clusters.tsv",
    separator="\t",
    include_header=False
)
"""
}