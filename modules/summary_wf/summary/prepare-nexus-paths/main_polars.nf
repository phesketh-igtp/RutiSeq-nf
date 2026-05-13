process PREPARE_NEXUS_PATHS{

    conda params.snippy_env 

    tag "${lineage}; t=${distance}"

    publishDir "${params.outDir}/results/networks/${lineage}/", 
        mode: 'copy', 
        overwrite: true

    errorStrategy 'ignore'

    input:
        tuple val(lineage),
            val(distance),
            path(snp_fasta),
            path(snps_tab),
            path(clusters_tab)

    output:
        path("nexus.tuple.csv"), emit: nexus_tuple
        path(clusters_tab),      emit: clusters_tab

    script:

    def snp_fasta_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.fasta"
    def snp_tab_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Amend/${lineage}_joint_cf${params.mtbseq_mincovf}_cr${params.mtbseq_mincovr}_fr${params.mtbseq_minfreq}_ph${params.mtbseq_minphred20}_samples*_amended_u${params.mtbseq_unambig}_phylo_w${params.mtbseq_window}.tab"
    def clusters_tab_path="${params.outDir}/db/comparison/mtbseq/${lineage}/Groups/${params.runID}_${lineage}_d${distance}.processed.clusters.tsv"

    """
    #!/usr/bin/env pythons
    import polars as pl

    # Inputs (replace with your actual paths/values)
    clusters_tab = "clusters.tsv"
    lineage = "${lineage}"
    snp_fasta_path = "${snp_fasta_path}"
    snp_tab_path = "${snp_tab_path}"
    clusters_tab_path = "${clusters_tab_path}"

    # Load table (assuming tab-delimited, no header)
    df = pl.read_csv(clusters_tab_path, separator="\t", has_header=False)

    # Assign column names (based on your usage: col1 = genome ID, col4 = cluster ID)
    df = df.rename({
        "column_1": "genome",
        "column_4": "cluster"
    })

    # --------------------------
    # Step 1: Identify clusters containing the lineage
    # --------------------------
    # Equivalent of: grep -w lineage | cut -f4 | sort | uniq
    unique_clusters = (
        df.filter(pl.col("genome") == lineage)
        .select("cluster")
        .unique()
    )

    # --------------------------
    # Step 2: Keep clusters with >= 3 genomes
    # --------------------------
    # Count occurrences of each cluster
    cluster_counts = df.groupby("cluster").count()

    # Join with unique clusters and filter
    final_clusters = (
        unique_clusters.join(cluster_counts, on="cluster")
                    .filter(pl.col("count") >= 3)
                    .select("cluster")
    )

    # --------------------------
    # Step 3: Handle empty case
    # --------------------------
    if final_clusters.height == 0:
        print(f"No clusters found for lineage {lineage}. Creating empty nexus file.")
        pl.DataFrame().write_csv("nexus.tuple.csv")
    else:
        rows = []
        for cluster_id in final_clusters["cluster"]:
            # Get genomes for this cluster and join with ';'
            genomes = (
                df.filter(pl.col("cluster") == cluster_id)
                .select("genome")
                .to_series()
                .to_list()
            )
            genome_str = ";".join(genomes)
            rows.append({
                "lineage": lineage,
                "cluster": cluster_id,
                "snp_fasta_path": snp_fasta_path,
                "snp_tab_path": snp_tab_path,
                "clusters_tab_path": clusters_tab_path,
            })
    result_df = pl.DataFrame(rows)

    # --------------------------
    # Step 4: Remove clusters starting with "nX-"
    # --------------------------
    result_df = result_df.filter(~pl.col("cluster").str.starts_with("nX-"))

    # Write output CSV
    result_df.write_csv("nexus.tuple.csv")
    """
}

/*
@author: Poppy J Hesketh Best
@date: 2025-04-01
@version: 1.0.2
@description:
    This process prepares the paths for the NEXUS files for each cluster.
    It generates a CSV file with the paths to the NEXUS files and the
    corresponding tab files. The CSV file is used as input for the
    GENERATE_NEXUS process.
    The tuple has the following format:
    ["lineage", "clusterID", "fasta_path", "tab_path"]
@changelog:
    v1.0.0-2025-04-01: Initial version
    v1.0.1-2025-04-04: Added filtering to remove clusters with less than 3 genomes
    v1.0.2-2025-12-01: Updated the paths for the new system.
*/