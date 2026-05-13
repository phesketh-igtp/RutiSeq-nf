#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
from Bio import Phylo

SINGLETON_VALUES = {"singleton", "singleton/singleton/singleton"}
OUTGROUPS = ("H37Rv", "MTB_anc")

def read_clusters(path: Path) -> pl.DataFrame:
    return pl.read_csv(path, separator="\t", infer_schema_length=10000)


def read_and_root_tree(tree_path: Path, outgroup: str = "MTB_anc"):
    tree = Phylo.read(str(tree_path), "newick")
    tip_names = {t.name for t in tree.get_terminals() if t.name is not None}
    if outgroup in tip_names:
        tree.root_with_outgroup(outgroup)
    return tree


def tree_tip_labels(tree) -> List[str]:
    return [t.name for t in tree.get_terminals() if t.name is not None]


def add_tip_label_column(df: pl.DataFrame, sample_col: str = "SampleID") -> pl.DataFrame:
    if "Tip_lable" in df.columns:
        return df
    return df.with_columns(pl.col(sample_col).alias("Tip_lable"))


def ensure_outgroup_rows(df: pl.DataFrame, sample_col: str = "SampleID") -> pl.DataFrame:
    schema = df.schema
    rows = []
    for og in OUTGROUPS:
        row = {c: None for c in df.columns}
        row[sample_col] = og
        if "Tip_lable" in df.columns:
            row["Tip_lable"] = og
        rows.append(row)
    add_df = pl.DataFrame(rows, schema=schema)
    return pl.concat([df, add_df], how="vertical").unique(subset=[sample_col])


def replace_singletons_with_nulls(df: pl.DataFrame, skip_cols: Optional[set[str]] = None) -> pl.DataFrame:
    skip_cols = skip_cols or set()
    exprs = []
    for c in df.columns:
        if c in skip_cols:
            exprs.append(pl.col(c))
        else:
            exprs.append(
                pl.when(pl.col(c).cast(pl.Utf8).is_in(list(SINGLETON_VALUES)))
                .then(None)
                .otherwise(pl.col(c))
                .alias(c)
            )
    return df.select(exprs)

def drop_lineage_for_matrix(df: pl.DataFrame) -> pl.DataFrame:
    if "lineage" in df.columns:
        return df.drop("lineage")
    return df

def stable_color(key: str) -> str:
    h = hashlib.md5(key.encode("utf-8")).hexdigest()
    r = int(h[0:2], 16)
    g = int(h[2:4], 16)
    b = int(h[4:6], 16)
    # brighten a bit so colors aren't too dark
    r = (r + 96) // 2
    g = (g + 96) // 2
    b = (b + 96) // 2
    return f"#{r:02x}{g:02x}{b:02x}"

def create_tree_palette_placeholder(
    tree_clusters_df: pl.DataFrame,
    lineage_short: str,
    id_col: str = "SampleID",
) -> Dict[str, Dict[str, str]]:
    palette: Dict[str, Dict[str, str]] = {}
    for col in tree_clusters_df.columns:
        if col == id_col:
            continue
        values = (
            tree_clusters_df.select(pl.col(col).cast(pl.Utf8))
            .drop_nulls()
            .unique()
            .to_series()
            .to_list()
        )
        mapping = {}
        for v in values:
            v_str = str(v)
            mapping[v_str] = stable_color(f"{lineage_short}|{col}|{v_str}")
        palette[col] = mapping
    return palette

def write_outputs(
    lineage_id: str,
    outdir: Path,
    tree,
    tree_clusters: pl.DataFrame,
    tree_clusters_df: pl.DataFrame,
    palette: Dict[str, Dict[str, str]],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    tree_newick = outdir / f"{lineage_id}.contree.newick"
    Phylo.write(tree, str(tree_newick), "newick")

    clusters_parquet = outdir / f"{lineage_id}.tree.clusters.parquet"
    matrix_parquet = outdir / f"{lineage_id}.tree.clusters.matrix.parquet"
    tree_clusters.write_parquet(clusters_parquet)
    tree_clusters_df.write_parquet(matrix_parquet)

    palette_json = outdir / f"{lineage_id}.palette.json"
    palette_json.write_text(json.dumps(palette, indent=2, sort_keys=True), encoding="utf-8")

    workspace = {
        "lineageID": lineage_id,
        "files": {
            "tree_newick": tree_newick.name,
            "tree_clusters_parquet": clusters_parquet.name,
            "tree_clusters_matrix_parquet": matrix_parquet.name,
            "palette_json": palette_json.name,
        },
        "outgroups": list(OUTGROUPS),
    }
    (outdir / f"{lineage_id}.contree.workspace.json").write_text(
        json.dumps(workspace, indent=2),
        encoding="utf-8",
    )

def main():
    ap = argparse.ArgumentParser(description="Polars port of plot_ML-phylogeny.R (data prep only)")
    ap.add_argument("--lineageID", required=True, help="e.g. lineage4.8")
    ap.add_argument("--clusters", default="clusters.tsv", help="Path to clusters.tsv")
    ap.add_argument("--tree", default="snp.contree", help="Path to Newick tree (snp.contree)")
    ap.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    args = ap.parse_args()

    lineage_id = args.lineageID
    lineage_short = lineage_id.replace("lineage", "L")

    clusters_path = Path(args.clusters)
    tree_path = Path(args.tree)
    outdir = Path(args.outdir)

    # Read inputs
    clusters = read_clusters(clusters_path)
    tree = read_and_root_tree(tree_path, outgroup="MTB_anc")

    # R: filtered_clusters <- clusters |> filter(lineage == lineageID)
    filtered_clusters = clusters
    if "lineage" in clusters.columns:
        filtered_clusters = clusters.filter(pl.col("lineage") == lineage_id)

    # Build tree-specific metadata
    tips = tree_tip_labels(tree)

    # R: tree.clusters <- clusters %>% filter(SampleID %in% tree.tips) |> distinct()
    if "SampleID" not in clusters.columns:
        raise ValueError("clusters.tsv must contain a 'SampleID' column")

    tree_clusters = (
        clusters
        .filter(pl.col("SampleID").is_in(tips))
        .unique()
    )

    # R: mutate(Tip_lable = SampleID)
    tree_clusters = add_tip_label_column(tree_clusters, sample_col="SampleID")

    # Add outgroup rows (H37Rv + MTB_anc)
    tree_clusters = ensure_outgroup_rows(tree_clusters, sample_col="SampleID")

    # Replace singletons with null in both versions
    tree_clusters = replace_singletons_with_nulls(tree_clusters, skip_cols={"SampleID"})

    tree_clusters_df = drop_lineage_for_matrix(tree_clusters)

    # If you want to mimic R's "column_to_rownames(SampleID)", keep SampleID as key column.
    # Downstream code can join on SampleID instead of relying on rownames.

    # Palette (placeholder for create_tree_palette)
    palette = create_tree_palette_placeholder(tree_clusters_df, lineage_short=lineage_short, id_col="SampleID")

    # Write outputs (Python workspace replacing RData)
    write_outputs(
        lineage_id=lineage_id,
        outdir=outdir,
        tree=tree,
        tree_clusters=tree_clusters,
        tree_clusters_df=tree_clusters_df,
        palette=palette,
    )


if __name__ == "__main__":
    main()