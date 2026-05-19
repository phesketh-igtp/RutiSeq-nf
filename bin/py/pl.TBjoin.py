#!/usr/bin/env python3

"""
tbjoin_polars.py

First-draft Python/Polars rewrite of MTBSeq TBjoin.

Purpose
-------
This script is designed to replace the MTBSeq TBjoin stage while preserving
the MTBSeq-compatible Join table required by downstream MTBSeq modules such
as TBamend.

Main performance improvement
----------------------------
The original Perl implementation repeatedly rewrites the growing Join table
once per sample. This script instead:

    1. Builds the global variant scaffold.
    2. Processes each sample into a compact per-sample call file.
    3. Writes the final wide MTBSeq-compatible Join table once.

This should scale much better for large numbers of genomes.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import polars as pl

# pip install numpy polars

# =============================================================================
# Constants
# =============================================================================

POSITION_COLUMNS = [
    "pos",
    "insindex",
    "ref",
    "A",
    "C",
    "G",
    "T",
    "N",
    "GAP",
    "a",
    "c",
    "g",
    "t",
    "n",
    "gap",
    "A_qual_20",
    "C_qual_20",
    "G_qual_20",
    "T_qual_20",
    "N_qual_20",
    "GAP_qual_20",
]

NUMERIC_POSITION_COLUMNS = [c for c in POSITION_COLUMNS if c != "ref"]

CALL_FIELDS = [
    "Type",
    "Allel",
    "CovFor",
    "CovRev",
    "Qual20",
    "Freq",
    "Cov",
    "Subst",
]

DEFAULT_CALL = ["-", "-", "0", "0", "0", "0", "0", "-"]

CODON_TABLE = {
    "TCA": "Ser", "TCC": "Ser", "TCG": "Ser", "TCT": "Ser",
    "TTC": "Phe", "TTT": "Phe", "TTA": "Leu", "TTG": "Leu",
    "TAC": "Tyr", "TAT": "Tyr", "TAA": "_",   "TAG": "_",
    "TGA": "_",   "TGC": "Cys", "TGT": "Cys", "TGG": "Trp",
    "CTA": "Leu", "CTC": "Leu", "CTG": "Leu", "CTT": "Leu",
    "CCA": "Pro", "CCC": "Pro", "CCG": "Pro", "CCT": "Pro",
    "CAC": "His", "CAT": "His", "CAA": "Gln", "CAG": "Gln",
    "CGA": "Arg", "CGC": "Arg", "CGG": "Arg", "CGT": "Arg",
    "ATA": "Ile", "ATC": "Ile", "ATT": "Ile", "ATG": "Met",
    "ACA": "Thr", "ACC": "Thr", "ACG": "Thr", "ACT": "Thr",
    "AAC": "Asn", "AAT": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGC": "Ser", "AGT": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GTA": "Val", "GTC": "Val", "GTG": "Val", "GTT": "Val",
    "GCA": "Ala", "GCC": "Ala", "GCG": "Ala", "GCT": "Ala",
    "GAC": "Asp", "GAT": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGA": "Gly", "GGC": "Gly", "GGG": "Gly", "GGT": "Gly",
}


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class Gene:
    id: str
    name: str
    start: int
    stop: int
    frame: str
    product: str
    description: str
    cogfunction: str
    cogcats: str
    status_region_name: str
    status_function_name: str
    type: str


@dataclass(frozen=True)
class ScaffoldRow:
    pos: int
    insindex: int
    ref: str
    gene: str
    gene_name: str
    product: str


@dataclass
class VariantCall:
    ref: str
    type: str
    allele: str
    cov_f: int
    cov_r: int
    qual20: int
    freq: str
    cov: int
    subst: str
    gene: str
    product: str
    count: int


# =============================================================================
# General helpers
# =============================================================================

def fmt_param_value(x) -> str:
    """
    Format parameters like Perl string interpolation.

    In particular:
        75.0 -> 75

    This prevents output names like:
        fr75.0

    when MTBSeq expects:
        fr75
    """
    try:
        xf = float(x)
        if xf.is_integer():
            return str(int(xf))
        return str(x)
    except Exception:
        return str(x)


def build_output_names(
    group_name: str,
    micovf: int,
    micovr: int,
    mifreq: float,
    miphred20: int,
    nsamples: int,
) -> Tuple[str, str]:
    param_string = (
        f"_cf{fmt_param_value(micovf)}"
        f"_cr{fmt_param_value(micovr)}"
        f"_fr{fmt_param_value(mifreq)}"
        f"_ph{fmt_param_value(miphred20)}"
    )
    base = f"{group_name}_joint{param_string}_samples{nsamples}"
    return f"{base}.tab", f"{base}.log"


def parse_fasta(path: Path) -> str:
    """
    Equivalent to TBtools::parse_fasta.

    - Ignores FASTA headers.
    - Concatenates all sequence lines.
    - Returns a single genome string.
    """
    seq = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.startswith(">"):
                seq.append(line)
    return "".join(seq)


def reverse_complement(dna: str) -> str:
    return dna[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def codon2aa(codon: str) -> str:
    codon = codon.upper()
    return CODON_TABLE.get(codon, "Bad")


def perl_freq(numer: float, denom: float) -> str:
    """
    Match Perl sprintf("%.2f", numer / denom * 100) where possible.

    If denom is zero, return 0 rather than crashing.
    The original Perl code assumes this should not happen for insertions,
    but guarding here is safer.
    """
    if denom == 0:
        return "0"
    return f"{(numer / denom * 100):.2f}"


def format_freq_for_join(freq) -> str:
    """
    Match print_joint_table frequency formatting:

        $freq1 = sprintf("%.2f", $freq1) unless($freq1 == 0);

    So:
        0       -> "0"
        100     -> "100.00"
        100.00  -> "100.00"
    """
    try:
        x = float(freq)
    except Exception:
        return str(freq)

    if x == 0:
        return "0"

    return f"{x:.2f}"


def extract_sample_id(var_file: Path) -> str:
    """
    Equivalent to Perl:

        $file =~ /^(.*)\\.gatk_position_variants_.*\\.tab$/;
        my $id = $1;
    """
    name = var_file.name
    m = re.match(r"^(.*)\.gatk_position_variants_.*\.tab$", name)
    if m:
        return m.group(1)
    return name


def position_file_for_sample(sample_id: str) -> str:
    return f"{sample_id}.gatk_position_table.tab"


# =============================================================================
# Annotation parsing
# =============================================================================

def parse_annotation(path: Path) -> Tuple[Dict[str, Gene], Dict[int, str]]:
    """
    Parse MTBSeq custom annotation file.

    Expected columns:
        0  id
        1  name
        2  start
        3  stop
        4  frame
        5  product
        6  description
        7  cogfunction
        8  cogcats
        9  status_region_name
        10 status_function_name
        11 gene_type

    Handles:
        - comment lines beginning with #
        - blank lines
        - plain header lines such as:
              id name start stop ...
        - overlapping annotations by overwriting with later entries, as in Perl
    """
    genes: Dict[str, Gene] = {}
    annotation_by_pos: Dict[int, str] = {}

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\r\n")

            if not line:
                continue

            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                continue

            # Skip header rows such as:
            # id    name    start    stop ...
            try:
                start = int(fields[2])
                stop = int(fields[3])
            except ValueError:
                continue

            if len(fields) < 12:
                fields += [" "] * (12 - len(fields))

            gene = Gene(
                id=fields[0],
                name=fields[1],
                start=start,
                stop=stop,
                frame=fields[4],
                product=fields[5],
                description=fields[6],
                cogfunction=fields[7],
                cogcats=fields[8],
                status_region_name=fields[9],
                status_function_name=fields[10],
                type=fields[11],
            )

            genes[gene.id] = gene

            lo = min(gene.start, gene.stop)
            hi = max(gene.start, gene.stop)

            for pos in range(lo, hi + 1):
                annotation_by_pos[pos] = gene.id

    return genes, annotation_by_pos

# =============================================================================
# Variant file parsing and global scaffold
# =============================================================================

def parse_variant_files(
    call_out: Path,
    var_files: List[Path],
    micovf: int,
    micovr: int,
    mifreq: float,
    miphred20: int,
    snp_vars: int,
) -> Tuple[pl.DataFrame, Set[Tuple[int, int, str]], List[str]]:
    """
    Parse per-sample variant files.

    This reproduces TBtools::parse_variants for the parts needed by TBjoin.

    Returns:
        variant_index:
            Polars DataFrame with columns:
                pos, insindex, ref

        original_variant_keys:
            set of:
                (pos, insindex, sample_id)

            This replaces the Perl $strain existence marker.

        ids:
            sample IDs in Perl-compatible sorted order.
    """
    rows = []
    original_keys: Set[Tuple[int, int, str]] = set()
    ids: List[str] = []

    for var_file in sorted(var_files, key=lambda p: p.name):
        sample_id = extract_sample_id(var_file)
        ids.append(sample_id)

        path = call_out / var_file.name

        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")

            try:
                next(reader)
            except StopIteration:
                continue

            for fields in reader:
                if not fields:
                    continue

                fields += [""] * max(0, 17 - len(fields))

                pos = int(fields[0])
                insindex = int(fields[1])
                ref_tmp = fields[2]
                vtype = fields[3]
                allele = fields[4]

                try:
                    covf = float(fields[5])
                    covr = float(fields[6])
                    qual20 = float(fields[7])
                    freq1 = float(fields[8])
                except ValueError:
                    continue

                if snp_vars == 1 and vtype != "SNP":
                    continue

                unambiguous = (
                    covf >= micovf
                    and covr >= micovr
                    and freq1 >= mifreq
                    and qual20 >= miphred20
                    and (
                        re.search(r"[ACGTacgt]", allele) is not None
                        or "GAP" in allele
                    )
                )

                if not unambiguous:
                    continue

                rows.append((pos, insindex, ref_tmp))
                original_keys.add((pos, insindex, sample_id))

    if rows:
        variant_index = (
            pl.DataFrame(
                rows,
                schema=["pos", "insindex", "ref"],
                orient="row",
            )
            .unique(subset=["pos", "insindex"], keep="first")
            .sort(["pos", "insindex"])
        )
    else:
        variant_index = pl.DataFrame(
            {
                "pos": pl.Series([], dtype=pl.Int64),
                "insindex": pl.Series([], dtype=pl.Int64),
                "ref": pl.Series([], dtype=pl.Utf8),
            }
        )

    return variant_index, original_keys, ids


def build_scaffold(
    variant_index: pl.DataFrame,
    genes: Dict[str, Gene],
    annotation_by_pos: Dict[int, str],
) -> List[ScaffoldRow]:
    """
    Equivalent to the row part of print_joint_table_scaffold.
    """
    scaffold: List[ScaffoldRow] = []

    for row in variant_index.sort(["pos", "insindex"]).iter_rows(named=True):
        pos = int(row["pos"])
        insindex = int(row["insindex"])
        ref = row["ref"] if row["ref"] is not None else " "

        gene_id = annotation_by_pos.get(pos, " ")
        gene_name = " "
        product = " "

        if gene_id in genes:
            gene_name = genes[gene_id].name
            product = genes[gene_id].product

        scaffold.append(
            ScaffoldRow(
                pos=pos,
                insindex=insindex,
                ref=ref,
                gene=gene_id,
                gene_name=gene_name,
                product=product,
            )
        )

    return scaffold


def build_required_positions(variant_index: pl.DataFrame) -> pl.DataFrame:
    """
    For insertion calls, TBtools::call_variants needs the parent row:
        pos, insindex == 0

    even if the final scaffold row is:
        pos, insindex != 0
    """
    base = variant_index.select(["pos", "insindex"])

    insertion_parents = (
        variant_index
        .filter(pl.col("insindex") != 0)
        .select("pos")
        .with_columns(pl.lit(0, dtype=pl.Int64).alias("insindex"))
        .select(["pos", "insindex"])
    )

    return (
        pl.concat([base, insertion_parents])
        .unique(subset=["pos", "insindex"])
        .sort(["pos", "insindex"])
    )


# =============================================================================
# Position table handling
# =============================================================================

def scan_position_table(path: Path) -> pl.LazyFrame:
    """
    Scan an MTBSeq position table.

    MTBSeq headers are inconsistent:
        #Pos
        RefBase
        Aqual_20
        Cqual_20
        G_qual_20
        Tqual_20
        Nqual_20
        GAPqual_20

    We override names by position to avoid relying on exact header spelling.
    """
    lf = pl.scan_csv(
        path,
        separator="\t",
        has_header=True,
        new_columns=POSITION_COLUMNS,
        infer_schema=False,
        truncate_ragged_lines=True,
        ignore_errors=False,
    )

    cast_exprs = []

    for col in NUMERIC_POSITION_COLUMNS:
        cast_exprs.append(
            pl.col(col)
            .cast(pl.Int64, strict=False)
            .fill_null(0)
            .alias(col)
        )

    cast_exprs.append(
        pl.col("ref")
        .fill_null(" ")
        .cast(pl.Utf8)
        .alias("ref")
    )

    return lf.with_columns(cast_exprs)


def add_basic_count_columns(lf: pl.LazyFrame) -> pl.LazyFrame:
    return (
        lf.with_columns(
            [
                (pl.col("A") + pl.col("a")).alias("adenosin"),
                (pl.col("C") + pl.col("c")).alias("cytosin"),
                (pl.col("G") + pl.col("g")).alias("guanosin"),
                (pl.col("T") + pl.col("t")).alias("thymin"),
                (pl.col("N") + pl.col("n")).alias("nucleosin"),
                (pl.col("GAP") + pl.col("gap")).alias("gaps"),
            ]
        )
        .with_columns(
            (
                pl.col("adenosin")
                + pl.col("cytosin")
                + pl.col("guanosin")
                + pl.col("thymin")
                + pl.col("nucleosin")
                + pl.col("gaps")
            ).alias("cov")
        )
    )


def update_breadth_arrays(
    lf: pl.LazyFrame,
    any_cov: np.ndarray,
    unambig: np.ndarray,
    nothing: np.ndarray,
    micovf: int,
    micovr: int,
    miphred20: int,
    mifreq: float,
) -> None:
    """
    Reproduce the position_stats part of TBtools::parse_position_table.

    Only insindex == 0 contributes to breadth.
    """
    lf2 = add_basic_count_columns(lf).filter(pl.col("insindex") == 0)

    lf2 = lf2.with_columns(
        [
            pl.when(pl.col("cov") > 0)
            .then(pl.col("adenosin") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("A_freq"),

            pl.when(pl.col("cov") > 0)
            .then(pl.col("cytosin") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("C_freq"),

            pl.when(pl.col("cov") > 0)
            .then(pl.col("guanosin") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("G_freq"),

            pl.when(pl.col("cov") > 0)
            .then(pl.col("thymin") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("T_freq"),

            pl.when(pl.col("cov") > 0)
            .then(pl.col("nucleosin") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("N_freq"),

            pl.when(pl.col("cov") > 0)
            .then(pl.col("gaps") / pl.col("cov") * 100)
            .otherwise(0)
            .alias("GAP_freq"),
        ]
    )

    unamb_expr = (
        (
            (pl.col("A") >= micovf)
            & (pl.col("a") >= micovr)
            & (pl.col("A_qual_20") >= miphred20)
            & (pl.col("A_freq") >= mifreq)
        )
        | (
            (pl.col("C") >= micovf)
            & (pl.col("c") >= micovr)
            & (pl.col("C_qual_20") >= miphred20)
            & (pl.col("C_freq") >= mifreq)
        )
        | (
            (pl.col("G") >= micovf)
            & (pl.col("g") >= micovr)
            & (pl.col("G_qual_20") >= miphred20)
            & (pl.col("G_freq") >= mifreq)
        )
        | (
            (pl.col("T") >= micovf)
            & (pl.col("t") >= micovr)
            & (pl.col("T_qual_20") >= miphred20)
            & (pl.col("T_freq") >= mifreq)
        )
        | (
            (pl.col("GAP") >= micovf)
            & (pl.col("gap") >= micovr)
            & (pl.col("GAP_qual_20") >= miphred20)
            & (pl.col("GAP_freq") >= mifreq)
        )
    )

    stats = (
        lf2.select(
            [
                "pos",
                (pl.col("cov") == 0).alias("is_nothing"),
                (pl.col("cov") > 0).alias("is_any"),
                unamb_expr.alias("is_unambig"),
            ]
        )
        .collect()
    )

    if stats.height == 0:
        return

    positions = stats["pos"].to_numpy()

    is_any = stats["is_any"].to_numpy()
    is_nothing = stats["is_nothing"].to_numpy()
    is_unambig = stats["is_unambig"].to_numpy()

    any_cov[positions[is_any]] += 1
    nothing[positions[is_nothing]] += 1
    unambig[positions[is_unambig]] += 1


# =============================================================================
# Variant calling
# =============================================================================

def choose_majority_base(row: dict) -> Tuple[str, int, int, int, int]:
    """
    Reproduce the actual Perl tie behaviour.

    Perl comment says:
        A < C < G < T < N < GAP

    But the code uses independent if statements, not elsif.
    Therefore later bases overwrite earlier bases.

    Effective tie priority:
        GAP > N > T > G > C > A
    """
    counts = {
        "A": row["adenosin"],
        "C": row["cytosin"],
        "G": row["guanosin"],
        "T": row["thymin"],
        "N": row["nucleosin"],
        "GAP": row["gaps"],
    }

    maximum = max(counts.values())

    allele = " "
    count = 0
    qual20 = 0
    cov_f = 0
    cov_r = 0

    for base in ["A", "C", "G", "T", "N", "GAP"]:
        if counts[base] == maximum:
            allele = base
            count = int(counts[base])

            if base == "A":
                qual20 = int(row["A_qual_20"])
                cov_f = int(row["A"])
                cov_r = int(row["a"])
            elif base == "C":
                qual20 = int(row["C_qual_20"])
                cov_f = int(row["C"])
                cov_r = int(row["c"])
            elif base == "G":
                qual20 = int(row["G_qual_20"])
                cov_f = int(row["G"])
                cov_r = int(row["g"])
            elif base == "T":
                qual20 = int(row["T_qual_20"])
                cov_f = int(row["T"])
                cov_r = int(row["t"])
            elif base == "N":
                qual20 = int(row["N_qual_20"])
                cov_f = int(row["N"])
                cov_r = int(row["n"])
            elif base == "GAP":
                qual20 = int(row["GAP_qual_20"])
                cov_f = int(row["GAP"])
                cov_r = int(row["gap"])

    return allele, count, qual20, cov_f, cov_r


def choose_lowfreq_base(
    row: dict,
    micovf: int,
    micovr: int,
    miphred20: int,
    mifreq: float,
) -> Optional[Tuple[str, int, int, int, int]]:
    """
    Low-frequency mode equivalent to TBtools::call_variants.

    Choose highest-frequency non-reference base passing thresholds.
    Tie behaviour follows actual Perl overwrite order.
    """
    ref = row["ref"]
    cov = row["cov"]

    candidate_counts = []

    def maybe_add(base: str, total: int, fwd: int, rev: int, q20: int) -> None:
        if base != ref and fwd >= micovf and rev >= micovr and q20 >= miphred20:
            if cov > 0 and total / cov * 100 >= mifreq:
                candidate_counts.append(total)

    maybe_add("A", row["adenosin"], row["A"], row["a"], row["A_qual_20"])
    maybe_add("C", row["cytosin"], row["C"], row["c"], row["C_qual_20"])
    maybe_add("G", row["guanosin"], row["G"], row["g"], row["G_qual_20"])
    maybe_add("T", row["thymin"], row["T"], row["t"], row["T_qual_20"])
    maybe_add("N", row["nucleosin"], row["N"], row["n"], row["N_qual_20"])
    maybe_add("GAP", row["gaps"], row["GAP"], row["gap"], row["GAP_qual_20"])

    if not candidate_counts:
        return None

    maximum_low = max(candidate_counts)

    allele = None
    count = 0
    qual20 = 0
    cov_f = 0
    cov_r = 0

    counts = {
        "A": row["adenosin"],
        "C": row["cytosin"],
        "G": row["guanosin"],
        "T": row["thymin"],
        "N": row["nucleosin"],
        "GAP": row["gaps"],
    }

    for base in ["A", "C", "G", "T", "N", "GAP"]:
        if counts[base] == maximum_low:
            allele = base
            count = int(counts[base])

            if base == "A":
                qual20 = int(row["A_qual_20"])
                cov_f = int(row["A"])
                cov_r = int(row["a"])
            elif base == "C":
                qual20 = int(row["C_qual_20"])
                cov_f = int(row["C"])
                cov_r = int(row["c"])
            elif base == "G":
                qual20 = int(row["G_qual_20"])
                cov_f = int(row["G"])
                cov_r = int(row["g"])
            elif base == "T":
                qual20 = int(row["T_qual_20"])
                cov_f = int(row["T"])
                cov_r = int(row["t"])
            elif base == "N":
                qual20 = int(row["N_qual_20"])
                cov_f = int(row["N"])
                cov_r = int(row["n"])
            elif base == "GAP":
                qual20 = int(row["GAP_qual_20"])
                cov_f = int(row["GAP"])
                cov_r = int(row["gap"])

    if allele is None:
        return None

    return allele, count, qual20, cov_f, cov_r


def substitution_string(
    pos: int,
    allele: str,
    vtype: str,
    genes: Dict[str, Gene],
    annotation_by_pos: Dict[int, str],
    ref_genome: str,
) -> Tuple[str, int]:
    """
    Reproduce the SNP substitution string from TBtools::call_variants.

    Returns:
        substitution string
        real_substitution flag
    """
    gene_id = annotation_by_pos.get(pos, " ")
    gene = genes.get(gene_id)

    if vtype != "SNP" or gene is None or gene.type != "CDS":
        return " ", 0

    start = gene.start
    stop = gene.stop

    if start < stop:
        length = stop - start + 1
        gene_pos = pos - start + 1
        rank = int((gene_pos + 2) / 3)

        dna = ref_genome[start - 1 : start - 1 + length]

        codon_start = (rank - 1) * 3
        codon = dna[codon_start : codon_start + 3]
        aa = codon2aa(codon)

        snp_dna = dna[: gene_pos - 1] + allele + dna[gene_pos:]
        snp_codon = snp_dna[codon_start : codon_start + 3]

        snp_aa = "bad"
        if not re.search(r"[Nn]", snp_codon):
            snp_aa = codon2aa(snp_codon)

        real_sub = 0 if aa == snp_aa else 1
        return f"{aa}{rank}{snp_aa} ({codon}/{snp_codon})", real_sub

    if start > stop:
        length = start - stop + 1
        gene_pos = abs(pos - start) + 1
        rank = int((gene_pos + 2) / 3)

        dna = ref_genome[stop - 1 : stop - 1 + length]
        rev_dna = reverse_complement(dna)
        rev_allele = reverse_complement(allele)

        codon_start = (rank - 1) * 3
        codon = rev_dna[codon_start : codon_start + 3]
        aa = codon2aa(codon)

        snp_dna = rev_dna[: gene_pos - 1] + rev_allele + rev_dna[gene_pos:]
        snp_codon = snp_dna[codon_start : codon_start + 3]

        snp_aa = "bad"
        if not re.search(r"[Nn]", snp_codon):
            snp_aa = codon2aa(snp_codon)

        real_sub = 0 if aa == snp_aa else 1
        return f"{aa}{rank}{snp_aa} ({codon}/{snp_codon})", real_sub

    return " ", 0


def call_one_row(
    row: dict,
    parent_cov_by_pos: Dict[int, int],
    genes: Dict[str, Gene],
    annotation_by_pos: Dict[int, str],
    ref_genome: str,
    micovf: int,
    micovr: int,
    miphred20: int,
    mifreq: float,
    lowfreq_vars: int,
) -> VariantCall:
    """
    Row-wise equivalent of TBtools::call_variants for a single row.
    """
    ref_tmp = row["ref"]
    pos = int(row["pos"])
    insindex = int(row["insindex"])
    cov = int(row["cov"])

    allele = " "
    freq1 = "0"
    count1 = 0
    qual20 = 0
    cov_f = 0
    cov_r = 0
    vtype = "none"

    if lowfreq_vars == 1 and insindex == 0 and cov > 0:
        low = choose_lowfreq_base(row, micovf, micovr, miphred20, mifreq)

        if low is not None:
            allele, count1, qual20, cov_f, cov_r = low
            freq1 = perl_freq(count1, cov)

            if allele == "GAP":
                vtype = "Del"
            elif allele != ref_tmp:
                vtype = "SNP"

    if cov > 0 and allele == " ":
        allele, count1, qual20, cov_f, cov_r = choose_majority_base(row)
        freq1 = perl_freq(count1, cov)

        if allele == "GAP":
            vtype = "Del"
        elif allele != ref_tmp:
            vtype = "SNP"

    if cov <= 0:
        vtype = "Unc"
        allele = "U"

    if insindex != 0:
        vtype = "Ins"
        parent_cov = parent_cov_by_pos.get(pos, 0)
        qual20 = cov_f + cov_r

        if allele == "A":
            freq1 = perl_freq(row["adenosin"], parent_cov)
        elif allele == "C":
            freq1 = perl_freq(row["cytosin"], parent_cov)
        elif allele == "G":
            freq1 = perl_freq(row["guanosin"], parent_cov)
        elif allele == "T":
            freq1 = perl_freq(row["thymin"], parent_cov)
        elif allele == "N":
            freq1 = perl_freq(row["nucleosin"], parent_cov)

    subst, _real_sub = substitution_string(
        pos=pos,
        allele=allele,
        vtype=vtype,
        genes=genes,
        annotation_by_pos=annotation_by_pos,
        ref_genome=ref_genome,
    )

    gene_id = annotation_by_pos.get(pos, " ")
    product = genes[gene_id].product if gene_id in genes else " "

    return VariantCall(
        ref=ref_tmp,
        type=vtype,
        allele=allele,
        cov_f=int(cov_f),
        cov_r=int(cov_r),
        qual20=int(qual20),
        freq=str(freq1),
        cov=int(cov),
        subst=subst,
        gene=gene_id,
        product=product,
        count=int(count1),
    )


# =============================================================================
# Per-sample processing
# =============================================================================

def process_sample(
    sample_id: str,
    pos_out: Path,
    sample_calls_dir: Path,
    scaffold: List[ScaffoldRow],
    required_positions: pl.DataFrame,
    original_variant_keys: Set[Tuple[int, int, str]],
    genes: Dict[str, Gene],
    annotation_by_pos: Dict[int, str],
    ref_genome: str,
    any_cov: np.ndarray,
    unambig: np.ndarray,
    nothing: np.ndarray,
    micovf: int,
    micovr: int,
    miphred20: int,
    mifreq: float,
    lowfreq_vars: int,
) -> Path:
    """
    Process one sample.

    Output:
        work/sample_calls/<sample>.calls.tsv

    Each line corresponds exactly to one scaffold row and contains:
        Type Allel CovFor CovRev Qual20 Freq Cov Subst
    """
    position_path = pos_out / position_file_for_sample(sample_id)
    out_path = sample_calls_dir / f"{sample_id}.calls.tsv"

    if not position_path.exists():
        with out_path.open("w") as out:
            for _ in scaffold:
                out.write("\t".join(DEFAULT_CALL) + "\n")
        return out_path

    lf = scan_position_table(position_path)

    update_breadth_arrays(
        lf=lf,
        any_cov=any_cov,
        unambig=unambig,
        nothing=nothing,
        micovf=micovf,
        micovr=micovr,
        miphred20=miphred20,
        mifreq=mifreq,
    )

    filtered = (
        add_basic_count_columns(lf)
        .join(required_positions.lazy(), on=["pos", "insindex"], how="inner")
        .collect()
    )

    parent_cov_by_pos: Dict[int, int] = {}

    for row in filtered.filter(pl.col("insindex") == 0).iter_rows(named=True):
        parent_cov_by_pos[int(row["pos"])] = int(row["cov"])

    call_map: Dict[Tuple[int, int], VariantCall] = {}

    for row in filtered.iter_rows(named=True):
        key = (int(row["pos"]), int(row["insindex"]))

        call_map[key] = call_one_row(
            row=row,
            parent_cov_by_pos=parent_cov_by_pos,
            genes=genes,
            annotation_by_pos=annotation_by_pos,
            ref_genome=ref_genome,
            micovf=micovf,
            micovr=micovr,
            miphred20=miphred20,
            mifreq=mifreq,
            lowfreq_vars=lowfreq_vars,
        )

    with out_path.open("w") as out:
        for srow in scaffold:
            key = (srow.pos, srow.insindex)
            call = call_map.get(key)

            if call is None:
                fields = DEFAULT_CALL.copy()
            else:
                allele = call.allele

                # Perl behaviour:
                #   $allel1 = lc($allel1) unless(exists $strains->{$pos.$index.$id});
                if (srow.pos, srow.insindex, sample_id) not in original_variant_keys:
                    allele = allele.lower()

                subst = call.subst
                # Match Perl print_joint_table behaviour:
                #   my $subs = $values[8];
                #   $subs = "-" unless($subs);
                #
                # In Perl, a single space " " is truthy, so it is preserved.
                # Only an actually empty/missing value should become "-".
                if subst is None or subst == "":
                    subst = "-"

                fields = [
                    call.type or "-",
                    allele or "-",
                    str(call.cov_f),
                    str(call.cov_r),
                    str(call.qual20),
                    format_freq_for_join(call.freq),
                    str(call.cov),
                    subst,
                ]

            out.write("\t".join(fields) + "\n")

    return out_path


# =============================================================================
# Join table writing
# =============================================================================

def write_group_files(
    scaffold: List[ScaffoldRow],
    ids: List[str],
    sample_calls_dir: Path,
    group_dir: Path,
    group_size: int,
    ) -> List[Path]:
    """
    Build intermediate grouped wide call files.

    This avoids opening thousands of per-sample files simultaneously when
    writing the final Join table.
    """
    group_paths: List[Path] = []

    for group_idx, start in enumerate(range(0, len(ids), group_size), start=1):
        group_ids = ids[start : start + group_size]
        group_path = group_dir / f"group_{group_idx:04d}.wide.tsv"
        group_paths.append(group_path)

        handles = []

        try:
            for sid in group_ids:
                handles.append(
                    (sid, (sample_calls_dir / f"{sid}.calls.tsv").open())
                )

            with group_path.open("w") as out:
                for _ in scaffold:
                    parts: List[str] = []

                    for _sid, handle in handles:
                        line = handle.readline()

                        if not line:
                            parts.extend(DEFAULT_CALL)
                        else:
                            fields = line.rstrip("\r\n").split("\t")

                            if len(fields) < 8:
                                fields += DEFAULT_CALL[len(fields):]

                            parts.extend(fields[:8])

                    out.write("\t".join(parts) + "\n")

        finally:
            for _sid, handle in handles:
                handle.close()

    return group_paths


def write_final_join_table(
    join_path: Path,
    scaffold: List[ScaffoldRow],
    ids: List[str],
    group_paths: List[Path],
) -> None:
    """
    Write final MTBSeq-compatible Join table.

    Header intentionally mimics TBtools::print_joint_table_scaffold:

        print OUT "#\\t\\t\\t\\t\\t\\t$header\\n";
        print OUT "#Position\\tInsindex\\tRef\\tGene\\tGeneName\\tAnnotation\\t$second_header\\n";
    """
    first_header = "#\t\t\t\t\t\t" + "\t\t\t\t\t\t\t\t".join(ids)

    second_header = "#Position\tInsindex\tRef\tGene\tGeneName\tAnnotation\t"
    for _sid in ids:
        second_header += "Type\tAllel\tCovFor\tCovRev\tQual20\tFreq\tCov\tSubst\t"

    group_handles = []

    try:
        for group_path in group_paths:
            group_handles.append(group_path.open())

        with join_path.open("w") as out:
            out.write(first_header + "\n")
            out.write(second_header + "\n")

            for srow in scaffold:
                base = [
                    str(srow.pos),
                    str(srow.insindex),
                    srow.ref,
                    srow.gene,
                    srow.gene_name,
                    srow.product,
                ]

                call_parts: List[str] = []

                for handle in group_handles:
                    line = handle.readline()

                    if line:
                        call_parts.extend(line.rstrip("\r\n").split("\t"))

                out.write("\t".join(base + call_parts) + "\n")

    finally:
        for handle in group_handles:
            handle.close()


# =============================================================================
# Breadth log writing
# =============================================================================

def write_breadth_log(
    path: Path,
    genome_length: int,
    ids: List[str],
    any_cov: np.ndarray,
    unambig: np.ndarray,
    nothing: np.ndarray,
) -> None:
    """
    Reproduce TBtools::print_position_stats output.

    Misspellings are intentionally preserved:
        unambigous
        AnyCoverager
    """
    n = len(ids)

    n95 = int(n / 100 * 95)
    n90 = int(n / 100 * 90)
    n75 = int(n / 100 * 75)

    any_slice = any_cov[1 : genome_length + 1]
    unambig_slice = unambig[1 : genome_length + 1]
    nothing_slice = nothing[1 : genome_length + 1]

    positions_covered_100 = int(np.sum(any_slice == n))
    positions_covered_95 = int(np.sum(any_slice >= n95))
    positions_covered_90 = int(np.sum(any_slice >= n90))
    positions_covered_75 = int(np.sum(any_slice >= n75))

    positions_unambig_100 = int(np.sum(unambig_slice == n))
    positions_unambig_95 = int(np.sum(unambig_slice >= n95))
    positions_unambig_90 = int(np.sum(unambig_slice >= n90))
    positions_unambig_75 = int(np.sum(unambig_slice >= n75))

    positions_uncovered = int(np.sum(nothing_slice == n))

    with path.open("w") as out:
        out.write(f"# Number of datasets:\t{n}\n")
        out.write(f"# Genome length:\t{genome_length}\n")
        out.write(f"# Positions covered in 100% of datasets:\t{positions_covered_100}\n")
        out.write(f"# Positions covered in 95% of datasets:\t{positions_covered_95}\n")
        out.write(f"# Positions covered in 90% of datasets:\t{positions_covered_90}\n")
        out.write(f"# Positions covered in 75% of datasets:\t{positions_covered_75}\n")
        out.write(f"# Positions unambigous in 100% of datasets:\t{positions_unambig_100}\n")
        out.write(f"# Positions unambigous in 95% of datasets:\t{positions_unambig_95}\n")
        out.write(f"# Positions unambigous in 90% of datasets:\t{positions_unambig_90}\n")
        out.write(f"# Positions unambigous in 75% of datasets:\t{positions_unambig_75}\n")
        out.write(f"# Positions uncovered in all datasets:\t{positions_uncovered}\n")
        out.write("# IDS: \n")
        out.write("\n".join(ids) + "\n")
        out.write("\n")
        out.write("# Position\tInsindex\tUnambigousCoverage\tAnyCoverager\tNoCoverage\n")

        for pos in range(1, genome_length + 1):
            out.write(
                f"{pos}\t0\t{int(unambig[pos])}\t{int(any_cov[pos])}\t{int(nothing[pos])}\n"
            )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Polars-based replacement for MTBSeq TBjoin."
    )

    parser.add_argument("--var-dir", required=True, type=Path)
    parser.add_argument("--pos-out", required=True, type=Path)
    parser.add_argument("--call-out", required=True, type=Path)
    parser.add_argument("--join-out", required=True, type=Path)
    parser.add_argument("--group-name", required=True)

    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA filename inside --var-dir",
    )

    parser.add_argument(
        "--refg",
        required=True,
        help="Annotation filename inside --var-dir",
    )

    parser.add_argument("--mincovf", type=int, default=4)
    parser.add_argument("--mincovr", type=int, default=4)
    parser.add_argument("--minphred20", type=int, default=4)
    parser.add_argument("--minfreq", type=float, default=75)

    parser.add_argument("--snp-vars", type=int, default=0)
    parser.add_argument("--lowfreq-vars", type=int, default=0)

    parser.add_argument(
        "--variant-glob",
        default="*.gatk_position_variants_*.tab",
        help="Glob pattern inside --call-out for variant files.",
    )

    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Intermediate directory. Default: JOIN_OUT/tbjoin_polars_work",
    )

    parser.add_argument(
        "--group-size",
        type=int,
        default=100,
        help="Number of samples per intermediate wide group file.",
    )

    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep intermediate files.",
    )

    args = parser.parse_args()

    args.join_out.mkdir(parents=True, exist_ok=True)

    work_dir = args.work_dir or (args.join_out / "tbjoin_polars_work")
    sample_calls_dir = work_dir / "sample_calls"
    group_dir = work_dir / "groups"

    if work_dir.exists():
        shutil.rmtree(work_dir)

    sample_calls_dir.mkdir(parents=True, exist_ok=True)
    group_dir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Parsing reference FASTA...")
    ref_genome = parse_fasta(args.var_dir / args.ref)
    genome_length = len(ref_genome)
    print(f"[INFO] Reference length: {genome_length}")

    print("[INFO] Parsing annotation...")
    genes, annotation_by_pos = parse_annotation(args.var_dir / args.refg)
    print(f"[INFO] Parsed genes: {len(genes)}")

    print("[INFO] Finding variant files...")
    var_files = sorted(args.call_out.glob(args.variant_glob), key=lambda p: p.name)
    print(f"[INFO] Variant files: {len(var_files)}")

    print("[INFO] Building global variant index...")
    variant_index, original_variant_keys, ids = parse_variant_files(
        call_out=args.call_out,
        var_files=var_files,
        micovf=args.mincovf,
        micovr=args.mincovr,
        mifreq=args.minfreq,
        miphred20=args.minphred20,
        snp_vars=args.snp_vars,
    )

    print(f"[INFO] Samples: {len(ids)}")
    print(f"[INFO] Global variant rows: {variant_index.height}")

    join_file, breadth_file = build_output_names(
        group_name=args.group_name,
        micovf=args.mincovf,
        micovr=args.mincovr,
        mifreq=args.minfreq,
        miphred20=args.minphred20,
        nsamples=len(ids),
    )

    join_path = args.join_out / join_file
    breadth_path = args.join_out / breadth_file

    print("[INFO] Building scaffold...")
    scaffold = build_scaffold(
        variant_index=variant_index,
        genes=genes,
        annotation_by_pos=annotation_by_pos,
    )

    required_positions = build_required_positions(variant_index)

    print(f"[INFO] Scaffold rows: {len(scaffold)}")
    print(f"[INFO] Required internal rows: {required_positions.height}")

    any_cov = np.zeros(genome_length + 1, dtype=np.uint16)
    unambig = np.zeros(genome_length + 1, dtype=np.uint16)
    nothing = np.zeros(genome_length + 1, dtype=np.uint16)

    print("[INFO] Processing samples...")

    for i, sample_id in enumerate(ids, start=1):
        print(f"[INFO] [{i}/{len(ids)}] {sample_id}")

        process_sample(
            sample_id=sample_id,
            pos_out=args.pos_out,
            sample_calls_dir=sample_calls_dir,
            scaffold=scaffold,
            required_positions=required_positions,
            original_variant_keys=original_variant_keys,
            genes=genes,
            annotation_by_pos=annotation_by_pos,
            ref_genome=ref_genome,
            any_cov=any_cov,
            unambig=unambig,
            nothing=nothing,
            micovf=args.mincovf,
            micovr=args.mincovr,
            miphred20=args.minphred20,
            mifreq=args.minfreq,
            lowfreq_vars=args.lowfreq_vars,
        )

    print("[INFO] Writing grouped intermediate wide files...")
    group_paths = write_group_files(
        scaffold=scaffold,
        ids=ids,
        sample_calls_dir=sample_calls_dir,
        group_dir=group_dir,
        group_size=args.group_size,
    )

    print("[INFO] Writing final MTBSeq-compatible Join table...")
    write_final_join_table(
        join_path=join_path,
        scaffold=scaffold,
        ids=ids,
        group_paths=group_paths,
    )

    print("[INFO] Writing breadth log...")
    write_breadth_log(
        path=breadth_path,
        genome_length=genome_length,
        ids=ids,
        any_cov=any_cov,
        unambig=unambig,
        nothing=nothing,
    )

    if not args.keep_work:
        shutil.rmtree(work_dir)

    print("[INFO] Done.")
    print(f"[INFO] Join table: {join_path}")
    print(f"[INFO] Breadth log: {breadth_path}")


if __name__ == "__main__":
    main()