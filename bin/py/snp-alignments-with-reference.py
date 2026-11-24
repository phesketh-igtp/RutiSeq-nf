#!/usr/bin/env python3
"""
Run SNP-sites on a multi-FASTA alignment to:

1. Create a SNP-only FASTA alignment.
2. Create a VCF containing variable positions.
3. Write constant-site counts using ``snp-sites -C``.
4. Optionally extract the same variable columns from one or more aligned
   reference sequences.
5. Create a combined SNP alignment containing samples and references.

Important:
    The input FASTA and optional reference FASTA must use the same alignment
    coordinate system. Every sequence must have exactly the same aligned length.

Usage
-----
Without a reference:

    python run_snp_sites_and_extract_references.py \
        input_alignment.fasta \
        --prefix results/my_analysis

With one or more aligned references:

    python run_snp_sites_and_extract_references.py \
        input_alignment.fasta \
        references.fasta \
        --prefix results/my_analysis

Outputs
-------
    results/my_analysis.snps.fasta
    results/my_analysis.vcf
    results/my_analysis.constant_sites.txt
    results/my_analysis.combined_snps.fasta

When a reference FASTA is supplied:

    results/my_analysis.reference_snps.fasta

When no reference is supplied, combined_snps.fasta contains the same records
as the SNP-only input alignment.
"""

from __future__ import annotations

import argparse
import bz2
import gzip
import lzma
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import TextIO


def open_text_auto(path: Path) -> TextIO:
    """Open plain, gzip, bzip2, or xz/LZMA text, detected by file magic."""
    with path.open("rb") as raw_handle:
        magic = raw_handle.read(6)

    if magic.startswith(b"\x1f\x8b"):
        return gzip.open(path, "rt", encoding="utf-8")
    if magic.startswith(b"BZh"):
        return bz2.open(path, "rt", encoding="utf-8")
    if magic.startswith(b"\xfd7zXZ\x00"):
        return lzma.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_fasta(path: Path) -> OrderedDict[str, str]:
    """Read a FASTA file and return records in input order."""
    records: OrderedDict[str, str] = OrderedDict()
    name: str | None = None
    chunks: list[str] = []

    with open_text_auto(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()

                name = line[1:].strip().split()[0]
                if not name:
                    raise ValueError(
                        f"Empty FASTA identifier in {path}, line {line_number}"
                    )
                if name in records:
                    raise ValueError(
                        f"Duplicate FASTA identifier {name!r} in {path}"
                    )
                chunks = []
            else:
                if name is None:
                    raise ValueError(
                        f"Sequence before first FASTA header in {path}, "
                        f"line {line_number}"
                    )
                chunks.append("".join(line.split()))

    if name is not None:
        records[name] = "".join(chunks).upper()

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def write_fasta(
    records: OrderedDict[str, str],
    path: Path,
    width: int = 80,
) -> None:
    """Write FASTA records with sequences wrapped to the requested width."""
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            if sequence:
                for start in range(0, len(sequence), width):
                    handle.write(sequence[start:start + width] + "\n")
            else:
                handle.write("\n")


def validate_alignment(
    records: OrderedDict[str, str],
    label: str,
) -> int:
    """Ensure every sequence has the same aligned length."""
    lengths = {len(sequence) for sequence in records.values()}

    if len(lengths) != 1:
        details = ", ".join(
            f"{name}={len(sequence)}"
            for name, sequence in records.items()
        )
        raise ValueError(
            f"{label} sequences do not have equal lengths: {details}"
        )

    return next(iter(lengths))


def run_command(command: list[str]) -> None:
    """Run a command and raise a readable error if it fails."""
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        command_text = " ".join(command)
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {command_text}"
        ) from exc


def run_command_to_file(
    command: list[str],
    output_path: Path,
) -> None:
    """Run a command and write its standard output directly to a file."""
    try:
        with output_path.open("w", encoding="utf-8") as output_handle:
            subprocess.run(
                command,
                check=True,
                stdout=output_handle,
            )
    except subprocess.CalledProcessError as exc:
        command_text = " ".join(command)
        raise RuntimeError(
            f"Command failed with exit code {exc.returncode}: {command_text}"
        ) from exc


def read_vcf_positions(path: Path) -> list[int]:
    """Read one-based alignment positions from a SNP-sites VCF."""
    positions: list[int] = []

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(
                    f"Malformed VCF record in {path}, line {line_number}"
                )

            try:
                position = int(fields[1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid VCF POS in {path}, line {line_number}: "
                    f"{fields[1]!r}"
                ) from exc

            positions.append(position)

    return positions


def unique_reference_name(name: str, existing: set[str]) -> str:
    """Return a unique reference name that does not clash with sample names."""
    if name not in existing:
        return name

    candidate = f"reference__{name}"
    suffix = 2
    while candidate in existing:
        candidate = f"reference__{name}_{suffix}"
        suffix += 1
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Create a SNP-only alignment and VCF with SNP-sites, write "
            "constant-site counts using SNP-sites -C, and optionally extract "
            "the corresponding aligned columns from a reference FASTA."
        )
    )
    parser.add_argument(
        "alignment",
        type=Path,
        help=(
            "Input aligned multi-FASTA; may be plain text, gzip, bzip2, "
            "or xz/LZMA compressed"
        ),
    )
    parser.add_argument(
        "reference",
        nargs="?",
        type=Path,
        default=None,
        help=(
            "Optional aligned reference FASTA; may be plain text, gzip, "
            "bzip2, or xz/LZMA compressed"
        ),
    )
    parser.add_argument(
        "--prefix",
        type=Path,
        default=Path("snp_sites_output"),
        help="Output prefix (default: snp_sites_output)",
    )
    parser.add_argument(
        "--acgt-only",
        action="store_true",
        help=(
            "Pass -c to SNP-sites, retaining only variable columns that "
            "contain exclusively A, C, G, and T"
        ),
    )
    args = parser.parse_args()

    if shutil.which("snp-sites") is None:
        raise RuntimeError(
            "snp-sites was not found in PATH. Install it before running this script."
        )
    if not args.alignment.is_file():
        raise FileNotFoundError(f"Alignment file not found: {args.alignment}")
    if args.reference is not None and not args.reference.is_file():
        raise FileNotFoundError(f"Reference file not found: {args.reference}")

    input_records = parse_fasta(args.alignment)
    input_length = validate_alignment(input_records, "Input alignment")

    reference_records: OrderedDict[str, str] | None = None
    if args.reference is not None:
        reference_records = parse_fasta(args.reference)
        reference_length = validate_alignment(
            reference_records,
            "Reference FASTA",
        )
        if input_length != reference_length:
            raise ValueError(
                "Input and reference sequences must share the same aligned "
                "coordinate system: "
                f"input length={input_length}, "
                f"reference length={reference_length}"
            )

    prefix = args.prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)

    snp_alignment_path = Path(f"{prefix}.snps.fasta")
    vcf_path = Path(f"{prefix}.vcf")
    constant_sites_path = Path(f"{prefix}.constant_sites.txt")
    reference_snps_path = Path(f"{prefix}.reference_snps.fasta")
    combined_path = Path(f"{prefix}.combined_snps.fasta")

    common_options = ["-c"] if args.acgt_only else []

    # Normalize the input to a temporary uncompressed FASTA so each SNP-sites
    # invocation receives exactly the same content.
    with tempfile.TemporaryDirectory(prefix="snp_sites_") as temp_dir:
        normalized_alignment = Path(temp_dir) / "alignment.fasta"
        write_fasta(input_records, normalized_alignment)

        run_command([
            "snp-sites",
            *common_options,
            "-m",
            "-o",
            str(snp_alignment_path),
            str(normalized_alignment),
        ])

        run_command([
            "snp-sites",
            *common_options,
            "-v",
            "-o",
            str(vcf_path),
            str(normalized_alignment),
        ])

        # SNP-sites -C writes constant-site counts to standard output.
        run_command_to_file(
            [
                "snp-sites",
                *common_options,
                "-C",
                str(normalized_alignment),
            ],
            constant_sites_path,
        )

    # VCF positions are needed to validate the outputs and extract equivalent
    # columns from the optional aligned reference sequences.
    positions = read_vcf_positions(vcf_path)

    if len(positions) != len(set(positions)):
        raise RuntimeError(
            "The SNP-sites VCF contains duplicate positions. "
            "Reference-column extraction would be ambiguous."
        )

    for position in positions:
        if position < 1 or position > input_length:
            raise ValueError(
                f"VCF position {position} is outside the alignment length "
                f"{input_length}"
            )

    snp_records = parse_fasta(snp_alignment_path)
    snp_alignment_length = validate_alignment(
        snp_records,
        "SNP-only alignment",
    )
    expected_snp_length = len(positions)

    if snp_alignment_length != expected_snp_length:
        raise RuntimeError(
            "SNP FASTA length does not match the number of VCF records: "
            f"SNP FASTA length={snp_alignment_length}, "
            f"VCF records={expected_snp_length}. "
            "The outputs may have been generated with different filtering behaviour."
        )

    combined_records: OrderedDict[str, str] = OrderedDict(snp_records)

    if reference_records is not None:
        extracted_refs: OrderedDict[str, str] = OrderedDict()
        used_names = set(snp_records)

        for name, sequence in reference_records.items():
            output_name = unique_reference_name(name, used_names)
            used_names.add(output_name)
            extracted_refs[output_name] = "".join(
                sequence[position - 1]
                for position in positions
            )

        write_fasta(extracted_refs, reference_snps_path)
        combined_records.update(extracted_refs)

    write_fasta(combined_records, combined_path)

    print(f"SNP-only input alignment: {snp_alignment_path}")
    print(f"VCF:                      {vcf_path}")
    print(f"Constant-site counts:     {constant_sites_path}")
    if reference_records is not None:
        print(f"Reference SNP alignment:  {reference_snps_path}")
    print(f"Combined SNP alignment:   {combined_path}")
    print(f"Input alignment length:   {input_length}")
    print(f"Variable sites extracted: {len(positions)}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)