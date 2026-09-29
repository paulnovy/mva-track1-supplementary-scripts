#!/usr/bin/env python3
"""Exact algebraic bulk RD/BAF identifiability examples (no patient reanalysis)."""

from __future__ import annotations

import os
import argparse
import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Mixture:
    name: str
    fractions: dict[str, float]


COPIES = {
    "AB": (1, 1),
    "AAB": (2, 1),
    "ABB": (1, 2),
    "A": (1, 0),
    "B": (0, 1),
}


def totals(fractions: dict[str, float]) -> tuple[float, float, float, float]:
    if abs(sum(fractions.values()) - 1.0) > 1e-12:
        raise ValueError("cell fractions must sum to one")
    a = sum(fractions.get(state, 0.0) * COPIES[state][0] for state in COPIES)
    b = sum(fractions.get(state, 0.0) * COPIES[state][1] for state in COPIES)
    cn = a + b
    return a, b, cn / 2.0, b / cn


def coherent_gain(f: float) -> Mixture:
    return Mixture(f"coherent gain, f={f:g}", {"AB": 1 - f, "AAB": f})


def coherent_loss(f: float) -> Mixture:
    return Mixture(f"coherent loss, f={f:g}", {"AB": 1 - f, "A": f})


def enumerate_cells(counts: dict[str, int]) -> tuple[float, float, float, float]:
    n = sum(counts.values())
    if n <= 0:
        raise ValueError("at least one cell is required")
    # Deliberately sum integer chromosome copies rather than delegating to the
    # fractional-mixture formula. This is an independent arithmetic check.
    a_total = sum(counts.get(state, 0) * COPIES[state][0] for state in COPIES)
    b_total = sum(counts.get(state, 0) * COPIES[state][1] for state in COPIES)
    return a_total / n, b_total / n, (a_total + b_total) / (2 * n), b_total / (a_total + b_total)


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins", type=Path, help="Optional local bin-file provenance only; no mosaicism calibration")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    args.out.mkdir(mode=0o700, parents=True, exist_ok=False)

    f_values = (0.0, 0.05, 0.10, 0.20, 0.40)
    mixtures = [coherent_gain(f) for f in f_values] + [coherent_loss(f) for f in f_values]
    mixtures += [
        Mixture(
            "equal AAB/ABB/A/B mixture, f=0.40",
            {"AB": 0.60, "AAB": 0.10, "ABB": 0.10, "A": 0.10, "B": 0.10},
        ),
        Mixture(
            "near-cancellation, g=0.11 l=0.09 balanced haplotypes",
            {"AB": 0.80, "AAB": 0.055, "ABB": 0.055, "A": 0.045, "B": 0.045},
        ),
        Mixture(
            "near-cancellation, g=0.11 l=0.09 monosomy retaining A only",
            {"AB": 0.80, "AAB": 0.055, "ABB": 0.055, "A": 0.09},
        ),
        Mixture(
            "imbalanced, g=0.11 AAB and l=0.09 A",
            {"AB": 0.80, "AAB": 0.11, "A": 0.09},
        ),
    ]
    rows = []
    for mixture in mixtures:
        a, b, rd, baf = totals(mixture.fractions)
        rows.append({
            "scenario": mixture.name,
            "AB": mixture.fractions.get("AB", 0), "AAB": mixture.fractions.get("AAB", 0),
            "ABB": mixture.fractions.get("ABB", 0), "A": mixture.fractions.get("A", 0),
            "B": mixture.fractions.get("B", 0), "aneuploid_cell_fraction": 1 - mixture.fractions.get("AB", 0),
            "expected_A_copies_per_cell": a, "expected_B_copies_per_cell": b,
            "expected_relative_RD": rd, "expected_BAF_B": baf,
        })
    cols = list(rows[0])
    write_tsv(args.out / "bulk_mixture_expectations.tsv", rows, cols)

    enumerated = [
        ("balanced f=0.40 (N=400)", {"AB": 240, "AAB": 40, "ABB": 40, "A": 40, "B": 40}),
        ("balanced near-cancellation g=0.11 l=0.09 (N=200)", {"AB": 160, "AAB": 11, "ABB": 11, "A": 9, "B": 9}),
    ]
    check_rows = []
    for name, counts in enumerated:
        a, b, rd, baf = enumerate_cells(counts)
        check_rows.append({"scenario": name, "counts": ";".join(f"{k}={v}" for k, v in counts.items()),
                           "A_total_per_cell": a, "B_total_per_cell": b,
                           "relative_RD": rd, "BAF_B": baf})
    write_tsv(args.out / "fixed_cell_enumeration_check.tsv", check_rows, list(check_rows[0]))

    if args.bins is not None:
        with args.bins.open("rb") as handle:
            bins_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
        with args.bins.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            clean = [row for row in reader if row["clean"] == "1"]
        chr_count = len({row["chrom"] for row in clean})
        context = [
            {"input": str(args.bins), "sha256": bins_sha256, "clean_1mb_autosomal_bins": len(clean),
             "autosomes_observed": chr_count,
             "use": "file-structure context only; no RD/BAF values were inspected, shifted, or used as a threshold"}
        ]
        write_tsv(args.out / "input_provenance.tsv", context, list(context[0]))

    report = """# Bulk-mixture identifiability demonstration

## Question and answer

This is an exact cell-mixture algebra demonstration, not a patient-specific mosaicism
analysis. It shows that bulk chromosome-wide read depth (RD) and unsigned/signed B-allele
frequency (BAF) marginals can be identical to diploid AB even when a nonzero fraction of
cells is aneuploid.

For a diploid AB population mixed with equal fractions `f/4` each of AAB, ABB, A, and B,
the expected totals per cell are A=1 and B=1: RD=1 and BAF(B)=0.5 for every f. Thus at
f=0.40, 40% of cells are aneuploid while these bulk marginals are exactly diploid. The
fixed-cell N=400 enumeration independently gives the same values.

## Coherent comparison

For a coherent trisomy-like AAB fraction f, expected RD is `1+f/2` and BAF(B) is
`1/(2+f)`. For a coherent monosomy-like A fraction f, RD is `1-f/2` and BAF(B) is
`(1-f)/(2-f)`. The accompanying table lists f=0, 0.05, 0.10, 0.20, and 0.40.

Balanced gains and losses can cancel. With g=0.11 gain and l=0.09 loss, each split equally
between homologues, A=B=1.01, RD=1.01, BAF(B)=0.5, despite 20% abnormal cells. The
N=200 exact enumeration verifies that result. The imbalanced scenarios show that the
residual RD/BAF pattern depends on which chromosome copy/haplotype is affected.

## Scope and limits

An optional bins input records file-structure counts only in input_provenance.tsv.
No per-bin RD/BAF/noise value was examined, shifted,
or used as a threshold. No BAM, BCF, VCF, API, caller, threshold, segment, or previous
3–5 Mb injection calibration was used or rerun.

This is **not evidence that this patient has any such mixture**, is not a realistic estimate
of aneuploid-cell burden, and is not a sensitivity/limit-of-detection benchmark. It assumes
idealized exact cell states and equal sequencing representation; it omits allele-specific
mapping, read-count sampling, phasing error, GC/mappability, CNVs, lineage structure,
tissue composition, and biology of MVA. It demonstrates nonidentifiability of these bulk
marginals under stated mixtures only. It cannot establish MVA, infer causal order, or show
that aneuploidy is primary rather than secondary.

The two TSV tables are the numerical deliverable. Use plot_bulk_mixture_examples.py
separately to render the synthetic figure when matplotlib is installed.
"""
    (args.out / "BULK_MIXTURE_IDENTIFIABILITY.md").write_text(report)
    for path in args.out.iterdir():
        path.chmod(0o600)


if __name__ == "__main__":
    main()
