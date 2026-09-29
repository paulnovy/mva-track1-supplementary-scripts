#!/usr/bin/env python3
"""Extract a small named-gene panel from a bgzip/gzip VCF without altering calls."""
import os
import argparse
import csv
import gzip
from pathlib import Path
from classify_canonical_panel_rows import changed_span, overlaps


def alt_present(gt: str) -> bool:
    alleles = gt.replace("|", "/").split("/")
    return any(a not in {"0", "."} for a in alleles)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path,
                        help="TSV: symbol, ENSG, chromosome, start, end, strand, description")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    intervals = []
    with args.coordinates.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row:
                intervals.append((row[0], row[2].removeprefix("chr"), int(row[3]), int(row[4])))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    header = ["gene", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample", "gt", "alt_present"]
    with gzip.open(args.vcf, "rt") as source, args.out.open("x", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0].removeprefix("chr")
            pos = int(fields[1])
            if len(fields) != 10:
                raise ValueError("Use a single-sample VCF")
            spans = [changed_span(pos, fields[3], alt) for alt in fields[4].split(",")]
            for gene, gene_chrom, start, end in intervals:
                if chrom == gene_chrom and any(overlaps(span, (start, end)) for span in spans):
                    format_keys = fields[8].split(":")
                    sample_values = fields[9].split(":")
                    gt = sample_values[format_keys.index("GT")] if "GT" in format_keys and format_keys.index("GT") < len(sample_values) else "."
                    writer.writerow([gene, *fields[:10], gt, str(alt_present(gt)).lower()])


if __name__ == "__main__":
    main()
