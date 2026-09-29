#!/usr/bin/env python3
"""Classify called small-variant alleles against primary and secondary transcripts.

The input transcript JSON is Ensembl ``lookup/id?expand=1`` output.  This
script deliberately preserves the original per-source VCF evidence while
classifying the minimal reference span of each *called* ALT allele.
"""
from __future__ import annotations

import os
import argparse
import csv
import json
from pathlib import Path

import pysam

from classify_canonical_panel_rows import called_alts, changed_span, overlaps


def parse_assignment(value: str) -> tuple[str, Path]:
    try:
        symbol, filename = value.split("=", 1)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected SYMBOL=PATH") from error
    return symbol, Path(filename)


def transcript_geometry(transcript: dict, flank: int) -> tuple[list[tuple[int, int]], list[tuple[str, tuple[int, int]]]]:
    translation = transcript.get("Translation")
    if not translation:
        return [], []
    cds_start, cds_end = sorted((translation["start"], translation["end"]))
    exons = sorted((exon["start"], exon["end"]) for exon in transcript.get("Exon", []))
    coding = [(max(start, cds_start), min(end, cds_end)) for start, end in exons]
    coding = [interval for interval in coding if interval[0] <= interval[1]]
    splice = []
    for (_, left_end), (right_start, _) in zip(exons, exons[1:]):
        if transcript["strand"] == 1:
            splice.extend((("donor", (left_end + 1, left_end + flank)), ("acceptor", (right_start - flank, right_start - 1))))
        else:
            splice.extend((("acceptor", (left_end + 1, left_end + flank)), ("donor", (right_start - flank, right_start - 1))))
    return coding, splice


def region(spans: list[tuple[int, int]], geometry: tuple[list[tuple[int, int]], list[tuple[str, tuple[int, int]]]]) -> str:
    coding, splice = geometry
    if any(overlaps(span, interval) for span in spans for interval in coding):
        return "coding"
    roles = sorted({role for span in spans for role, interval in splice if overlaps(span, interval)})
    return "splice_" + "_".join(roles) if roles else "other"


def load_gene(symbol: str, path: Path, flank: int) -> tuple[str, int, int, list[tuple[str, str, bool, tuple]]]:
    data = json.loads(path.read_text())
    transcripts = [entry for entry in data.get("Transcript", []) if entry.get("biotype") == "protein_coding" and entry.get("Translation")]
    primary_ids = {
        entry["id"] for entry in transcripts
        if entry.get("is_mane_select") or any(item.get("type") == "MANE_Select" for item in entry.get("MANE", []))
    }
    if not primary_ids:
        primary_ids = {entry["id"] for entry in transcripts if entry.get("is_canonical")}
    if not primary_ids:
        raise ValueError(f"{symbol}: no MANE/canonical translated transcript")
    return (
        str(data["seq_region_name"]).removeprefix("chr"),
        data["start"], data["end"],
        [(entry["id"], entry.get("display_name", entry["id"]), entry["id"] in primary_ids, transcript_geometry(entry, flank)) for entry in transcripts],
    )


def sample_value(record: pysam.VariantRecord, sample: str, key: str) -> str:
    value = record.samples[sample].get(key)
    if value is None:
        return "."
    return ",".join("." if item is None else str(item) for item in value) if isinstance(value, tuple) else str(value)


def gt_value(record: pysam.VariantRecord, sample: str) -> str:
    value = record.samples[sample].get("GT")
    if value is None:
        return "."
    separator = "|" if record.samples[sample].phased else "/"
    return separator.join("." if item is None else str(item) for item in value)


def filter_value(record: pysam.VariantRecord) -> str:
    """Preserve VCF's distinct PASS and unfiltered-dot states."""
    values = list(record.filter.keys())
    return ";".join(values) if values else "."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", action="append", type=parse_assignment, required=True, help="label=VCF")
    parser.add_argument("--gene", action="append", type=parse_assignment, required=True, help="symbol=Ensembl expanded JSON")
    parser.add_argument("--splice-flank", type=int, default=8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    if args.splice_flank < 1:
        parser.error("splice flank must be positive")
    genes = {symbol: load_gene(symbol, path, args.splice_flank) for symbol, path in args.gene}
    fields = ["source", "gene", "chrom", "pos", "id", "ref", "alt", "qual", "filter", "format", "sample", "gt", "dp", "ad", "gq", "called_alt_indices", "called_alt", "allele_span", "primary_transcript", "primary_region", "secondary_coding_transcripts", "secondary_splice_flank%d_transcripts" % args.splice_flank]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source, vcf_path in args.source:
            with pysam.VariantFile(str(vcf_path)) as vcf:
                if len(vcf.header.samples) != 1:
                    raise ValueError("Use a single-sample VCF; samples must not be silently pooled")
                sample = next(iter(vcf.header.samples))
                for record in vcf:
                    chrom = record.contig.removeprefix("chr")
                    for symbol, (gene_chrom, start, end, transcripts) in genes.items():
                        if chrom != gene_chrom:
                            continue
                        gt = gt_value(record, sample)
                        indices, alts = called_alts(",".join(record.alts or ()), gt)
                        if not indices:
                            continue
                        spans = [changed_span(record.pos, record.ref, alt) for alt in alts]
                        if not any(overlaps(span, (start, end)) for span in spans):
                            continue
                        primary = [(tid, name, geometry) for tid, name, is_primary, geometry in transcripts if is_primary]
                        secondary = [(tid, name, geometry) for tid, name, is_primary, geometry in transcripts if not is_primary]
                        primary_regions = [f"{tid}:{region(spans, geometry)}" for tid, _, geometry in primary]
                        secondary_coding = [tid for tid, _, geometry in secondary if region(spans, geometry) == "coding"]
                        secondary_splice = [tid for tid, _, geometry in secondary if region(spans, geometry).startswith("splice_")]
                        writer.writerow({
                            "source": source, "gene": symbol, "chrom": record.contig, "pos": record.pos, "id": record.id or ".",
                            "ref": record.ref, "alt": ",".join(record.alts or ()), "qual": record.qual if record.qual is not None else ".",
                            "filter": filter_value(record), "format": ":".join(record.format.keys()), "sample": sample,
                            "gt": gt, "dp": sample_value(record, sample, "DP"), "ad": sample_value(record, sample, "AD"), "gq": sample_value(record, sample, "GQ"),
                            "called_alt_indices": ",".join(map(str, indices)), "called_alt": ",".join(alts), "allele_span": ",".join(f"{a}-{b}" for a, b in spans),
                            "primary_transcript": ",".join(tid for tid, _, _ in primary), "primary_region": ",".join(primary_regions),
                            "secondary_coding_transcripts": ",".join(secondary_coding), "secondary_splice_flank%d_transcripts" % args.splice_flank: ",".join(secondary_splice),
                        })


if __name__ == "__main__":
    main()
