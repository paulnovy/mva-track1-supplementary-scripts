#!/usr/bin/env python3
"""Select one phenotype-linked candidate and cross-review existing CNV segments."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine-result", type=Path, required=True)
    parser.add_argument("--cnv-segments", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    engine = json.loads(args.engine_result.read_text(encoding="utf-8"))
    candidate = next((c for c in engine.get("candidates", [])
                      if c.get("gene_symbol") != "BUB1B"
                      and "HP:0000121" in c.get("contributing_hpo_terms", [])
                      and c.get("variants")), None)
    if not candidate:
        atomic_json(args.output, {
            "schema": "mva1_targeted_existing_cnv_review_v1",
            "status": "not_run_no_concrete_new_nephrocalcinosis_candidate",
            "global_rerun": False,
        })
        return 3
    variants = [v for v in candidate["variants"] if v.get("contig") and v.get("start")]
    if not variants:
        return 3
    contig = str(variants[0]["contig"])
    starts = [int(v["start"]) for v in variants if str(v["contig"]) == contig]
    region = {"contig": contig, "start": max(1, min(starts) - 100_000), "end": max(starts) + 100_000}
    overlaps = []
    with args.cnv_segments.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row.get("chrom") or row.get("chr") or row.get("CHROM")
            try:
                start = int(float(row.get("start") or row.get("START") or ""))
                end = int(float(row.get("end") or row.get("END") or ""))
            except ValueError:
                continue
            if chrom in {contig, contig.removeprefix("chr"), "chr" + contig.removeprefix("chr")} and start <= region["end"] and end >= region["start"]:
                overlaps.append(row)
    atomic_json(args.output, {
        "schema": "mva1_targeted_existing_cnv_review_v1",
        "status": "complete",
        "candidate": {"rank": candidate.get("rank"), "gene_symbol": candidate.get("gene_symbol")},
        "checkpoint": "top_non_BUB1B_candidate_with_nephrocalcinosis_HPO_evidence_and_contributing_variant",
        "region": region,
        "existing_cnv_segment_overlap_count": len(overlaps),
        "existing_cnv_segment_overlaps": overlaps,
        "global_rerun": False,
        "limitations": ["exploratory_existing_RD_BAF_segments_only", "candidate_window_is_not_a_gene_level_CNV_exclusion"],
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
