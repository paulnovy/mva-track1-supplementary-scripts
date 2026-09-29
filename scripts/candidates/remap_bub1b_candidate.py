#!/usr/bin/env python3
"""Build and score reference/alternate local haplotypes for one exact insertion."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("hidden", HERE / "analyze_bub1b_hidden_alleles.py")
hidden = importlib.util.module_from_spec(spec); spec.loader.exec_module(hidden)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + "\n")


def parse_paf(path: Path) -> dict[str, int]:
    scores = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {field[:2]: field[5:] for field in fields[12:] if len(field) > 5 and field[2:5] == ":i:"}
            score = int(tags.get("AS", fields[9]))
            scores[fields[0]] = max(score, scores.get(fields[0], -10**9))
    return scores


def prepare(args: argparse.Namespace) -> None:
    args.output.mkdir(parents=True, exist_ok=True)
    rows = hidden.read_sam(args.sam, "chr15", args.region_start, args.region_end)
    selected_names = {row["qhash"] for row in rows if row["rname"] == "chr15" and not (row["flag"] & 4)
                      and row["pos"] <= args.position + args.radius and row["end"] >= args.position - args.radius}
    best = {}
    for row in rows:
        if row["qhash"] not in selected_names or row["flag"] & (4 | 256 | 512 | 1024 | 2048) or row["sequence"] == "*":
            continue
        end_name = "1" if row["flag"] & 64 else ("2" if row["flag"] & 128 else "0")
        key = f"{row['qhash']}/{end_name}"
        rank = (row["mapq"], row["aligned_query"])
        if key not in best or rank > best[key][0]:
            best[key] = (rank, row["sequence"])
    reference_start, reference = hidden.load_fasta(args.reference)
    window_start = args.position - args.radius
    window_end = args.position + args.radius
    ref_haplotype = hidden.ref_slice(reference_start, reference, window_start, window_end)
    insertion_offset = args.position - window_start + 1
    alt_haplotype = ref_haplotype[:insertion_offset] + args.insertion.upper() + ref_haplotype[insertion_offset:]
    write_fasta(args.output / "reference-haplotype.fa", [(f"ref_chr15_{window_start}_{window_end}", ref_haplotype)])
    write_fasta(args.output / "alternate-haplotype.fa", [(f"alt_chr15_{args.position}_ins{len(args.insertion)}", alt_haplotype)])
    write_fasta(args.output / "haplotypes.fa", [(f"ref_chr15_{window_start}_{window_end}", ref_haplotype),
                                                 (f"alt_chr15_{args.position}_ins{len(args.insertion)}", alt_haplotype)])
    write_fasta(args.output / "reads.fa", [(name, value[1]) for name, value in sorted(best.items())])
    write_fasta(args.output / "insertion.fa", [(f"chr15_{args.position}_inserted_sequence", args.insertion.upper())])
    metadata = {"status": "prepared", "position_1based_left_flank": args.position,
                "allele_boundary": f"chr15:{args.position}|{args.position + 1}",
                "insertion_length": len(args.insertion), "insertion_sequence": args.insertion.upper(),
                "window_1based_inclusive": [window_start, window_end], "read_ends": len(best),
                "independent_templates": len(selected_names),
                "read_identifiers": "SHA-256(template name)/read-end only"}
    (args.output / "preparation.json").write_text(json.dumps(metadata, indent=2) + "\n")


def summarize(args: argparse.Namespace) -> None:
    ref = parse_paf(args.output / "reads-to-reference.paf")
    alt = parse_paf(args.output / "reads-to-alternate.paf")
    read_rows = []
    for name in sorted(set(ref) | set(alt)):
        delta = alt.get(name, -10**9) - ref.get(name, -10**9)
        read_rows.append({"read_end_hash": name, "reference_as": ref.get(name), "alternate_as": alt.get(name),
                          "alternate_minus_reference_as": delta})
    by_template = defaultdict(list)
    for row in read_rows:
        by_template[row["read_end_hash"].split("/")[0]].append(row["alternate_minus_reference_as"])
    template_delta = {name: sum(delta for delta in values if abs(delta) < 10**8) for name, values in by_template.items()}
    counts = {
        "alt_better_by_at_least_10": sum(value >= 10 for value in template_delta.values()),
        "ref_better_by_at_least_10": sum(value <= -10 for value in template_delta.values()),
        "difference_within_9": sum(abs(value) < 10 for value in template_delta.values()),
    }
    donor_lines = sum(1 for line in (args.output / "insertion-to-grch38.paf").open() if line.strip())
    payload = json.loads((args.output / "preparation.json").read_text())
    payload.update({"status": "complete", "aligner": "minimap2 2.27-r1193",
                    "template_score_comparison": counts, "templates_scored": len(template_delta),
                    "inserted_sequence_grch38_paf_alignments": donor_lines,
                    "interpretation": "Alternate-preferential local alignments support the represented insertion allele; ties are expected for read ends that do not cross the breakpoint. The inserted low-complexity sequence alone may have no reportable minimap2 hit even when short motifs recur in GRCh38.",
                    "outputs_sha256": {path.name: hidden.sha256(path) for path in args.output.iterdir()
                                       if path.is_file() and path.name != "remap-summary.json"}})
    (args.output / "remap-summary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--sam", required=True, type=Path); prep.add_argument("--reference", required=True, type=Path)
    prep.add_argument("--region-start", required=True, type=int); prep.add_argument("--region-end", required=True, type=int)
    prep.add_argument("--position", required=True, type=int); prep.add_argument("--insertion", required=True)
    prep.add_argument("--radius", type=int, default=500); prep.add_argument("--output", required=True, type=Path)
    finish = sub.add_parser("summarize"); finish.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args) if args.command == "prepare" else summarize(args)


if __name__ == "__main__":
    main()
