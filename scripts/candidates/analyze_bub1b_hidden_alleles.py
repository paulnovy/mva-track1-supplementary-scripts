#!/usr/bin/env python3
"""Summarize locus-wide short-read evidence for hidden BUB1B alleles.

The input SAM is deliberately regional but was created with samtools view -P so
that mapped and unmapped mates outside the interval are retained. Read names are
hashed in every published summary; raw sequence-bearing files remain private.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
TARGETS = {40_209_701: ("T", "G", "p.Leu737Ter"), 40_220_612: ("T", "G", "p.Asn1002Lys")}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def name_hash(name: str) -> str:
    return hashlib.sha256(name.encode()).hexdigest()


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def entropy(sequence: str) -> float:
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    counts = Counter(sequence)
    return round(-sum((n / len(sequence)) * math.log2(n / len(sequence)) for n in counts.values()), 4)


def longest_exact(query: str, target: str, floor: int = 10, ceiling: int = 40) -> int:
    query = query.upper(); target = target.upper()
    for width in range(min(len(query), len(target), ceiling), floor - 1, -1):
        if any(query[i:i + width] in target for i in range(len(query) - width + 1)):
            return width
    return 0


def mean_bq(quality: str) -> float | None:
    if not quality or quality == "*":
        return None
    return round(statistics.fmean(ord(char) - 33 for char in quality), 2)


def parse_alignment(fields: list[str], chrom: str, start: int, end: int) -> dict:
    flag = int(fields[1]); pos = int(fields[3]); cigar = fields[5]
    sequence = fields[9]; quality = fields[10]
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    ref_cursor = pos; query_cursor = 0; events = []
    aligned_query = 0
    for index, (length, op) in enumerate(operations):
        if op in "M=X":
            ref_cursor += length; query_cursor += length; aligned_query += length
        elif op in "DN":
            ref_cursor += length
        elif op == "I":
            inserted = sequence[query_cursor:query_cursor + length]
            inserted_q = quality[query_cursor:query_cursor + length] if quality != "*" else ""
            events.append({"kind": "insertion", "left": ref_cursor - 1, "right": ref_cursor,
                           "length": length, "sequence": inserted, "mean_bq": mean_bq(inserted_q)})
            query_cursor += length; aligned_query += length
        elif op == "S":
            clipped = sequence[query_cursor:query_cursor + length]
            clipped_q = quality[query_cursor:query_cursor + length] if quality != "*" else ""
            if index == 0:
                left, right, side = pos - 1, pos, "leading"
            else:
                left, right, side = ref_cursor - 1, ref_cursor, "trailing"
            events.append({"kind": "softclip", "side": side, "left": left, "right": right,
                           "length": length, "sequence": clipped, "mean_bq": mean_bq(clipped_q)})
            query_cursor += length
        elif op == "H" or op == "P":
            pass
    ref_end = ref_cursor - 1 if not (flag & 4) and cigar != "*" else 0
    tags = {}
    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return {
        "qname": fields[0], "qhash": name_hash(fields[0]), "flag": flag, "rname": fields[2],
        "pos": pos, "end": ref_end, "mapq": int(fields[4]), "cigar": cigar,
        "mate_rname": fields[6], "mate_pos": int(fields[7]), "tlen": int(fields[8]),
        "sequence": sequence, "quality": quality, "events": events, "aligned_query": aligned_query,
        "sa": tags.get("SA", ""), "reverse": bool(flag & 16),
        "on_target": fields[2] == chrom and pos <= end and ref_end >= start and not (flag & 4),
    }


def read_sam(path: Path, chrom: str, start: int, end: int) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 11:
                rows.append(parse_alignment(fields, chrom, start, end))
    return rows


def eligible(row: dict, mapq: int = 20) -> bool:
    return row["on_target"] and row["mapq"] >= mapq and row["aligned_query"] >= 30 and not (
        row["flag"] & (4 | 256 | 512 | 1024 | 2048)
    )


def load_fasta(path: Path) -> tuple[int, str]:
    with path.open() as handle:
        header = next(handle).strip()
        sequence = "".join(line.strip() for line in handle)
    match = re.search(r":(\d+)-(\d+)$", header)
    if not match:
        raise ValueError(f"Regional FASTA header lacks coordinates: {header}")
    return int(match.group(1)), sequence.upper()


def ref_slice(reference_start: int, reference: str, start: int, end: int) -> str:
    return reference[start - reference_start:end - reference_start + 1]


def load_depth(path: Path) -> dict[int, int]:
    values = {}
    with path.open() as handle:
        for line in handle:
            _, position, depth = line.rstrip("\n").split("\t")
            values[int(position)] = int(depth)
    return values


def compress_positions(positions: list[int]) -> list[dict]:
    if not positions:
        return []
    result = []; left = previous = positions[0]
    for position in positions[1:]:
        if position != previous + 1:
            result.append({"start": left, "end": previous, "bases": previous - left + 1})
            left = position
        previous = position
    result.append({"start": left, "end": previous, "bases": previous - left + 1})
    return result


def depth_summary(q0: dict[int, int], q20: dict[int, int], q60: dict[int, int], start: int, end: int) -> dict:
    values0 = [q0.get(pos, 0) for pos in range(start, end + 1)]
    values20 = [q20.get(pos, 0) for pos in range(start, end + 1)]
    values60 = [q60.get(pos, 0) for pos in range(start, end + 1)]
    low = [pos for pos in range(start, end + 1) if q20.get(pos, 0) < 10]
    mapping_uncertain = [pos for pos in range(start, end + 1)
                         if q0.get(pos, 0) >= 10 and q20.get(pos, 0) / q0[pos] < 0.5]
    return {
        "bases": len(values20),
        "mapq0_baseq20": {"mean": round(statistics.fmean(values0), 4), "median": statistics.median(values0),
                           "minimum": min(values0), "maximum": max(values0)},
        "mapq20_baseq20": {"mean": round(statistics.fmean(values20), 4), "median": statistics.median(values20),
                            "minimum": min(values20), "maximum": max(values20),
                            "bases_ge_10_fraction": round(sum(x >= 10 for x in values20) / len(values20), 6),
                            "bases_ge_20_fraction": round(sum(x >= 20 for x in values20) / len(values20), 6)},
        "mapq60_baseq20": {"mean": round(statistics.fmean(values60), 4), "median": statistics.median(values60),
                            "minimum": min(values60), "maximum": max(values60)},
        "depth_lt_10_intervals": compress_positions(low),
        "mapq20_less_than_half_of_mapq0_intervals": compress_positions(mapping_uncertain),
    }


def build_signal_clusters(rows: list[dict], alu: str, reference_start: int, reference: str) -> list[dict]:
    observations = []
    for row in rows:
        if not eligible(row):
            continue
        for event in row["events"]:
            threshold = 10 if event["kind"] == "insertion" else 15
            if event["length"] < threshold or (event["mean_bq"] is not None and event["mean_bq"] < 20):
                continue
            observations.append({**event, "qhash": row["qhash"], "mapq": row["mapq"],
                                 "reverse": row["reverse"], "start": row["pos"], "sa": bool(row["sa"])})
    observations.sort(key=lambda item: item["left"])
    grouped = []
    for item in observations:
        if not grouped or item["left"] - grouped[-1][-1]["left"] > 5:
            grouped.append([item])
        else:
            grouped[-1].append(item)
    clusters = []
    for items in grouped:
        templates = {item["qhash"] for item in items}
        insertion_templates = {item["qhash"] for item in items if item["kind"] == "insertion"}
        clip_templates = {item["qhash"] for item in items if item["kind"] == "softclip"}
        seq_counts = Counter((item["sequence"].upper(), item["length"]) for item in items if item["kind"] == "insertion")
        leading = {item["qhash"] for item in items if item.get("side") == "leading"}
        trailing = {item["qhash"] for item in items if item.get("side") == "trailing"}
        sequences = [item["sequence"] for item in items]
        left = min(item["left"] for item in items); right = max(item["right"] for item in items)
        repeat_context = ref_slice(reference_start, reference, max(reference_start, left - 50), left + 50)
        top_sequence = None
        if seq_counts:
            (sequence, length), count = seq_counts.most_common(1)[0]
            top_sequence = {"sequence": sequence, "length": length, "observations": count,
                            "independent_templates": len({item["qhash"] for item in items
                                                           if item["kind"] == "insertion" and item["sequence"].upper() == sequence})}
        credible = len(templates) >= 3 and len({item["start"] for item in items}) >= 2 and (
            len(insertion_templates) >= 3 or len(clip_templates) >= 3
        )
        clusters.append({
            "cluster_id": f"chr15:{left}-{right}", "breakpoint_interval_1based": [left, right],
            "boundary_convention": "insertion or clipping lies between 1-based left and right flanking bases",
            "observations": len(items), "independent_templates": len(templates),
            "insertion_templates": len(insertion_templates), "softclip_templates": len(clip_templates),
            "leading_clip_templates": len(leading), "trailing_clip_templates": len(trailing),
            "mapq60_templates": len({item["qhash"] for item in items if item["mapq"] >= 60}),
            "forward_templates": len({item["qhash"] for item in items if not item["reverse"]}),
            "reverse_templates": len({item["qhash"] for item in items if item["reverse"]}),
            "sa_tag_templates": len({item["qhash"] for item in items if item["sa"]}),
            "distinct_alignment_starts": len({item["start"] for item in items}),
            "lengths": dict(sorted(Counter(item["length"] for item in items).items())),
            "top_exact_insertion": top_sequence,
            "sequence_features": {
                "poly_a_or_t_10_templates": len({item["qhash"] for item in items
                                                  if "A" * 10 in item["sequence"].upper() or "T" * 10 in item["sequence"].upper()}),
                "aluya5_exact_ge_15_templates": len({item["qhash"] for item in items
                                                      if max(longest_exact(item["sequence"], alu),
                                                             longest_exact(item["sequence"], revcomp(alu))) >= 15}),
                "median_sequence_entropy_bits": round(statistics.median(entropy(seq) for seq in sequences), 4),
                "reference_101bp_entropy_bits": entropy(repeat_context),
            },
            "passes_candidate_signal_threshold": credible,
        })
    return sorted(clusters, key=lambda item: (-item["independent_templates"], item["breakpoint_interval_1based"][0]))


def discordant_summary(rows: list[dict]) -> dict:
    anchored = [row for row in rows if eligible(row) and row["flag"] & 1]
    categories = defaultdict(set); mate_seen = defaultdict(set); anchor_bins = defaultdict(set)
    by_hash = defaultdict(list)
    for row in rows:
        by_hash[row["qhash"]].append(row)
    for row in anchored:
        if row["flag"] & 8:
            categories["mate_unmapped"].add(row["qhash"])
        elif row["mate_rname"] not in {"=", "chr15"}:
            categories["mate_other_contig"].add(row["qhash"])
        elif abs(row["tlen"]) > 1000:
            categories["absolute_template_length_gt_1000"].add(row["qhash"])
        else:
            continue
        anchor_bins[(row["pos"] // 200) * 200].add(row["qhash"])
        for mate in by_hash[row["qhash"]]:
            if mate is not row and (mate["flag"] & 128 != row["flag"] & 128):
                mate_seen[row["qhash"]].add((mate["rname"], mate["pos"], bool(mate["flag"] & 4)))
    anomalous = set().union(*categories.values()) if categories else set()
    return {
        "categories_independent_templates": {key: len(value) for key, value in sorted(categories.items())},
        "all_anomalous_independent_templates": len(anomalous),
        "anomalous_templates_with_a_retrieved_mate_alignment": sum(bool(mate_seen[key]) for key in anomalous),
        "top_200bp_anchor_bins": [{"start": key, "end": key + 199, "independent_templates": len(value)}
                                  for key, value in sorted(anchor_bins.items(), key=lambda item: (-len(item[1]), item[0]))[:20]],
        "retrieval_note": "samtools view -P performs a second pass to retrieve complete pairs outside the requested region; supplementary alignments outside the region are represented only when fetched with a selected pair or by SA tags.",
    }


def target_base_counts(rows: list[dict]) -> dict:
    result = {}
    for target, (ref, alt, label) in TARGETS.items():
        per_template = {}
        for row in rows:
            if not eligible(row) or not (row["pos"] <= target <= row["end"]):
                continue
            ref_cursor = row["pos"]; query_cursor = 0; base = None; bq = None
            for length, op in CIGAR_RE.findall(row["cigar"]):
                length = int(length)
                if op in "M=X":
                    if ref_cursor <= target < ref_cursor + length:
                        offset = query_cursor + target - ref_cursor
                        base = row["sequence"][offset].upper()
                        bq = ord(row["quality"][offset]) - 33 if row["quality"] != "*" else 99
                        break
                    ref_cursor += length; query_cursor += length
                elif op in "DN": ref_cursor += length
                elif op in "IS": query_cursor += length
            if base and bq >= 20:
                per_template.setdefault(row["qhash"], base)
        counts = Counter(per_template.values())
        result[str(target)] = {"label": label, "ref": ref, "alt": alt,
                               "fragment_counts": dict(sorted(counts.items())), "total": sum(counts.values())}
    return result


def load_delly(path: Path) -> list[dict]:
    rows = []
    if not path.exists():
        return rows
    with path.open() as handle:
        header = next(handle, "").rstrip("\n").split("\t")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            if row.get("pos", "").isdigit():
                rows.append(row)
    return rows


def annotation_for(left: int, gene_start: int, gene_end: int, exons: list[tuple[int, int]]) -> dict:
    if left < gene_start:
        return {"context": "5prime_flank", "distance_to_gene_bases": gene_start - left}
    if left > gene_end:
        return {"context": "3prime_flank", "distance_to_gene_bases": left - gene_end}
    exon = next(((a, b) for a, b in exons if a <= left <= b), None)
    return {"context": "exon" if exon else "intron", "mane_exon_interval": list(exon) if exon else None}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True, type=Path)
    parser.add_argument("--depth-q0", required=True, type=Path)
    parser.add_argument("--depth-q20", required=True, type=Path)
    parser.add_argument("--depth-q60", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--mane", required=True, type=Path)
    parser.add_argument("--alu", required=True, type=Path)
    parser.add_argument("--delly", required=True, type=Path)
    parser.add_argument("--chrom", default="chr15")
    parser.add_argument("--start", required=True, type=int)
    parser.add_argument("--end", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    rows = read_sam(args.sam, args.chrom, args.start, args.end)
    reference_start, reference = load_fasta(args.reference)
    mane = json.loads(args.mane.read_text())
    gene_start, gene_end = int(mane["start"]), int(mane["end"])
    transcript = next(tx for tx in mane["Transcript"] if tx["id"] == "ENST00000287598")
    exons = sorted((int(exon["start"]), int(exon["end"])) for exon in transcript["Exon"])
    alu_payload = json.loads(args.alu.read_text()); alu = alu_payload["consensus_sequence"]
    clusters = build_signal_clusters(rows, alu, reference_start, reference)
    for cluster in clusters:
        cluster["annotation"] = annotation_for(cluster["breakpoint_interval_1based"][0], gene_start, gene_end, exons)
        cluster["matching_delly_ids"] = [row.get("id") for row in load_delly(args.delly)
                                          if row.get("pos", "").isdigit() and
                                          abs(int(row["pos"]) - cluster["breakpoint_interval_1based"][0]) <= 10]
    candidates = [cluster for cluster in clusters if cluster["passes_candidate_signal_threshold"]]
    payload = {
        "schema": "bub1b_hidden_allele_search_v1", "status": "preliminary",
        "scope": {"chrom": args.chrom, "start": args.start, "end": args.end,
                  "coordinate_system": "1-based inclusive", "bases": args.end - args.start + 1,
                  "gene": {"symbol": "BUB1B", "ensembl": mane["id"], "start": gene_start, "end": gene_end,
                           "strand": mane["strand"], "mane": "ENST00000287598.11/NM_001211.6"},
                  "flanks": {"five_prime_bases": gene_start - args.start,
                             "three_prime_bases": args.end - gene_end}},
        "filters": {"primary_nonduplicate_qcpass": True, "minimum_mapq": 20, "minimum_aligned_query": 30,
                    "minimum_softclip": 15, "minimum_cigar_insertion": 10, "minimum_event_mean_baseq": 20,
                    "cluster_radius_bases": 5, "candidate_minimum_independent_templates": 3,
                    "template_deduplication": "SHA-256 of QNAME; overlapping mates count once"},
        "alignment_counts": {"sam_records": len(rows), "target_overlapping_records": sum(row["on_target"] for row in rows),
                             "eligible_primary_records": sum(eligible(row) for row in rows),
                             "supplementary_target_records": sum(row["on_target"] and bool(row["flag"] & 2048) for row in rows),
                             "sa_tagged_target_records": sum(row["on_target"] and bool(row["sa"]) for row in rows),
                             "unmapped_fetched_records": sum(bool(row["flag"] & 4) for row in rows)},
        "depth": depth_summary(load_depth(args.depth_q0), load_depth(args.depth_q20), load_depth(args.depth_q60), args.start, args.end),
        "known_snv_extraction_controls": target_base_counts(rows),
        "discordant_and_mate_evidence": discordant_summary(rows),
        "candidate_clusters": candidates,
        "all_signal_clusters": clusters,
        "delly_records": load_delly(args.delly),
        "limitations": [
            "This is a bounded locus-wide short-read search, not a calibrated exclusion of all MEIs or structural variants.",
            "The no-alt UCSC-style GRCh38_standard reference lacks vendor decoys; repetitive donor sequence can remain ambiguous.",
            "A candidate threshold is an evidence triage rule, not pathogenicity or genotype probability.",
            "The two known SNVs validate extraction but are not MEI sensitivity controls.",
        ],
        "source_sha256": {key: sha256(path) for key, path in {
            "regional_sam": args.sam, "depth_q0": args.depth_q0, "depth_q20": args.depth_q20,
            "depth_q60": args.depth_q60, "regional_reference": args.reference,
            "mane_annotation": args.mane, "aluya5_consensus": args.alu, "delly_records": args.delly}.items()},
    }
    (args.output / "search.preliminary.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
