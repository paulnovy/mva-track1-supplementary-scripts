#!/usr/bin/env python3
"""Summarize protected targeted SAM/depth evidence for BUB1B MEI and CLCNKB SV calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
BUB_BOUNDARY = 40_196_543  # between chr15:40196542 and 40196543 (GRCh38)
CLCN_LEFT = 16_050_024     # reconstructed split boundary, between 16050023/24
CLCN_RIGHT = 16_059_763    # reconstructed split boundary, between 16059762/63


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def digest_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def longest_exact_query(query: str, reference: str, floor: int = 8) -> int:
    query = query.upper()
    reference = reference.upper()
    for width in range(min(len(query), len(reference), 40), floor - 1, -1):
        if any(query[index:index + width] in reference for index in range(len(query) - width + 1)):
            return width
    return 0


def parse_cigar(cigar: str, start: int) -> tuple[int, int, int, list[tuple[int, int]]]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    ref = start
    insertions: list[tuple[int, int]] = []
    for length, op in operations:
        if op in "MDN=X":
            ref += length
        elif op == "I":
            insertions.append((ref, length))
    leading = operations[0][0] if operations and operations[0][1] == "S" else 0
    trailing = operations[-1][0] if operations and operations[-1][1] == "S" else 0
    return ref - 1, leading, trailing, insertions


def sam_rows(path: Path, alu: str) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & 4:
                continue
            pos = int(fields[3]); mapq = int(fields[4]); cigar = fields[5]
            end, leading, trailing, insertions = parse_cigar(cigar, pos)
            tags = {part[:2]: part[5:] for part in fields[11:] if len(part) >= 6 and part[2:5] in {":Z:", ":i:"}}
            def tag_int(name: str):
                value = tags.get(name, "")
                return int(value) if value.lstrip("-").isdigit() else None
            base = {
                "qname_hash": digest_text(fields[0]), "flag": flag, "rname": fields[2], "pos": pos,
                "end": end, "mapq": mapq, "cigar": cigar, "mate_rname": fields[6],
                "mate_pos": int(fields[7]) if fields[7].lstrip("-").isdigit() else 0,
                "tlen": int(fields[8]) if fields[8].lstrip("-").isdigit() else 0,
                "reverse": bool(flag & 16), "duplicate": bool(flag & 1024), "qcfail": bool(flag & 512),
                "secondary": bool(flag & 256), "supplementary": bool(flag & 2048),
                "sa": tags.get("SA", ""), "nm": int(tags["NM"]) if tags.get("NM", "").isdigit() else None,
                "s1": tag_int("s1"), "s2": tag_int("s2"), "tp": tags.get("tp", ""),
                "insertions": insertions, "clips": [],
            }
            sequence = fields[9]
            for side, length, boundary, clipped in (
                ("leading", leading, pos, sequence[:leading]),
                ("trailing", trailing, end + 1, sequence[-trailing:] if trailing else ""),
            ):
                if length:
                    base["clips"].append({
                        "side": side, "length": length, "boundary": boundary,
                        "poly_a_or_t_10": "A" * 10 in clipped.upper() or "T" * 10 in clipped.upper(),
                        "alu_exact_match": longest_exact_query(clipped, alu),
                        "alu_antisense_exact_match": longest_exact_query(clipped, revcomp(alu)),
                    })
            rows.append(base)
    return rows


def eligible(row: dict, mapq: int = 0) -> bool:
    return not (row["duplicate"] or row["qcfail"] or row["secondary"] or row["supplementary"]) and row["mapq"] >= mapq


def clips_near(rows: list[dict], boundary: int, radius: int = 5, minimum: int = 10) -> list[tuple[dict, dict]]:
    return [
        (row, clip) for row in rows if eligible(row)
        for clip in row["clips"] if clip["length"] >= minimum and abs(clip["boundary"] - boundary) <= radius
    ]


def clcn_directional_split_support(rows: list[dict], radius: int = 5, minimum: int = 10) -> list[tuple[dict, dict]]:
    return [
        (row, clip) for row in rows if eligible(row)
        for clip in row["clips"] if clip["length"] >= minimum and (
            (clip["side"] == "trailing" and abs(row["end"] - (CLCN_LEFT - 1)) <= radius)
            or (clip["side"] == "leading" and abs(row["pos"] - CLCN_RIGHT) <= radius)
        )
    ]


def summarize_support(items: list[tuple[dict, dict]]) -> dict:
    fragments = {row["qname_hash"] for row, _ in items}
    return {
        "alignments": len(items), "unique_fragments": len(fragments),
        "mapq_ge_20_fragments": len({row["qname_hash"] for row, _ in items if row["mapq"] >= 20}),
        "mapq_60_fragments": len({row["qname_hash"] for row, _ in items if row["mapq"] >= 60}),
        "forward_fragments": len({row["qname_hash"] for row, _ in items if not row["reverse"]}),
        "reverse_fragments": len({row["qname_hash"] for row, _ in items if row["reverse"]}),
        "sa_tag_fragments": len({row["qname_hash"] for row, _ in items if row["sa"]}),
        "poly_a_or_t_10_fragments": len({row["qname_hash"] for row, clip in items if clip["poly_a_or_t_10"]}),
        "alu_match_ge_15_fragments": len({row["qname_hash"] for row, clip in items if max(clip["alu_exact_match"], clip["alu_antisense_exact_match"]) >= 15}),
        "boundary_counts": dict(sorted(Counter(clip["boundary"] for _, clip in items).items())),
        "side_counts": dict(Counter(clip["side"] for _, clip in items)),
    }


def depth_values(path: Path) -> dict[int, int]:
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            _, position, depth = line.rstrip("\n").split("\t")
            values[int(position)] = int(depth)
    return values


def interval_depth(values: dict[int, int], start: int, end: int) -> dict:
    observations = [values.get(position, 0) for position in range(start, end + 1)]
    return {
        "start": start, "end": end, "bases": len(observations),
        "mean": round(statistics.fmean(observations), 4), "median": statistics.median(observations),
        "minimum": min(observations), "maximum": max(observations),
        "bases_ge_10_fraction": round(sum(value >= 10 for value in observations) / len(observations), 6),
        "bases_ge_20_fraction": round(sum(value >= 20 for value in observations) / len(observations), 6),
    }


def window_metrics(rows: list[dict], depths: dict[int, int], start: int, end: int) -> dict:
    selected = [row for row in rows if eligible(row, 20) and row["pos"] <= end and row["end"] >= start]
    return {
        **interval_depth(depths, start, end),
        "primary_nonduplicate_mapq20_fragments": len({row["qname_hash"] for row in selected}),
        "generic_clip_ge_10_fragments": len({
            row["qname_hash"] for row in selected for clip in row["clips"]
            if clip["length"] >= 10 and start <= clip["boundary"] <= end
        }),
        "sa_tag_fragments": len({row["qname_hash"] for row in selected if row["sa"]}),
        "abnormal_mate_fragments": len({
            row["qname_hash"] for row in selected
            if row["mate_rname"] not in {"=", "chr15"} or row["mate_pos"] == 0 or abs(row["tlen"]) > 1_000
        }),
    }


def clip_clusters(rows: list[dict], minimum: int = 10) -> list[dict]:
    grouped: dict[int, list[tuple[dict, dict]]] = defaultdict(list)
    for row in rows:
        if not eligible(row):
            continue
        for clip in row["clips"]:
            if clip["length"] >= minimum:
                grouped[clip["boundary"]].append((row, clip))
    return [
        {"boundary": boundary, **summarize_support(items)}
        for boundary, items in sorted(grouped.items(), key=lambda pair: (-len({x[0]["qname_hash"] for x in pair[1]}), pair[0]))
    ]


def sa_links(items: list[tuple[dict, dict]], target_boundary: int, radius: int = 15) -> dict:
    linked = set(); partners = Counter(); paralog = set()
    for row, _ in items:
        for entry in row["sa"].split(";"):
            if not entry:
                continue
            fields = entry.split(",")
            if len(fields) < 6:
                continue
            chrom, raw_pos, strand, cigar, raw_mapq, raw_nm = fields[:6]
            try:
                pos = int(raw_pos); mapq = int(raw_mapq)
            except ValueError:
                continue
            partners[(chrom, pos, strand, cigar, mapq)] += 1
            if chrom == "chr1" and abs(pos - target_boundary) <= radius:
                linked.add(row["qname_hash"])
            if chrom == "chr1" and 16_015_000 <= pos <= 16_040_000:
                paralog.add(row["qname_hash"])
    return {
        "unique_fragments_linked_to_reciprocal_breakpoint": len(linked),
        "unique_fragments_with_sa_in_clcnka_paralog_window": len(paralog),
        "top_sa_partners": [
            {"chrom": key[0], "position": key[1], "strand": key[2], "cigar": key[3], "mapq": key[4], "alignments": count}
            for key, count in partners.most_common(12)
        ],
    }


def paired_breakpoint_support(rows: list[dict], left: int, right: int, radius: int = 750) -> tuple[dict, set[str]]:
    by_fragment: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if eligible(row):
            by_fragment[row["qname_hash"]].append(row)
    selected: dict[str, list[dict]] = {}
    for name, alignments in by_fragment.items():
        for row in alignments:
            mate_reverse = bool(row["flag"] & 32)
            inward = (
                abs(row["pos"] - left) <= radius and abs(row["mate_pos"] - right) <= radius
                and not row["reverse"] and mate_reverse
            ) or (
                abs(row["pos"] - right) <= radius and abs(row["mate_pos"] - left) <= radius
                and row["reverse"] and not mate_reverse
            )
            if row["mate_rname"] in {"=", "chr1"} and row["mate_pos"] and inward:
                selected[name] = alignments
                break
    both_mapq20 = set()
    ambiguous = set()
    for name, alignments in selected.items():
        first = max((row["mapq"] for row in alignments if row["flag"] & 64), default=-1)
        second = max((row["mapq"] for row in alignments if row["flag"] & 128), default=-1)
        if first >= 20 and second >= 20:
            both_mapq20.add(name)
        if any(row["s1"] is not None and row["s2"] is not None and row["s1"] - row["s2"] <= 10 for row in alignments):
            ambiguous.add(name)
    return {
        "unique_fragments": len(selected),
        "both_alignments_mapq_ge_20_fragments": len(both_mapq20),
        "minimap2_s1_minus_s2_le_10_fragments": len(ambiguous),
        "radius_bases": radius,
        "orientation": "inward-facing FR",
    }, set(selected)


def summarize_realignment(path: Path) -> dict:
    hits: dict[str, list[tuple[str, int]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            target = "CLCNKA" if fields[5].startswith("chr1:16021000-") else "CLCNKB"
            hits[fields[0]].append((target, int(fields[11])))
    best = Counter(); any_target = Counter()
    for alignments in hits.values():
        for target in {target for target, _ in alignments}:
            any_target[target] += 1
        top = max(mapq for _, mapq in alignments)
        targets = {target for target, mapq in alignments if mapq == top}
        best[next(iter(targets)) if len(targets) == 1 else "tie"] += 1
    return {
        "read_ends": len(hits),
        "any_clcnka": any_target["CLCNKA"], "any_clcnkb": any_target["CLCNKB"],
        "best_clcnka": best["CLCNKA"], "best_clcnkb": best["CLCNKB"], "best_tie": best["tie"],
    }


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bub-sam", type=Path, required=True)
    parser.add_argument("--bub-depth", type=Path, required=True)
    parser.add_argument("--clcn-sam", type=Path, required=True)
    parser.add_argument("--clcn-depth", type=Path, required=True, help="minimum mapping quality 20; no base-quality floor")
    parser.add_argument("--clcn-depth-mapq20", type=Path, required=True)
    parser.add_argument("--clcn-depth-mapq60", type=Path, required=True)
    parser.add_argument("--dfam", type=Path, required=True)
    parser.add_argument("--delly-records", type=Path, required=True)
    parser.add_argument("--clcn-realignment", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True, mode=0o700)

    alu = json.loads(args.dfam.read_text(encoding="utf-8"))["consensus_sequence"]
    bub = sam_rows(args.bub_sam, alu)
    clcn = sam_rows(args.clcn_sam, alu)
    bub_support = clips_near(bub, BUB_BOUNDARY)
    bub_primary_fragments = {row["qname_hash"] for row in bub if eligible(row)}
    bub_primary_mapq20_fragments = {row["qname_hash"] for row in bub if eligible(row, 20)}
    bub_discordant_fragments = {
        row["qname_hash"] for row in bub if eligible(row)
        and (row["mate_rname"] not in {"=", "chr15"} or row["mate_pos"] == 0 or abs(row["tlen"]) > 1_000)
    }
    bub_spanning = {
        row["qname_hash"] for row in bub if eligible(row, 20)
        and row["pos"] <= BUB_BOUNDARY - 1 and row["end"] >= BUB_BOUNDARY
        and not any(abs(site - BUB_BOUNDARY) <= 5 for site, _ in row["insertions"])
    }
    bub_insertions = {
        row["qname_hash"] for row in bub if eligible(row)
        for site, length in row["insertions"] if abs(site - BUB_BOUNDARY) <= 5 and length >= 10
    }
    left_support = clips_near(clcn, CLCN_LEFT, radius=10)
    right_support = clips_near(clcn, CLCN_RIGHT, radius=10)
    exact_split = clcn_directional_split_support(clcn)
    exact_split_fragments = {row["qname_hash"] for row, _ in exact_split}
    exact_split_mapq20_fragments = {row["qname_hash"] for row, _ in exact_split if row["mapq"] >= 20}
    left_fragments = {row["qname_hash"] for row, _ in left_support}
    right_fragments = {row["qname_hash"] for row, _ in right_support}

    bub_depth = depth_values(args.bub_depth)
    clcn_depth = depth_values(args.clcn_depth)
    clcn_depth_mapq20 = depth_values(args.clcn_depth_mapq20)
    clcn_depth_mapq60 = depth_values(args.clcn_depth_mapq60)
    intervals = {
        "bub1b_breakpoint_window": interval_depth(bub_depth, BUB_BOUNDARY - 500, BUB_BOUNDARY + 500),
        "clcn_window_a": (16_035_712, 16_058_679),
        "clcn_window_b": (16_037_378, 16_060_547),
        "clcn_target_interval": (16_050_022, 16_059_762),
        "clcnka_gene": (16_022_036, 16_034_050),
        "clcnkb_gene": (16_043_782, 16_057_326),
        "control_left": (15_990_000, 16_005_000),
        "control_right": (16_065_000, 16_080_000),
        "control_distant": (17_000_000, 17_015_000),
    }
    intervals.pop("bub1b_breakpoint_window")
    bub_windows = {
        "upstream_control": window_metrics(bub, bub_depth, BUB_BOUNDARY - 1500, BUB_BOUNDARY - 501),
        "target": window_metrics(bub, bub_depth, BUB_BOUNDARY - 500, BUB_BOUNDARY + 499),
        "downstream_control": window_metrics(bub, bub_depth, BUB_BOUNDARY + 500, BUB_BOUNDARY + 1499),
    }
    depth = {"bub1b_windows_mapq20_no_baseq_floor": bub_windows}
    for label, (start, end) in intervals.items():
        depth[label] = {
            "mapq20_no_baseq_floor": interval_depth(clcn_depth, start, end),
            "baseq20_mapq20": interval_depth(clcn_depth_mapq20, start, end),
            "mapq60_no_baseq_floor": interval_depth(clcn_depth_mapq60, start, end),
        }
    controls = [depth[key]["mapq20_no_baseq_floor"]["mean"] for key in ("control_left", "control_right", "control_distant")]
    control_mean = statistics.fmean(controls)
    depth_ratio = depth["clcn_target_interval"]["mapq20_no_baseq_floor"]["mean"] / control_mean if control_mean else None
    pair_metrics = {}; pair_sets = {}
    for ident, left, right in (
        ("clcn_window_a", 16_035_711, 16_058_679),
        ("clcn_window_b", 16_037_377, 16_060_547),
        ("clcn_target", 16_050_021, 16_059_762),
    ):
        pair_metrics[ident], pair_sets[ident] = paired_breakpoint_support(clcn, left, right)

    payload = {
        "schema": "mva1_targeted_sv_mei_evidence_v1", "assembly": "GRCh38",
        "filters": {
            "clip_minimum_bases": 10, "breakpoint_radius_bases": {"bub1b": 5, "clcnkb": 10},
            "primary_nonduplicate_non_qcfail": True,
            "depth_profiles": ["MAPQ>=20", "baseQ>=20 and MAPQ>=20", "MAPQ>=60"],
            "fragment_identity": "SHA-256 of read name; original read names not emitted",
        },
        "bub1b_published_alu": {
            "literature": {"doi": "10.1038/hgv.2017.21", "pmcid": "PMC5462940"},
            "event": "AluYa5 antisense insertion with 16-nt TSD in intron 8 polypyrimidine tract",
            "tsd": "GTAATTTTGCTTCTTT", "mapped_tsd_grch38": "chr15:40196527-40196542",
            "mapped_acceptor_ag_grch38": "chr15:40196543-40196544",
            "interrogated_boundary": "chr15:40196542|40196543",
            "mapping_basis": "TSD exact sequence plus MANE Select ENST00000287598/NM_001211.6 exon 9 start",
            "target_and_control_windows": bub_windows,
            "clip_support": summarize_support(bub_support),
            "window_primary_nonduplicate_fragments": len(bub_primary_fragments),
            "window_primary_nonduplicate_mapq20_fragments": len(bub_primary_mapq20_fragments),
            "window_discordant_or_unmapped_mate_fragments": len(bub_discordant_fragments),
            "reference_spanning_mapq20_primary_nonduplicate_fragments": len(bub_spanning),
            "cigar_insertions_ge_10bp_within_5bp_fragments": len(bub_insertions),
            "top_clip_clusters_in_2kb_window": clip_clusters(bub)[:20],
            "status": "metrics_computed; interpretation_required",
            "limit": (
                "Absence of a short-read clip/insertion cluster is not a validated genome-wide MEI sensitivity estimate "
                "and does not exclude a different BUB1B insertion or a difficult-to-map/truncated event."
            ),
        },
        "clcnkb_deletion": {
            "target_labels": list(pair_metrics),
            "interrogated_boundaries": ["chr1:16050023|16050024", "chr1:16059762|16059763"],
            "exact_directional_split_support": {
                "alignments": len(exact_split), "unique_fragments": len(exact_split_fragments),
                "mapq_ge_20_alignments": sum(row["mapq"] >= 20 for row, _ in exact_split),
                "mapq_ge_20_unique_fragments": len(exact_split_mapq20_fragments),
                "predicate": "trailing clip >=10 bp ending at 16050023 +/-5, or leading clip >=10 bp starting at 16059763 +/-5",
            },
            "supporting_read_realignment": summarize_realignment(args.clcn_realignment),
            "left_breakpoint_support": summarize_support(left_support),
            "right_breakpoint_support": summarize_support(right_support),
            "fragments_seen_at_both_local_clip_clusters": len(left_fragments & right_fragments),
            "paired_breakpoint_proxy_support": pair_metrics,
            "paired_proxy_fragment_set_intersections": {
                "window_a_window_b": len(pair_sets["clcn_window_a"] & pair_sets["clcn_window_b"]),
                "window_a_target": len(pair_sets["clcn_window_a"] & pair_sets["clcn_target"]),
                "window_b_target": len(pair_sets["clcn_window_b"] & pair_sets["clcn_target"]),
            },
            "target_paired_proxy_fragments_contained_in_split_set": len(pair_sets["clcn_target"] & exact_split_fragments),
            "left_sa_specificity": sa_links(left_support, CLCN_RIGHT),
            "right_sa_specificity": sa_links(right_support, CLCN_LEFT),
            "target_to_mean_control_depth_ratio_mapq20": round(depth_ratio, 5) if depth_ratio is not None else None,
            "status": "metrics_computed; dosage_and_paralog_interpretation_required",
            "limit": (
                "Overlapping call records must not be assumed to be independent alleles. Short-read alignments in the "
                "homologous CLCNKA/CLCNKB locus are not orthogonal validation or a standalone dosage-loss assessment. "
                "DELLY input is fingerprinted only; its genotypes are not parsed by this script."
            ),
        },
        "depth": depth,
        "source_sha256": {
            "bub_sam": sha256(args.bub_sam), "bub_depth": sha256(args.bub_depth),
            "clcn_sam": sha256(args.clcn_sam), "clcn_depth": sha256(args.clcn_depth),
            "clcn_depth_mapq20": sha256(args.clcn_depth_mapq20),
            "clcn_depth_mapq60": sha256(args.clcn_depth_mapq60),
            "clcn_realignment": sha256(args.clcn_realignment),
            "dfam": sha256(args.dfam), "delly_records": sha256(args.delly_records),
        },
    }
    atomic_json(args.output / "targeted-evidence.json", payload)
    with (args.output / "evidence-status.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ("candidate", "evidence_question", "status", "key_observation", "remaining_limit")
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([
            {
                "candidate": "BUB1B published Kato AluYa5 breakpoint",
                "evidence_question": "same exact breakpoint in this WGS",
                "status": payload["bub1b_published_alu"]["status"],
                "key_observation": (
                    f"{len(bub_support)} qualifying clips; {len(bub_insertions)} >=10-bp CIGAR insertions; "
                    f"{len(bub_spanning)} reference-spanning MAPQ>=20 primary nonduplicate fragments"
                ),
                "remaining_limit": payload["bub1b_published_alu"]["limit"],
            },
            {
                "candidate": "CLCNKB target interval",
                "evidence_question": "local breakpoint support",
                "status": payload["clcnkb_deletion"]["status"],
                "key_observation": (
                    f"left/right qualifying clip fragments {len(left_fragments)}/{len(right_fragments)}; "
                    f"target-to-control MAPQ>=20 depth ratio {payload['clcnkb_deletion']['target_to_mean_control_depth_ratio_mapq20']}"
                ),
                "remaining_limit": payload["clcnkb_deletion"]["limit"],
            },
            {
                "candidate": "CLCNKB flanking target windows",
                "evidence_question": "independent second allele",
                "status": "independent_alleles_not_assessed",
                "key_observation": "Pair-support counts and fragment intersections are in paired-support-summary.tsv and targeted-evidence.json.",
                "remaining_limit": "Fragment overlap alone does not establish an independent structural allele or exact haplotype.",
            },
            {
                "candidate": "Other BUB1B mobile-element insertion",
                "evidence_question": "genome-wide or whole-gene calibrated MEI exclusion",
                "status": "uncallable by this breakpoint-specific assay",
                "key_observation": "Only the sequence-mapped Kato intron-8 AluYa5 breakpoint was interrogated with targeted controls",
                "remaining_limit": "A different insertion, truncated element, or poorly anchored event requires a validated broader MEI workflow or orthogonal assay.",
            },
        ])
    os.chmod(args.output / "evidence-status.tsv", 0o600)
    with (args.output / "depth-summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("interval", "profile", "start", "end", "mean_depth", "median_depth", "bases_ge_20_fraction"))
        for label, profiles in depth.items():
            if label == "bub1b_windows_mapq20_no_baseq_floor":
                for window, metrics in profiles.items():
                    writer.writerow((f"BUB1B_{window}", "MAPQ>=20", metrics["start"], metrics["end"], metrics["mean"], metrics["median"], metrics["bases_ge_20_fraction"]))
            else:
                for profile, metrics in profiles.items():
                    writer.writerow((label, profile, metrics["start"], metrics["end"], metrics["mean"], metrics["median"], metrics["bases_ge_20_fraction"]))
    os.chmod(args.output / "depth-summary.tsv", 0o600)
    with (args.output / "paired-support-summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("event", "proxy_fragments", "both_mapq20_fragments", "s1_minus_s2_le_10_fragments"))
        for ident, metrics in pair_metrics.items():
            writer.writerow((ident, metrics["unique_fragments"], metrics["both_alignments_mapq_ge_20_fragments"], metrics["minimap2_s1_minus_s2_le_10_fragments"]))
    os.chmod(args.output / "paired-support-summary.tsv", 0o600)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
