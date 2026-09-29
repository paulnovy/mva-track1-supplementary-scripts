#!/usr/bin/env python3
"""Prepare a blocked, offline AlphaGenome API stage manifest.

This utility deliberately stops at preparation.  It does not import an SDK,
read credentials or contact the API. It copies only the sensitive target
variant definitions, not sample identifiers, genotypes, or phenotypes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Any


class ValidationError(RuntimeError):
    pass


TARGET = "gdmscience.googleapis.com"
WINDOW = 1_048_576


def _variants(payload: dict[str, Any]) -> list[dict[str, Any]]:
    if payload.get("candidate_count") != 2:
        raise ValidationError("Atlas summary must declare candidate_count=2.")
    if payload.get("scope") != "local exact-allele AlphaGenome Atlas annotation of two BUB1B SNVs":
        raise ValidationError("Atlas summary has an unexpected scope.")
    records = payload.get("results")
    if not isinstance(records, list):
        raise ValidationError("Atlas summary must contain a results list.")
    if len(records) != 2:
        raise ValidationError("Atlas summary must contain exactly two variants.")
    result: list[dict[str, Any]] = []
    for item in records:
        if not isinstance(item, dict):
            raise ValidationError("Each Atlas variant must be an object.")
        if set(("id", "chrom", "pos", "ref", "alt")) - item.keys():
            raise ValidationError("Each Atlas result must contain id, chrom, pos, ref, and alt.")
        chrom = item["chrom"]
        position = item["pos"]
        ref = item["ref"]
        alt = item["alt"]
        if chrom == "15":
            chrom = "chr15"
        elif chrom != "chr15":
            raise ValidationError("Atlas results must target chromosome 15.")
        if not isinstance(position, int) or isinstance(position, bool):
            raise ValidationError("Each variant needs a 1-based chromosome and integer position.")
        if position < 1 or not isinstance(ref, str) or not isinstance(alt, str):
            raise ValidationError("Each variant needs valid REF and ALT alleles.")
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT" or ref == alt:
            raise ValidationError("The two targets must be simple, non-reference A/C/G/T SNVs.")
        result.append({"chrom": chrom, "pos": position, "ref": ref, "alt": alt})
    if len({(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in result}) != 2:
        raise ValidationError("Atlas summary contains duplicate target variants.")
    return sorted(result, key=lambda v: (v["chrom"], v["pos"], v["ref"], v["alt"]))


def _interval(chrom: str, pos: int) -> tuple[dict[str, Any], dict[str, Any]]:
    # SDK genome.Interval is zero-based half-open.  The official SDK center()
    # rounds this odd-width SNV interval up to pos, then resize(2**20) expands
    # equally: https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/data/genome.py
    start, end = pos - 1, pos
    center = pos
    query_start, query_end = center - WINDOW // 2, center + WINDOW // 2
    if query_start < 0:
        raise ValidationError("A centered 1 MiB interval would extend before chromosome start.")
    return ({"chrom": chrom, "start": start, "end": end, "coordinate_system": "zero_based_half_open"},
            {"chrom": chrom, "start": query_start, "end": query_end, "coordinate_system": "zero_based_half_open"})


def prepare(source: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise ValidationError(f"Output already exists; refusing overwrite: {output}")
    if not source.is_file():
        raise ValidationError(f"Atlas summary is not a regular file: {source}")
    try:
        raw = source.read_bytes()
        payload = json.loads(raw)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"Cannot read valid JSON Atlas summary: {source}") from exc
    if not isinstance(payload, dict):
        raise ValidationError("Atlas summary root must be an object.")
    if payload.get("status") != "complete":
        raise ValidationError("Atlas summary must have status=complete.")
    build = payload.get("reference_build", payload.get("genome_build"))
    if build != "GRCh38":
        raise ValidationError("Atlas summary must declare reference_build=GRCh38.")
    variants = _variants(payload)
    source_sha256 = hashlib.sha256(raw).hexdigest()
    manifest = {
        "schema": "alphagenome_api_stage_manifest_v1",
        "status": "BLOCKED",
        "phase": {"status": "not_started", "reason": "supported transport and renewed patient egress approval required"},
        "source": {"sha256": source_sha256, "format": "protected_local_atlas_summary_json"},
        "target": TARGET,
        "reference_build": "GRCh38",
        "scope": "proposed AlphaGenome API differential scoring of two BUB1B SNVs",
        "variants": [],
        "requested_scorers": ["RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"],
        "requested_score_scope": ["differential_splice", "differential_expression"],
        "resource_ceiling": {"cpu": 5, "memory_gib": 20},
        "execution": {
            "ordered_steps": [
                "resolve supported API transport for the exact target",
                "obtain renewed explicit patient egress approval",
                "submit the two coordinate-and-allele requests",
                "checkpoint differential splice and expression scores",
            ],
            "checkpoints": ["manifest_validated", "transport_authorized", "requests_submitted", "scores_checkpointed"],
        },
        "privacy": {"sensitive_genomic_coordinates_and_alleles": True, "ids_phenotypes_genotypes_excluded": True, "network_requests_made": False},
    }
    for v in variants:
        ref_interval, query_interval = _interval(v["chrom"], v["pos"])
        manifest["variants"].append({**v, "ref_interval": ref_interval, "query_interval_1mb": query_interval})
    output.mkdir(mode=0o700)
    manifest_path = output / "manifest.json"
    readme_path = output / "README.txt"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    readme_path.write_text(
        "AlphaGenome API stage is prepared but BLOCKED.\n"
        "This manifest is offline-only and contains two sensitive genomic GRCh38 SNV definitions.\n"
        "A real runner may be implemented only after a safe supported transport to "
        f"{TARGET} and renewed explicit patient egress approval are available.\n"
        "Do not bypass transport protection, add credentials to this directory, or treat this as queued/running.\n"
    )
    manifest_path.chmod(0o600)
    readme_path.chmod(0o600)
    return manifest


def main(argv: list[str] | None = None) -> int:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas-summary", required=True, type=Path, help="Protected existing local Atlas summary JSON")
    parser.add_argument("--out", required=True, type=Path, help="New protected output directory")
    args = parser.parse_args(argv)
    try:
        manifest = prepare(args.atlas_summary, args.out)
    except (OSError, ValidationError) as exc:
        print(f"ERROR: {exc}")
        return 2
    print(json.dumps({"status": manifest["status"], "phase": manifest["phase"]["status"], "output": str(args.out)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
