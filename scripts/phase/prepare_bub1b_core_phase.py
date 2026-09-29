#!/usr/bin/env python3
"""Prepare one supplemental core-callset BUB1B read-backed phase run.

Coordinates are read from a local protected evidence JSON at execution time and
are intentionally not embedded in this repository or printed by this tool.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


class ValidationError(RuntimeError):
    pass


def parse_target(item: object) -> tuple[str, int, str, str]:
    try:
        contig, position = str(item["coordinate"]).rsplit(":", 1)
        position = int(position.replace(",", ""))
        ref = str(item["ref"])
        alts = item["alt"]
        if not isinstance(alts, list) or len(alts) != 1:
            raise ValueError
        alt = str(alts[0])
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        raise ValidationError("Protected evidence does not contain a usable target coordinate.") from exc
    if not contig or position < 1 or not ref or not alt or any(base not in "ACGT" for base in ref + alt):
        raise ValidationError("Protected evidence contains an invalid target coordinate.")
    return contig, position, ref, alt


def main(argv: list[str] | None = None) -> int:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path, help="Protected local-evidence.json")
    parser.add_argument("--scope", required=True, type=Path, help="Prior protected/public BUB1B scope.json defining gene plus flanks")
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--core-vcf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="New protected output directory")
    parser.add_argument("--whatshap", default="whatshap")
    parser.add_argument("--min-qual", type=float, default=30.0, help="Inclusive core-call QUAL gate (default: 30)")
    args = parser.parse_args(argv)
    try:
        if args.min_qual < 0:
            raise ValidationError("--min-qual must be nonnegative.")
        payload = json.loads(args.evidence.read_text())
        candidates = payload.get("candidates")
        if not isinstance(candidates, list):
            raise ValidationError("Protected evidence lacks its candidate list.")
        # Use only targets independently observed in the core callset.  This is
        # deliberately not the historical source-only candidate VCF.
        targets = [parse_target(item) for item in candidates if item.get("core_bcftools")]
        if len(targets) < 2:
            raise ValidationError("Fewer than two protected targets are represented in the core callset; supplemental phase is unresolved.")
        if len(set(targets)) != len(targets):
            raise ValidationError("Protected evidence has duplicate target alleles; supplemental phase is unresolved.")
        contigs = {contig for contig, _, _, _ in targets}
        if len(contigs) != 1:
            raise ValidationError("Protected targets are not on one contig; no regional phase run is appropriate.")
        contig = next(iter(contigs))
        positions = [position for _, position, _, _ in targets]
        scope = json.loads(args.scope.read_text())
        region = str(scope.get("region", ""))
        try:
            scope_contig, interval = region.rsplit(":", 1)
            scope_start, scope_end = (int(value.replace(",", "")) for value in interval.split("-", 1))
        except (ValueError, AttributeError) as exc:
            raise ValidationError("Scope does not define a usable gene-plus-flanks region.") from exc
        if scope_contig != contig or scope_start < 1 or scope_end < scope_start or any(not scope_start <= position <= scope_end for position in positions):
            raise ValidationError("Scope does not contain every protected target on its contig.")
        tool = Path(__file__).with_name("prepare_regional_phasing.py")
        command = [sys.executable, str(tool), "run", "--bam", str(args.bam), "--vcf", str(args.core_vcf), "--fasta", str(args.fasta), "--region", region, "--out", str(args.out), "--whatshap", args.whatshap, "--min-qual", str(args.min_qual), "--min-dp", "10", "--dp-from-ad", "--min-gq", "20", "--min-pl-delta", "20", "--min-ad", "3"]
        for contig, position, ref, alt in targets:
            command.extend(["--target", f"{contig}:{position}:{ref}:{alt}"])
        result = subprocess.run(command, check=False)
        if result.returncode:
            return result.returncode
        summary = json.loads((args.out / "phase-block-summary.json").read_text())
        missing_targets = summary["targets"]["found_in_candidate_vcf"] != summary["targets"]["requested"]
        no_eligible = summary["candidate_records"] == 0
        outcome = "selection_no_eligible_candidates" if no_eligible else ("selection_target_allele_missing" if missing_targets else "phase_completed")
        result_summary = {"status": "unresolved" if no_eligible or missing_targets else "complete", "outcome": outcome, "format": "bub1b-core-supplemental-phase-v1", "candidate_records": summary["candidate_records"], "phase_block_count": summary["phase_block_count"], "targets_requested": summary["targets"]["requested"], "targets_found": summary["targets"]["found_in_candidate_vcf"], "targets_same_phase_set": summary["targets"]["targets_same_phase_set"], "target_pair_cis_or_trans": summary["targets"]["target_pair_cis_or_trans"], "scope_source": str(args.scope), "thresholds": {"qual": args.min_qual, "dp": 10, "depth_method": "sum_FORMAT/AD when FORMAT/DP absent", "gq": 20, "pl_best_vs_second": 20, "ad_when_available": 3}}
        (args.out / "summary.json").write_text(json.dumps(result_summary, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result_summary, sort_keys=True))
        return 3 if no_eligible or missing_targets else 0
    except (OSError, json.JSONDecodeError, ValidationError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
