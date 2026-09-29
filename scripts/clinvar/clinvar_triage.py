#!/usr/bin/env python3
"""Extract and compare exact ClinVar P/LP calls from two annotated VCFs.

Detailed sample-level values are written only to a caller-selected private JSON
file.  Standard output contains a deliberately limited aggregate with gene
symbols, never loci, alleles, genotypes, depths, or likelihoods.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import pysam


PLP = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}


def terms(value: Any) -> set[str]:
    if value is None:
        return set()
    values = value if isinstance(value, tuple) else (value,)
    return {str(item).lower().replace(" ", "_") for item in values if item not in (None, ".", "")}


def scalar(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    return value


def allele_key(record: pysam.VariantRecord) -> tuple[str, int, str, str, int]:
    if len(record.alts or ()) != 1:
        raise ValueError("Expected a biallelic normalized record")
    return record.contig, record.pos, record.ref, record.alts[0], int(record.info["ALLELEID"])


def genotype_state(gt: Any) -> str:
    if gt is None or any(item is None for item in gt):
        return "missing_or_partial"
    nonref = sum(item > 0 for item in gt)
    if not nonref:
        return "reference"
    if len(gt) == 1:
        return "hemizygous_alt"
    if nonref == len(gt):
        return "homozygous_alt"
    return "heterozygous_alt"


def called_pl_index(gt: Any) -> int | None:
    if gt is None or len(gt) != 2 or any(item is None for item in gt):
        return None
    a, b = sorted(gt)
    return b * (b + 1) // 2 + a


def call_details(record: pysam.VariantRecord, sample_name: str) -> dict[str, Any]:
    call = record.samples[sample_name]
    gt = call.get("GT")
    ad = scalar(call.get("AD"))
    explicit_dp = call.get("DP")
    ad_sum = sum(x for x in ad if isinstance(x, int) and x >= 0) if ad else None
    depth = explicit_dp if explicit_dp is not None else ad_sum
    alt_depth = ad[1] if ad and len(ad) > 1 and isinstance(ad[1], int) and ad[1] >= 0 else None
    allele_balance = alt_depth / ad_sum if alt_depth is not None and ad_sum else None
    pl = scalar(call.get("PL"))
    pl_index = called_pl_index(gt)
    called_pl_is_min = None
    pl_margin = None
    if pl and pl_index is not None and pl_index < len(pl) and all(isinstance(x, int) for x in pl):
        called_pl_is_min = pl[pl_index] == min(pl)
        ordered = sorted(pl)
        if len(ordered) > 1:
            pl_margin = ordered[1] - ordered[0]
    filters = list(record.filter.keys())
    filter_label = "PASS" if filters == ["PASS"] else ("dot" if not filters else ";".join(filters))
    state = genotype_state(gt)
    called_alt_present = any(item is not None and item > 0 for item in gt or ())
    checks = {
        "called_alt_present": called_alt_present,
        "complete_nonreference": state not in {"missing_or_partial", "reference"},
        "filter_pass": filter_label == "PASS",
        "depth_at_least_10": depth >= 10 if depth is not None else None,
        "gq_at_least_20": call.get("GQ") >= 20 if call.get("GQ") is not None else None,
        "called_pl_is_min": called_pl_is_min,
        "pl_margin_at_least_20": pl_margin >= 20 if pl_margin is not None else None,
    }
    if state == "heterozygous_alt":
        checks["heterozygous_allele_balance_0.20_to_0.80"] = (
            0.20 <= allele_balance <= 0.80 if allele_balance is not None else None
        )
    return {
        "genotype_state": state,
        "gt": list(gt) if gt is not None else None,
        "requires_original_genotype_review": state == "missing_or_partial" and called_alt_present,
        "genotype_caveat": ("A called ALT is retained for review. Recover original GT/allele mapping/AD/PL; split AD balance and PL are not complete-genotype evidence."
                            if state == "missing_or_partial" and called_alt_present else None),
        "phased": bool(call.phased),
        "filter": filter_label,
        "qual": record.qual,
        "depth": depth,
        "depth_source": "FORMAT/DP" if explicit_dp is not None else ("sum(FORMAT/AD)" if ad_sum is not None else None),
        "ref_depth": ad[0] if ad and isinstance(ad[0], int) and ad[0] >= 0 else None,
        "alt_depth": alt_depth,
        "allele_balance": allele_balance,
        "gq": call.get("GQ"),
        "pl_present": pl is not None,
        "called_pl_is_min": called_pl_is_min,
        "pl_margin": pl_margin,
        "checks": checks,
    }


def collect(path: Path) -> tuple[str, dict[tuple[str, int, str, str, int], dict[str, Any]]]:
    selected: dict[tuple[str, int, str, str, int], dict[str, Any]] = {}
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        if len(samples) != 1:
            raise ValueError(f"Expected exactly one sample column in {path}")
        sample_name = samples[0]
        for record in vcf:
            if "ALLELEID" not in record.info or not (terms(record.info.get("CLNSIG")) & PLP):
                continue
            details = call_details(record, sample_name)
            if not details["checks"]["called_alt_present"]:
                continue
            key = allele_key(record)
            if key in selected:
                raise ValueError("Duplicate exact candidate allele")
            geneinfo = str(record.info.get("GENEINFO") or "")
            selected[key] = {
                "locus": {"contig": key[0], "position": key[1], "ref": key[2], "alt": key[3]},
                "allele_id": key[4],
                "genes": sorted({item.split(":", 1)[0] for item in geneinfo.split("|") if item}),
                "clinvar": {
                    "significance": scalar(record.info.get("CLNSIG")),
                    "review_status": scalar(record.info.get("CLNREVSTAT")),
                    "conditions": scalar(record.info.get("CLNDN")),
                    "molecular_consequence": scalar(record.info.get("MC")),
                },
                "call": details,
            }
    return sample_name, selected


def raw_clinvar_details(path: Path, candidates: dict[tuple[str, int, str, str, int], dict[str, Any]]) -> None:
    with pysam.VariantFile(str(path)) as vcf:
        raw_contigs = set(vcf.header.contigs)
        for key, candidate in candidates.items():
            contig = key[0]
            if contig not in raw_contigs and contig.startswith("chr") and contig[3:] in raw_contigs:
                contig = contig[3:]
            matches = []
            for record in vcf.fetch(contig, key[1] - 1, key[1]):
                if record.info.get("ALLELEID") == key[4]:
                    matches.append(record)
            if not matches:
                candidate["clinvar_evidence"] = {"raw_record_found_by_allele_id": False}
                continue
            record = matches[0]
            scv_values = scalar(record.info.get("CLNSIGSCV")) or []
            if not isinstance(scv_values, list):
                scv_values = [scv_values]
            # ClinVar encodes multiple SCV accessions inside a pipe-delimited
            # String value even though the INFO header also declares Number=.
            scvs = [part for value in scv_values for part in str(value).split("|") if part]
            candidate["clinvar_evidence"] = {
                "raw_record_found_by_allele_id": True,
                "submission_accession_count": len(scvs),
                "conflicting_classifications": scalar(record.info.get("CLNSIGCONF")),
                "condition_database_ids": scalar(record.info.get("CLNDISDB")),
                "origin_codes": scalar(record.info.get("ORIGIN")),
                "variation_id": record.id,
                "hgvs": scalar(record.info.get("CLNHGVS")),
                "dbsnp_ids": scalar(record.info.get("RS")),
                "raw_significance": scalar(record.info.get("CLNSIG")),
                "raw_review_status": scalar(record.info.get("CLNREVSTAT")),
            }


def reciprocal_route_check(path: Path, candidate: dict[str, Any]) -> dict[str, Any]:
    locus = candidate["locus"]
    with pysam.VariantFile(str(path)) as vcf:
        samples = list(vcf.header.samples)
        records = list(vcf.fetch(locus["contig"], locus["position"] - 1, locus["position"]))
        exact = [r for r in records if r.pos == locus["position"] and r.ref == locus["ref"] and tuple(r.alts or ()) == (locus["alt"],)]
        result = {
            "same_position_record_count": len(records),
            "exact_allele_record_present": bool(exact),
        }
        if exact and samples:
            result["exact_record_genotype_state"] = genotype_state(exact[0].samples[samples[0]].get("GT"))
            result["exact_record_has_same_allele_id"] = exact[0].info.get("ALLELEID") == candidate["allele_id"]
        return result


def targeted_snv_pileup(path: Path, candidate: dict[str, Any]) -> dict[str, Any]:
    locus = candidate["locus"]
    if len(locus["ref"]) != 1 or len(locus["alt"]) != 1:
        return {"performed": False, "reason": "SNV-only check"}
    counts = {f"{base}_{strand}": 0 for base in ("ref", "alt", "other") for strand in ("forward", "reverse")}
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for column in bam.pileup(
            locus["contig"], locus["position"] - 1, locus["position"], truncate=True,
            stepper="all", min_base_quality=0,
        ):
            if column.reference_pos != locus["position"] - 1:
                continue
            for pileup_read in column.pileups:
                read = pileup_read.alignment
                if (
                    pileup_read.is_del or pileup_read.is_refskip or read.is_unmapped or read.is_secondary
                    or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality < 30
                ):
                    continue
                query_position = pileup_read.query_position
                if query_position is None or read.query_qualities[query_position] < 25:
                    continue
                base = read.query_sequence[query_position].upper()
                kind = "ref" if base == locus["ref"] else ("alt" if base == locus["alt"] else "other")
                strand = "reverse" if read.is_reverse else "forward"
                counts[f"{kind}_{strand}"] += 1
    ref_count = counts["ref_forward"] + counts["ref_reverse"]
    alt_count = counts["alt_forward"] + counts["alt_reverse"]
    total = sum(counts.values())
    return {
        "performed": True,
        "thresholds": "primary non-duplicate, non-QC-fail, MAPQ>=30, baseQ>=25; SNV observations",
        "counts": counts,
        "eligible_total": total,
        "alt_fraction_of_ref_alt": alt_count / (ref_count + alt_count) if ref_count + alt_count else None,
        "alt_both_strands": counts["alt_forward"] > 0 and counts["alt_reverse"] > 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--core", type=Path)
    parser.add_argument("--clinvar", required=True, type=Path)
    parser.add_argument("--existing", type=Path, help="Re-enrich an existing private extraction without rescanning VCFs")
    parser.add_argument("--bam", type=Path, help="Optional BAM for focused checks requested with --pileup-gene")
    parser.add_argument("--pileup-gene", action="append", default=[])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if args.existing:
        with args.existing.open(encoding="utf-8") as handle:
            result = json.load(handle)
        combined = result["candidates"]
        combined_by_key = {}
        for item in combined:
            locus = item["locus"]
            key = (locus["contig"], locus["position"], locus["ref"], locus["alt"], item["allele_id"])
            combined_by_key[key] = item
        raw_clinvar_details(args.clinvar, combined_by_key)
        result["inputs"]["clinvar"] = str(args.clinvar)
        source_path = Path(result["inputs"]["source"])
        core_path = Path(result["inputs"]["core"])
        for item in combined:
            if "source" not in item["routes"]:
                item.setdefault("reciprocal_route_checks", {})["source"] = reciprocal_route_check(source_path, item)
            if "core" not in item["routes"]:
                item.setdefault("reciprocal_route_checks", {})["core"] = reciprocal_route_check(core_path, item)
            if args.bam and set(item["genes"]) & set(args.pileup_gene):
                item["targeted_bam_check"] = targeted_snv_pileup(args.bam, item)
    else:
        if args.source is None or args.core is None:
            parser.error("--source and --core are required unless --existing is used")
        source_sample, source = collect(args.source)
        core_sample, core = collect(args.core)
        all_keys = sorted(set(source) | set(core))
        combined = []
        for key in all_keys:
            template = source.get(key) or core[key]
            item = {name: value for name, value in template.items() if name != "call"}
            item["routes"] = {}
            if key in source:
                item["routes"]["source"] = source[key]["call"]
            if key in core:
                item["routes"]["core"] = core[key]["call"]
            item["route_concordance"] = {
                "exact_allele_in_both": key in source and key in core,
                "genotype_state_equal": (
                    source[key]["call"]["genotype_state"] == core[key]["call"]["genotype_state"]
                    if key in source and key in core else None
                ),
            }
            combined.append(item)
        combined_by_key = {key: item for key, item in zip(all_keys, combined)}
        raw_clinvar_details(args.clinvar, combined_by_key)

        result = {
            "format": "clinvar-local-triage-v1",
            "privacy": "Contains individual-level variant and call details; keep protected locally.",
            "inputs": {
                "source": str(args.source),
                "core": str(args.core),
                "clinvar": str(args.clinvar),
                "source_sample_column_present": bool(source_sample),
                "core_sample_column_present": bool(core_sample),
            },
            "selection": "Exact prepared-ClinVar allele; explicitly called ALT (including partial GT); aggregate CLNSIG P/LP; partial genotypes require original allele/GT/AD/PL review and are not imputed or accepted for inheritance.",
            "counts": {
                "source": len(source),
                "core": len(core),
                "union": len(all_keys),
                "exact_in_both": sum(key in source and key in core for key in all_keys),
                "genotype_state_concordant": sum(
                    key in source and key in core and source[key]["call"]["genotype_state"] == core[key]["call"]["genotype_state"]
                    for key in all_keys
                ),
            },
            "candidates": combined,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    old_umask = os.umask(0o077)
    try:
        with args.output.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
    finally:
        os.umask(old_umask)

    safe = {
        "candidate_count": len(combined),
        "exact_in_both": result["counts"]["exact_in_both"],
        "genotype_state_concordant": result["counts"]["genotype_state_concordant"],
        "genes": sorted({gene for item in combined for gene in item["genes"]}),
    }
    print(json.dumps(safe, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
