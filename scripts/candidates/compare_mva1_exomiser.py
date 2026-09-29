#!/usr/bin/env python3
"""Compare baseline/full/axis Exomiser gene-model and allele rankings."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


PROFILES = ("baseline", "full", "axis")
FOCUS_GENES = {"BUB1B", "BUB1", "ATP6V1E1", "CLCNKA", "CLCNKB", "LZTR1", "FANCD2", "PIK3C2A"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader((line[1:] if line.startswith("#") else line for line in handle), delimiter="\t"))


def number(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(args.output, 0o700)

    genes: dict[str, list[dict[str, str]]] = {}
    variants: dict[str, list[dict[str, str]]] = {}
    source_files: dict[str, dict[str, str]] = {}
    for profile in PROFILES:
        native = args.runs / profile / "work" / "rare-disease" / "native"
        gene_path = native / f"MVA1-{profile}-exomiser.genes.tsv"
        variant_path = native / f"MVA1-{profile}-exomiser.variants.tsv"
        genes[profile] = table(gene_path)
        variants[profile] = table(variant_path)
        source_files[profile] = {
            "genes": str(gene_path), "genes_sha256": sha256(gene_path),
            "variants": str(variant_path), "variants_sha256": sha256(variant_path),
        }

    by_profile = {profile: {row["ID"]: row for row in rows} for profile, rows in genes.items()}
    ids = sorted(set().union(*(set(rows) for rows in by_profile.values())), key=lambda key: (
        min(int(by_profile[p].get(key, {}).get("RANK", "999999")) for p in PROFILES), key
    ))
    rank_fields = [
        "model_id", "gene_symbol", "moi",
        *[f"{profile}_{item}" for profile in PROFILES for item in ("rank", "combined_score", "phenotype_score", "variant_score")],
        "full_minus_baseline_rank", "axis_minus_full_rank",
    ]
    with (args.output / "rank-comparison.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rank_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for ident in ids:
            any_row = next(by_profile[p][ident] for p in PROFILES if ident in by_profile[p])
            out = {"model_id": ident, "gene_symbol": any_row["GENE_SYMBOL"], "moi": any_row["MOI"]}
            for profile in PROFILES:
                row = by_profile[profile].get(ident, {})
                out.update({
                    f"{profile}_rank": row.get("RANK", ""),
                    f"{profile}_combined_score": row.get("EXOMISER_GENE_COMBINED_SCORE", ""),
                    f"{profile}_phenotype_score": row.get("EXOMISER_GENE_PHENO_SCORE", ""),
                    f"{profile}_variant_score": row.get("EXOMISER_GENE_VARIANT_SCORE", ""),
                })
            for name, left, right in (("full_minus_baseline_rank", "full", "baseline"), ("axis_minus_full_rank", "axis", "full")):
                try:
                    out[name] = int(by_profile[left][ident]["RANK"]) - int(by_profile[right][ident]["RANK"])
                except KeyError:
                    out[name] = ""
            writer.writerow(out)
    os.chmod(args.output / "rank-comparison.tsv", 0o600)

    excluded_rows = [row for row in genes["full"] if row["GENE_SYMBOL"] != "BUB1B"]
    with (args.output / "full-bub1b-excluded-view.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ("display_rank_after_filter", "original_rank", "model_id", "gene_symbol", "moi", "combined_score", "phenotype_score", "variant_score")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for display_rank, row in enumerate(excluded_rows, 1):
            writer.writerow({
                "display_rank_after_filter": display_rank, "original_rank": row["RANK"], "model_id": row["ID"],
                "gene_symbol": row["GENE_SYMBOL"], "moi": row["MOI"],
                "combined_score": row["EXOMISER_GENE_COMBINED_SCORE"],
                "phenotype_score": row["EXOMISER_GENE_PHENO_SCORE"],
                "variant_score": row["EXOMISER_GENE_VARIANT_SCORE"],
            })
    os.chmod(args.output / "full-bub1b-excluded-view.tsv", 0o600)

    focus_models = []
    for ident in ids:
        any_row = next(by_profile[p][ident] for p in PROFILES if ident in by_profile[p])
        ranks = {profile: int(by_profile[profile][ident]["RANK"]) for profile in PROFILES if ident in by_profile[profile]}
        if any_row["GENE_SYMBOL"] in FOCUS_GENES or min(ranks.values()) <= 20:
            focus_models.append({
                "model_id": ident, "gene_symbol": any_row["GENE_SYMBOL"], "moi": any_row["MOI"],
                "profiles": {
                    profile: {
                        "rank": int(by_profile[profile][ident]["RANK"]),
                        "combined_score": number(by_profile[profile][ident]["EXOMISER_GENE_COMBINED_SCORE"]),
                        "phenotype_score": number(by_profile[profile][ident]["EXOMISER_GENE_PHENO_SCORE"]),
                        "variant_score": number(by_profile[profile][ident]["EXOMISER_GENE_VARIANT_SCORE"]),
                        "human_phenotype_evidence": by_profile[profile][ident].get("HUMAN_PHENO_EVIDENCE", ""),
                    } for profile in PROFILES if ident in by_profile[profile]
                },
            })

    allele_rows = []
    seen = set()
    for profile in PROFILES:
        for row in variants[profile]:
            if row["GENE_SYMBOL"] not in FOCUS_GENES and int(row["RANK"]) > 20:
                continue
            key = (profile, row["ID"], row["CONTIG"], row["START"], row["REF"], row["ALT"])
            if key in seen:
                continue
            seen.add(key)
            allele_rows.append({
                "profile": profile, "model_rank": int(row["RANK"]), "model_id": row["ID"],
                "gene_symbol": row["GENE_SYMBOL"], "moi": row["MOI"], "contig": row["CONTIG"],
                "start": int(row["START"]), "ref": row["REF"], "alt": row["ALT"],
                "genotype": row["GENOTYPE"], "functional_class": row["FUNCTIONAL_CLASS"], "hgvs": row["HGVS"],
                "contributing_variant": row["CONTRIBUTING_VARIANT"] == "1",
                "variant_score": number(row["EXOMISER_VARIANT_SCORE"]),
                "engine_acmg_classification": row["EXOMISER_ACMG_CLASSIFICATION"],
                "clinvar_interpretation": row["CLINVAR_PRIMARY_INTERPRETATION"],
                "max_frequency_source": row["MAX_FREQ_SOURCE"], "max_frequency": number(row["MAX_FREQ"]),
                "frequency_missing": not bool(row["MAX_FREQ_SOURCE"] or row["MAX_FREQ"]),
            })

    atomic_json(args.output / "comparison.json", {
        "schema": "mva1_exomiser_three_way_comparison_v1",
        "profiles": list(PROFILES), "source_files": source_files,
        "interpretation_contract": {
            "scores_are_probabilities": False, "scores_are_diagnoses": False,
            "axis_is_actual_phenotype_ablation_rerun": True,
            "bub1b_excluded_view_is_rerun": False,
            "bub1b_excluded_view_note": "Display-only filtering; all scores and original ranks are unchanged.",
            "candidate_rows_do_not_establish_phase": True,
            "two_rows_or_two_variants_do_not_establish_biallelic_in_trans": True,
        },
        "row_counts": {profile: {"gene_models": len(genes[profile]), "variants": len(variants[profile])} for profile in PROFILES},
        "variant_annotation_coverage": {
            profile: {
                "rows": len(variants[profile]),
                "functional_class_present": sum(bool(row["FUNCTIONAL_CLASS"]) for row in variants[profile]),
                "hgvs_present": sum(bool(row["HGVS"]) for row in variants[profile]),
                "population_frequency_present": sum(bool(row["MAX_FREQ_SOURCE"] and row["MAX_FREQ"]) for row in variants[profile]),
                "clinvar_interpretation_present": sum(row["CLINVAR_PRIMARY_INTERPRETATION"] not in {"", "NOT_PROVIDED"} for row in variants[profile]),
                "source": "Exomiser GRCh38 genome data 2512 fields, not annotations embedded in the input VCF",
            } for profile in PROFILES
        },
        "focus_models": focus_models, "focus_alleles": allele_rows,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
