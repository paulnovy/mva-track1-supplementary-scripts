#!/usr/bin/env python3
"""Resumable, bounded-space PLINK2 scoring for a local PGS catalog.

This keeps only one model's matched score rows on disk at a time.  PLINK2 is
still the sole scoring implementation; the Cockpit parser is used for the
per-model result contract.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import gc
import json
import os
import re
import shutil
import subprocess
import tempfile
import hashlib
from pathlib import Path

PGS_RE = re.compile(r"(PGS\d{6,})", re.IGNORECASE)


def chrom(v: str) -> str:
    x = v.strip().removeprefix("chr").removeprefix("CHR").upper()
    return "MT" if x == "M" else x


def autosome(v: str) -> bool:
    try:
        return 1 <= int(chrom(v)) <= 22
    except ValueError:
        return False


def load_pvar(path: Path, autosomes_only: bool) -> dict[tuple[str, str], tuple[str, set[str]]]:
    out: dict[tuple[str, str], tuple[str, set[str]]] = {}
    duplicate: set[tuple[str, str]] = set()
    with path.open(encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader((line for line in fh if not line.startswith("##")), delimiter="\t")
        for row in reader:
            c, pos, vid = row.get("#CHROM", row.get("CHROM", "")), row.get("POS", ""), row.get("ID", "")
            if not c or not pos or not vid or (autosomes_only and not autosome(c)):
                continue
            key = (chrom(c), pos)
            alleles = {str(row.get("REF") or "").upper()}
            alleles.update(a.upper() for a in str(row.get("ALT") or "").split(",") if a)
            if key in out:
                duplicate.add(key)
            else:
                out[key] = (vid, alleles)
    for key in duplicate:
        out.pop(key, None)
    return out


def source_header_rows(path: Path):
    fh = gzip.open(path, "rt", encoding="utf-8", errors="ignore", newline="")
    metadata: dict[str, str] = {}
    header = None
    for line in fh:
        if line.startswith("#"):
            raw = line.strip("#\n ")
            if "=" in raw:
                k, v = raw.split("=", 1)
                metadata[k.strip().lower()] = v.strip()
            continue
        header = line.rstrip("\n").split("\t")
        break
    if not header:
        fh.close()
        raise ValueError("score_header_missing")
    return metadata, csv.DictReader(fh, fieldnames=header, delimiter="\t"), fh


def run(cmd: list[str], log: Path) -> None:
    with log.open("ab") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=fh, check=False)
    if proc.returncode:
        raise RuntimeError(f"subprocess_failed_rc={proc.returncode}")


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def count_records(tool: str, path: Path) -> int:
    proc = subprocess.Popen([tool, "view", "-H", str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    assert proc.stdout is not None
    count = sum(1 for _ in proc.stdout)
    if proc.wait():
        raise RuntimeError("bcftools_record_count_failed")
    return count


def count_view(tool: str, path: Path, *extra: str) -> int:
    proc = subprocess.Popen([tool, "view", "-H", *extra, str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    assert proc.stdout is not None
    count = sum(1 for _ in proc.stdout)
    if proc.wait():
        raise RuntimeError("bcftools_filtered_record_count_failed")
    return count


def disk_floor(work: Path, minimum_bytes: int = 60 * 1024**3) -> None:
    if shutil.disk_usage(work).free < minimum_bytes:
        raise RuntimeError("disk_free_below_60GiB")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--input-vcf", type=Path, required=True, help="VCF/VCF.GZ or BCF input")
    ap.add_argument("--catalog-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--plink2", default="plink2")
    ap.add_argument("--helper-dir", type=Path, required=True)
    ap.add_argument("--reference-build-id", default="GRCh38_standard")
    ap.add_argument("--reference-fasta-path", default="")
    ap.add_argument("--contig-style", default="chr")
    ap.add_argument("--manifest-version", default="full-catalog-harmonized-GRCh38")
    ap.add_argument("--autosomes-only", action="store_true", default=True)
    ap.add_argument("--min-free-gib", type=int, default=60,
                    help="free-space floor; production default is 60 GiB")
    args = ap.parse_args()
    args.input_vcf = args.input_vcf.resolve()
    if args.input_vcf.suffix.lower() not in {".vcf", ".gz", ".bcf"}:
        raise SystemExit("input_must_be_vcf_vcfgz_or_bcf")
    w = args.work_dir.resolve()
    w.mkdir(parents=True, exist_ok=True)
    log = w / "bounded-prs.log"
    checkpoint = w / "checkpoint.json"
    result_path = w / f"{args.sample_id}.prs_results.json"
    pgen = w / f"{args.sample_id}.catalog_input"
    pvar = Path(str(pgen) + ".pvar")
    psam = Path(str(pgen) + ".psam")
    pgen_bin = Path(str(pgen) + ".pgen")
    input_sha = sha256_file(args.input_vcf)
    catalog_files = sorted(args.catalog_dir.glob("PGS*.txt.gz"))
    catalog_inventory = [{"name": f.name, "size": f.stat().st_size, "mtime_ns": f.stat().st_mtime_ns}
                         for f in catalog_files]
    identity = {"path": str(args.input_vcf), "size": args.input_vcf.stat().st_size,
                "mtime_ns": args.input_vcf.stat().st_mtime_ns, "sha256": input_sha}
    args_identity = {"sample_id": args.sample_id, "input": identity, "catalog_dir": str(args.catalog_dir.resolve()),
                     "catalog_inventory": catalog_inventory,
                     "reference_build_id": args.reference_build_id, "autosomes_only": True,
                     "plink_threads": 2, "plink_memory_mb": 4096, "min_free_gib": args.min_free_gib}
    if checkpoint.exists():
        old = json.loads(checkpoint.read_text())
        if old.get("run_identity") != args_identity:
            raise SystemExit("checkpoint_input_or_arguments_changed")
    filtered_input = w / ("filtered_input.bcf")
    if not filtered_input.exists():
        # Restrict to PASS/dot, biallelic records, and complete diploid calls;
        # missing calls are never silently treated as homozygous reference.
        run(["bcftools", "view", "-f", "PASS,.", "-m2", "-M2", "-i",
             'GT="0/0" || GT="0/1" || GT="1/0" || GT="1/1" || GT="0|0" || GT="0|1" || GT="1|0" || GT="1|1"',
             "-Ob", "-o", str(filtered_input), str(args.input_vcf)], log)
    input_records = count_records("bcftools", args.input_vcf)
    pass_records = count_view("bcftools", args.input_vcf, "-f", "PASS,.")
    multiallelic_records = count_view("bcftools", args.input_vcf, "-m3")
    missing_gt_records = count_view("bcftools", args.input_vcf, "-i", 'GT="." || GT="./." || GT=".|." || GT~"\\."')
    retained_records = count_records("bcftools", filtered_input)
    if not all(p.exists() and p.stat().st_size > 0 for p in (pgen_bin, pvar, psam)):
        for p in (pgen_bin, pvar, psam):
            p.unlink(missing_ok=True)
        run([args.plink2, "--bcf", str(filtered_input), "--chr", "1-22", "--vcf-half-call", "missing",
             "--max-alleles", "2", "--set-all-var-ids", "@:#", "--new-id-max-allele-len", "100", "missing",
             "--threads", "2", "--memory", "4096", "--make-pgen", "--out", str(pgen)], log)
    if not all(p.exists() and p.stat().st_size > 0 for p in (pgen_bin, pvar, psam)):
        raise RuntimeError("incomplete_pgen_output")
    variants = load_pvar(pvar, args.autosomes_only)
    files: dict[str, Path] = {}
    for f in catalog_files:
        m = PGS_RE.search(f.name)
        if m:
            files.setdefault(m.group(1).upper(), f)
    state = json.loads(checkpoint.read_text()) if checkpoint.exists() else {
        "schema": "bounded_prs_checkpoint_v1", "sample_id": args.sample_id,
        "catalog_count": len(files), "completed": {}, "failed": {}, "autosomes_only": True,
        "run_identity": args_identity, "input_records": input_records, "retained_records": retained_records,
    }
    completed = state.setdefault("completed", {})
    failed = state.setdefault("failed", {})
    parsed = []
    for index, (pgs_id, source) in enumerate(files.items(), 1):
        if pgs_id in completed:
            failed.pop(pgs_id, None)
            parsed.append(completed[pgs_id])
            continue
        score_path: Path | None = None
        fh = None
        try:
            disk_floor(w, args.min_free_gib * 1024**3)
            metadata, rows, fh = source_header_rows(source)
            if metadata.get("hmpos_build") != "GRCh38":
                fh.close()
                failed[pgs_id] = {"reason": "harmonized_build_missing_or_mismatched", "source_file": str(source), "index": index}
                atomic_json(checkpoint, state)
                continue
            weight_type = (metadata.get("weight_type") or "").lower()
            # Log-OR/log-HR are valid additive coefficients, unlike raw OR/HR.
            unsupported_markers = ("haplo", "diplo", "non-add", "interaction")
            raw_ratio = weight_type.strip() in {"or", "hr", "odds ratio", "hazard ratio"}
            if raw_ratio or any(marker in weight_type for marker in unsupported_markers):
                fh.close()
                failed[pgs_id] = {"reason": "unsupported_nonadditive_or_or_only_weight_type",
                                  "source_file": str(source), "index": index}
                atomic_json(checkpoint, state)
                continue
            source_rows = valid_rows = matched = excluded = malformed = nonfinite = duplicate = ambiguous = 0
            seen: set[str] = set()
            seen_weights: dict[str, tuple[str, str]] = {}
            ambiguous_ids: set[str] = set()
            seen_rows: dict[tuple[str, str], tuple[str, str]] = {}
            ambiguous_keys: set[tuple[str, str]] = set()
            with tempfile.NamedTemporaryFile("wt", dir=w, prefix="score-", suffix=".tsv", delete=False) as out:
                score_path = Path(out.name)
                writer = csv.writer(out, delimiter="\t", lineterminator="\n")
                writer.writerow(["ID", "EFFECT_ALLELE", pgs_id])
                for row in rows:
                    source_rows += 1
                    if any(str(row.get(k, '')).strip().lower() in {'true', '1', 'yes'}
                           for k in ('is_haplotype', 'is_diplotype', 'is_interaction')):
                        raise ValueError('unsupported_haplotype_diplotype_or_interaction_model')
                    c, pos = row.get("hm_chr") or "", row.get("hm_pos") or ""
                    effect, weight = (row.get("effect_allele") or "").upper(), row.get("effect_weight") or ""
                    if not c or not pos or not effect or not weight:
                        malformed += 1
                        continue
                    try:
                        value = float(weight)
                    except ValueError:
                        malformed += 1
                        continue
                    if not (value == value and abs(value) != float("inf")):
                        nonfinite += 1
                        continue
                    valid_rows += 1
                    if not autosome(c):
                        excluded += 1
                        continue
                    key = (chrom(c), pos)
                    present = variants.get(key)
                    if not present:
                        continue
                    # Only matched positions need ambiguity tracking. Keeping
                    # all 10M catalog rows would defeat bounded-memory scoring.
                    previous_row = seen_rows.get(key)
                    if previous_row is not None:
                        duplicate += 1
                        if previous_row != (effect, weight):
                            ambiguous_keys.add(key)
                    else:
                        seen_rows[key] = (effect, weight)
                    vid, alleles = present
                    if effect not in alleles:
                        continue
                    previous = seen_weights.get(vid)
                    if previous is not None:
                        if previous != (effect, weight):
                            ambiguous_ids.add(vid)
                        continue
                    seen_weights[vid] = (effect, weight)
                    seen.add(vid)
                    matched += 1
                    writer.writerow([vid, effect, weight])
            fh.close()
            fh = None
            for key in ambiguous_keys:
                present = variants.get(key)
                if present:
                    ambiguous_ids.add(present[0])
            if ambiguous_ids:
                ambiguous = len(ambiguous_ids)
                filtered_score = score_path.with_suffix(".filtered.tsv")
                with score_path.open(encoding="utf-8") as src, filtered_score.open("w", encoding="utf-8") as dst:
                    for line in src:
                        if line.startswith("ID\t"):
                            dst.write(line)
                        elif line.split("\t", 1)[0] not in ambiguous_ids:
                            dst.write(line)
                score_path.unlink()
                score_path = filtered_score
                matched -= ambiguous
            meta = w / f"{pgs_id}.metadata.json"
            metadata_doc = {
                "schema": "plink_pgs_catalog_metadata_v1", "target_build": "GRCh38",
                "catalog_count": 1, "runnable_count": int(matched > 0), "skipped_count": 0,
                "items": [{"pgs_id": pgs_id, "trait": metadata.get("trait_reported") or pgs_id,
                    "variant_count_total": valid_rows, "variant_count_matched": matched,
                    "source_row_count": source_rows, "malformed_row_count": malformed,
                    "nonfinite_weight_count": nonfinite, "autosome_excluded_variant_count": excluded,
                    "ambiguous_duplicate_count": ambiguous, "duplicate_rows_skipped": duplicate,
                    "source_file": str(source),
                    "score_build_id": "GRCh38", "harmonized_date": metadata.get("hmpos_date")}],
            }
            meta.write_text(json.dumps(metadata_doc, indent=2) + "\n", encoding="utf-8")
            del seen, seen_weights, seen_rows, ambiguous_keys, ambiguous_ids
            gc.collect()
            if matched:
                prefix = w / f"{pgs_id}.score"
                disk_floor(w, args.min_free_gib * 1024**3)
                run([args.plink2, "--pfile", str(pgen), "--score", str(score_path), "1", "2", "3",
                     "header-read", "no-mean-imputation", "ignore-dup-ids", "cols=+scoresums,-scoreavgs",
                     "--threads", "2", "--memory", "4096",
                     "--out", str(prefix)], log)
                sscore = Path(str(prefix) + ".sscore")
                parsed_path = w / f"{pgs_id}.result.json"
                run(["python3", str(args.helper_dir / "parse_plink_pgs_scores.py"), "--metadata", str(meta),
                     "--sscore", str(sscore), "--out", str(parsed_path), "--reference-build-id", args.reference_build_id,
                     "--reference-fasta-path", args.reference_fasta_path, "--contig-style", args.contig_style,
                     "--sample-build", "GRCh38", "--manifest-version", args.manifest_version], log)
                item = json.loads(parsed_path.read_text())["items"][0]
            else:
                item = {"pgs_id": pgs_id, "trait": metadata.get("trait_reported") or pgs_id,
                    "score_value": None, "percentile": None, "risk_band": None,
                    "overlap_pct": 0.0, "variant_count_total": valid_rows, "variant_count_matched": 0,
                    "match_rate": 0.0, "quality_label": "not_interpretable", "non_diagnostic": True,
                    "interpretation_status": "raw_score_only", "calibration_status": "missing",
                    "source_row_count": source_rows, "malformed_row_count": malformed,
                    "nonfinite_weight_count": nonfinite, "autosome_excluded_variant_count": excluded,
                    "ambiguous_duplicate_count": ambiguous}
            item["autosome_excluded_variant_count"] = excluded
            item["source_row_count"] = source_rows
            item["malformed_row_count"] = malformed
            item["nonfinite_weight_count"] = nonfinite
            item["ambiguous_duplicate_count"] = ambiguous
            item["source_file"] = str(source)
            item["score_file_retained"] = bool(matched)
            completed[pgs_id] = item
            failed.pop(pgs_id, None)
            atomic_json(checkpoint, state)
            parsed.append(item)
            score_path.unlink(missing_ok=True)
            if matched:
                # Keep sscore/QC and provenance; remove only temporary expanded score text.
                pass
        except RuntimeError:
            # A scoring/conversion subprocess failure is a genuine run error:
            # checkpoint what finished and stop for operator diagnosis.
            failed[pgs_id] = {"reason": "subprocess_failed", "source_file": str(source), "index": index}
            atomic_json(checkpoint, state)
            raise
        except ValueError as exc:
            # Bad catalog input is accounted for and does not prevent the
            # remaining independent models from completing.
            failed[pgs_id] = {"reason": str(exc), "source_file": str(source), "index": index}
            atomic_json(checkpoint, state)
            continue
        finally:
            if fh is not None:
                fh.close()
            if score_path is not None:
                score_path.unlink(missing_ok=True)
    atomic_json(checkpoint, state)
    atomic_json(result_path, {"status": "complete", "schema": "prs_result_batch_v1", "source_tool": "plink2_score_list_harmonized_bounded",
        "catalog_count": len(files), "score_column_count": sum(1 for x in parsed if x.get("variant_count_matched", 0)),
        "skipped_catalog_count": len(failed), "count": len(parsed), "all_catalogs_accounted": len(completed) + len(failed) == len(files),
        "input_records": input_records, "pass_or_dot_records": pass_records,
        "multiallelic_records": multiallelic_records, "missing_gt_records": missing_gt_records,
        "retained_records": retained_records,
        "autosomes_only": True, "raw_score_only": True, "items": parsed})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
