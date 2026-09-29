#!/usr/bin/env python3
"""Run strict, exact allele-level ClinVar annotation and aggregate-only QC.

The input VCF is never overwritten.  With --normalize-input this creates a
new biallelic, REF-checked copy before annotation.  Without it, the caller is
responsible for supplying an already normalized VCF; it is still REF-checked.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pysam

CLINVAR_COLUMNS = "INFO/ALLELEID,INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/MC"
CLINVAR_IDS = {item.split("/")[1] for item in CLINVAR_COLUMNS.split(",")}


class ValidationError(RuntimeError):
    pass


def run(command: list[str]) -> None:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        raise ValidationError(f"bcftools failed with exit status {exc.returncode}; no result was published.") from exc


def vcf_index(path: Path) -> Path:
    for suffix in (".csi", ".tbi"):
        candidate = Path(str(path) + suffix)
        if candidate.is_file():
            return candidate
    raise ValidationError("VCF index is required.")


def require_vcf(path: Path, label: str) -> None:
    if not path.is_file() or not str(path).endswith((".vcf.gz", ".vcf.bgz")):
        raise ValidationError(f"{label} must be an indexed bgzip VCF.")
    vcf_index(path)
    try:
        with pysam.VariantFile(str(path)) as vcf:
            next(iter(vcf.header.contigs), None)
    except Exception as exc:
        raise ValidationError(f"{label} cannot be opened by htslib.") from exc


def header_info_ids(path: Path) -> set[str]:
    with pysam.VariantFile(str(path)) as vcf:
        return set(vcf.header.info.keys())


def lineage_counts(path: Path) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            counts["records"] += 1
            counts["alt_alleles"] += len(record.alts or ())
    return dict(counts)


def partition_nonacgtn_ref(source: Path, temp: Path) -> tuple[Path, Path, dict[str, int]]:
    """Split only raw records whose REF cannot be reference-normalized.

    Quarantined records are copied verbatim into a BCF under the original
    header.  They are not repaired, annotated, or silently dropped.
    """
    eligible_plain = temp / "eligible-original.vcf"
    quarantined = temp / "quarantine-original.bcf"
    counts: Counter[str] = Counter()
    with pysam.VariantFile(str(source)) as input_vcf, \
            pysam.VariantFile(str(eligible_plain), "w", header=input_vcf.header) as eligible, \
            pysam.VariantFile(str(quarantined), "wb", header=input_vcf.header) as quarantine:
        for record in input_vcf:
            counts["raw_input_records"] += 1
            counts["raw_input_alt_alleles"] += len(record.alts or ())
            if not record.ref or any(base not in "ACGTN" for base in record.ref.upper()):
                counts["quarantined_nonacgtn_ref_records"] += 1
                counts["quarantined_nonacgtn_ref_alt_alleles"] += len(record.alts or ())
                quarantine.write(record)
            else:
                counts["eligible_records"] += 1
                counts["eligible_alt_alleles"] += len(record.alts or ())
                eligible.write(record)
    if counts["raw_input_records"] != counts["eligible_records"] + counts["quarantined_nonacgtn_ref_records"]:
        raise ValidationError("Raw input/quarantine record conservation failed.")
    if counts["raw_input_alt_alleles"] != counts["eligible_alt_alleles"] + counts["quarantined_nonacgtn_ref_alt_alleles"]:
        raise ValidationError("Raw input/quarantine ALT lineage conservation failed.")
    eligible_gz = Path(pysam.tabix_index(str(eligible_plain), preset="vcf", force=True))
    return eligible_gz, quarantined, dict(counts)


def validate_all_contigs_against_reference(vcf_path: Path, fasta: Path, label: str) -> None:
    """Validate populated contigs; empty quarantined declarations remain allowed."""
    indexed = subprocess.run(["bcftools", "index", "--stats", str(vcf_path)], capture_output=True, text=True, check=True)
    used = {line.split("\t")[0] for line in indexed.stdout.splitlines() if len(line.split("\t")) >= 3 and int(line.split("\t")[2]) > 0}
    with pysam.FastaFile(str(fasta)) as reference, pysam.VariantFile(str(vcf_path)) as vcf:
        lengths = dict(zip(reference.references, reference.lengths))
        for name, header_contig in vcf.header.contigs.items():
            if name not in used:
                continue
            if name not in lengths:
                raise ValidationError(f"{label} declares contig absent from the exact FASTA: {name!r}.")
            declared = header_contig.length
            if declared is not None and declared != lengths[name]:
                raise ValidationError(f"{label} contig length disagrees with the exact FASTA: {name!r}.")


def locus_records(iterator):
    """Yield sorted VCF records in bounded same-position groups."""
    pending = next(iterator, None)
    while pending is not None:
        locus = (pending.contig, pending.pos)
        group = [pending]
        pending = next(iterator, None)
        while pending is not None and (pending.contig, pending.pos) == locus:
            group.append(pending)
            pending = next(iterator, None)
        yield locus, group


def genotype_multiset(records, samples: list[str]) -> Counter:
    """Exact variant identity plus GT/phasing, preserving duplicate multiplicity."""
    result: Counter = Counter()
    for record in records:
        calls = tuple((sample, record.samples[sample].get("GT"), record.samples[sample].phased) for sample in samples)
        result[(record.ref, record.alts, calls)] += 1
    return result


def compare_gt_semantics(before: Path, after: Path) -> int:
    """Compare exact keyed record/GT multisets, allowing only same-locus reorder."""
    records = 0
    with pysam.VariantFile(str(before)) as left, pysam.VariantFile(str(after)) as right:
        left_samples, right_samples = list(left.header.samples), list(right.header.samples)
        if left_samples != right_samples:
            raise ValidationError("Annotation changed VCF sample columns.")
        left_groups, right_groups = locus_records(iter(left)), locus_records(iter(right))
        while True:
            left_group, right_group = next(left_groups, None), next(right_groups, None)
            if left_group is None and right_group is None:
                return records
            if left_group is None or right_group is None or left_group[0] != right_group[0]:
                raise ValidationError("Annotation changed record locations or record count.")
            if genotype_multiset(left_group[1], left_samples) != genotype_multiset(right_group[1], left_samples):
                raise ValidationError("Annotation changed variant identity or GT/phasing semantics.")
            records += len(left_group[1])


def terms(value) -> set[str]:
    if value is None:
        return set()
    values = value if isinstance(value, tuple) else (value,)
    return {str(item).lower().replace(" ", "_") for item in values if item not in (None, ".", "")}


def aggregate(path: Path) -> dict:
    counts: Counter[str] = Counter()
    review_status: Counter[str] = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            counts["input_retained_records"] += 1
            filter_value = str(record).split("\t")[6]
            counts["filter_dot"] += filter_value == "."
            counts["filter_pass"] += filter_value == "PASS"
            counts["filter_other"] += filter_value not in {".", "PASS"}
            if "ALLELEID" not in record.info:
                counts["unmatched_records"] += 1
                continue
            counts["matched_records"] += 1
            call_states = []
            alt_present = []
            for sample in vcf.header.samples:
                gt = record.samples[sample].get("GT")
                alt_present.append(any(allele is not None and allele > 0 for allele in gt or ()))
                if gt is None or any(allele is None for allele in gt):
                    call_states.append("partial_or_missing")
                elif any(allele > 0 for allele in gt):
                    call_states.append("complete_nonreference")
                else:
                    call_states.append("complete_reference")
            # A sites-only/zero-sample VCF is never a patient non-reference call.
            genotype_state = "complete_nonreference" if call_states and all(item == "complete_nonreference" for item in call_states) else ("partial_or_missing" if "partial_or_missing" in call_states else "complete_reference")
            counts[f"matched_genotype_{genotype_state}"] += 1
            if genotype_state == "complete_nonreference":
                filter_label = "dot" if filter_value == "." else ("pass" if filter_value == "PASS" else "other")
                counts[f"matched_complete_nonreference_filter_{filter_label}"] += 1
            significance = terms(record.info.get("CLNSIG"))
            reviews = terms(record.info.get("CLNREVSTAT"))
            if reviews:
                counts["matched_with_review_status"] += 1
                for value in reviews:
                    review_status[value] += 1
            else:
                counts["matched_without_review_status"] += 1
            conflict = any("conflict" in value for value in significance)
            plp = bool(significance & {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"})
            vus = bool(significance & {"uncertain_significance", "uncertain_significance/likely_benign", "uncertain_significance/likely_pathogenic"})
            # Keep historical complete-GT counters; a called ALT in a partial
            # GT also belongs in an explicit review stratum, without imputation.
            called_alt = bool(alt_present) and all(alt_present)
            if called_alt:
                counts["matched_called_alt"] += 1
                prefix = "matched_called_alt_"
                counts[prefix + "conflicting_classification"] += conflict
                counts[prefix + "pathogenic_or_likely_pathogenic"] += plp
                counts[prefix + "vus_or_uncertain"] += vus
                if genotype_state == "partial_or_missing":
                    counts["matched_partial_called_alt"] += 1
                    counts["matched_partial_called_alt_conflicting_classification"] += conflict
                    counts["matched_partial_called_alt_pathogenic_or_likely_pathogenic"] += plp
                    counts["matched_partial_called_alt_vus_or_uncertain"] += vus
            elif genotype_state == "partial_or_missing":
                counts["matched_partial_without_called_alt"] += 1
            if genotype_state == "complete_nonreference":
                counts["matched_complete_nonreference_conflicting_classification"] += conflict
                counts["matched_complete_nonreference_pathogenic_or_likely_pathogenic"] += plp
                counts["matched_complete_nonreference_nonconflicting_pathogenic_or_likely_pathogenic"] += plp and not conflict
                counts["matched_complete_nonreference_vus_or_uncertain"] += vus
                counts["matched_complete_nonreference_other_or_missing_clnsig"] += not (plp or vus or conflict)
    return {"counts": dict(sorted(counts.items())), "review_status_counts": dict(sorted(review_status.items())),
            "genotype_strata": "Complete-GT counters are preserved. Called-ALT counters include partial calls separately, requiring an explicit ALT in every sample; they are review counts, not validated diagnoses. Recover original allele mapping and GT/AD/PL before inheritance interpretation."}


def main(argv: list[str] | None = None) -> int:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--clinvar", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="New private result directory")
    parser.add_argument("--normalize-input", action="store_true", help="Strictly normalize an independent raw/core VCF first")
    args = parser.parse_args(argv)
    try:
        source, clinvar, fasta, out = (item.resolve() for item in (args.input, args.clinvar, args.fasta, args.out))
        require_vcf(source, "Input")
        require_vcf(clinvar, "ClinVar")
        if not fasta.is_file() or not Path(str(fasta) + ".fai").is_file():
            raise ValidationError("Exact FASTA and .fai are required.")
        # Prepared ClinVar deliberately retains unsupported-contig declarations
        # while quarantining their records, so validate populated records only.
        validate_all_contigs_against_reference(source, fasta, "Input")
        if out.exists() or not out.parent.is_dir():
            raise ValidationError("Output directory must not exist and its parent must exist.")
        collisions = CLINVAR_IDS & header_info_ids(source)
        if collisions:
            raise ValidationError("Input already declares ClinVar INFO fields; refusing to overwrite annotations.")
        source_lineage = lineage_counts(source)
        source_stat = source.stat()
        with tempfile.TemporaryDirectory(prefix=f".{out.name}.", dir=out.parent) as temporary:
            temp = Path(temporary)
            normalized = temp / "normalized.vcf.gz"
            quarantine = None
            partition_counts = None
            if args.normalize_input:
                eligible, quarantine, partition_counts = partition_nonacgtn_ref(source, temp)
                # -m -any is the supported split; --check-ref e aborts, never swaps.
                run(["bcftools", "norm", "--check-ref", "e", "-f", str(fasta), "-m", "-any", "--multi-overlaps", ".", "--old-rec-tag", "NORMALIZATION_ORIGINAL", "-Oz", "-o", str(normalized), str(eligible)])
                run(["bcftools", "index", "--csi", str(normalized)])
            else:
                # Strict validation only; do not silently renormalize this declared-normalized input.
                run(["bcftools", "norm", "--check-ref", "e", "-f", str(fasta), "-Ob", "-o", os.devnull, str(source)])
                normalized = source
            normalized_lineage = lineage_counts(normalized)
            if args.normalize_input and normalized_lineage["records"] != partition_counts["eligible_alt_alleles"]:
                raise ValidationError("Split record/ALT lineage count is not conserved; refusing to annotate.")
            annotated = temp / "annotated.vcf.gz"
            run(["bcftools", "annotate", "--pair-logic", "exact", "--annotations", str(clinvar), "--columns", CLINVAR_COLUMNS, "-Oz", "-o", str(annotated), str(normalized)])
            run(["bcftools", "index", "--csi", str(annotated)])
            before_records = lineage_counts(normalized)["records"]
            after_records = compare_gt_semantics(normalized, annotated)
            if before_records != after_records:
                raise ValidationError("Annotation did not preserve record count.")
            if source.stat().st_size != source_stat.st_size or source.stat().st_mtime_ns != source_stat.st_mtime_ns:
                raise ValidationError("Input changed during processing; refusing to publish a result.")
            report = aggregate(annotated)
            report.update({"status": "complete", "format": "clinvar-exact-v1", "normalization": "strict REF check; no REF swaps; biallelic split uses --multi-overlaps . only when requested", "match": "bcftools annotate --pair-logic exact", "source_provenance": {"path": str(source), "size_bytes": source_stat.st_size, "mtime_ns": source_stat.st_mtime_ns}, "lineage_counts": {"pre_normalization_source": source_lineage, "raw_partition": partition_counts, "post_normalization": normalized_lineage, "post_annotation": {"records": after_records}}, "input_records_pre_normalization": source_lineage["records"], "input_records_post_normalization": before_records, "input_records_retained_post_annotation": after_records, "genotype_semantics_preserved": True, "limitations": ["For --normalize-input, non-ACGTN REF records are preserved unchanged in quarantine-original.bcf and excluded from normalization/annotation.", "All remaining records use strict --check-ref e; other REF mismatches fail the run.", "Unmatched means no exact prepared-ClinVar allele, not benign.", "Complete nonreference and partial called-ALT classification strata are reported separately; original genotype/allele mapping is required for inheritance. Annotations remain on all retained records.", "Counts summarize classifications and review status; they are not diagnoses.", "FILTER dot is reported separately from PASS and other FILTER values."]})
            published = temp / "result"
            published.mkdir()
            if args.normalize_input:
                shutil.move(str(normalized), published / "normalized-core.vcf.gz")
                shutil.move(str(normalized) + ".csi", published / "normalized-core.vcf.gz.csi")
                shutil.move(str(quarantine), published / "quarantine-original.bcf")
            shutil.move(str(annotated), published / "annotated.vcf.gz")
            shutil.move(str(annotated) + ".csi", published / "annotated.vcf.gz.csi")
            (published / "aggregate-qc.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            os.replace(published, out)
        print(json.dumps({"status": "complete", "input_retained_records": before_records, "matched_records": report["counts"].get("matched_records", 0)}, sort_keys=True))
        return 0
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
