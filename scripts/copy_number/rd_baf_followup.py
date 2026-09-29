#!/usr/bin/env python3
"""Masked, sample-calibrated follow-up of completed 1 Mb RD/BAF summaries.

This is an exploratory screen, not a clinical mosaic-CNV caller. It deliberately
reuses completed MAPQ10 depth and genotype-selected BAF summaries; it never reads
alignments or emits variants.
"""
import argparse
import bisect
import collections
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics


AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))
PAR38 = {
    "chrX": ((10000, 2781479), (155701383, 156030895)),
    "chrY": ((10000, 2781479), (56887903, 57217415)),
}


def median(values):
    return statistics.median(values) if values else None


def quantile(values, p):
    values = sorted(values)
    if not values:
        return None
    pos = p * (len(values) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (pos - lo)


def mad(values, center=None):
    if not values:
        return None
    center = median(values) if center is None else center
    return median([abs(x - center) for x in values])


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_intervals(path, chrom_col, start_col, end_col, allowed=None):
    by_chrom = collections.defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom = fields[chrom_col]
            if allowed is not None and chrom not in allowed:
                continue
            start, end = int(fields[start_col]), int(fields[end_col])
            if end > start:
                by_chrom[chrom].append((start, end))
    for chrom, rows in by_chrom.items():
        rows.sort()
        merged = []
        for start, end in rows:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        by_chrom[chrom] = merged
    return by_chrom


def overlap_bp(intervals, start, end):
    if not intervals or end <= start:
        return 0
    starts = [x[0] for x in intervals]
    idx = max(0, bisect.bisect_left(starts, start) - 1)
    total = 0
    while idx < len(intervals) and intervals[idx][0] < end:
        left, right = intervals[idx]
        total += max(0, min(end, right) - max(start, left))
        idx += 1
    return total


def load_windows(path):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({
                "chrom": row["chrom"], "start": int(row["start"]), "end": int(row["end"]),
                "gc": float(row["gc"]), "depth": float(row["depth_x"]),
                "old_rd": float(row["relative_depth"]), "hets": int(row["het_count"]),
                "bafdev": float(row["median_abs_baf_minus_half"])
                if row["median_abs_baf_minus_half"] != "None" else None,
            })
    return rows


def annotate_masks(rows, umap, superdups, gaps, centromeres):
    for row in rows:
        span = row["end"] - row["start"]
        chrom = row["chrom"]
        row["unique_fraction"] = overlap_bp(umap.get(chrom), row["start"], row["end"]) / span
        row["segdup_fraction"] = overlap_bp(superdups.get(chrom), row["start"], row["end"]) / span
        row["gap_fraction"] = overlap_bp(gaps.get(chrom), row["start"], row["end"]) / span
        row["centromere_fraction"] = overlap_bp(centromeres.get(chrom), row["start"], row["end"]) / span
        reasons = []
        if row["unique_fraction"] < 0.85:
            reasons.append("low_unique_mappability")
        if row["segdup_fraction"] > 0.10:
            reasons.append("segmental_duplication")
        if row["gap_fraction"] > 0.01:
            reasons.append("assembly_gap")
        if row["centromere_fraction"] > 0:
            reasons.append("centromere")
        row["mask_reasons"] = reasons
        row["clean"] = not reasons


def gc_correct(rows):
    clean = [r for r in rows if r["clean"]]
    if len(clean) < 500:
        raise RuntimeError("too_few_clean_autosomal_bins")
    for row in clean:
        peers = [x for x in clean if x["chrom"] != row["chrom"]]
        peers.sort(key=lambda x: abs(x["gc"] - row["gc"]))
        expected = median([x["depth"] for x in peers[:200]])
        row["rd"] = row["depth"] / expected
    center = median([r["rd"] for r in clean])
    for row in clean:
        row["rd"] /= center
    return center


def contiguous_groups(rows, n):
    by_chrom = collections.defaultdict(list)
    for row in rows:
        if row["clean"]:
            by_chrom[row["chrom"]].append(row)
    groups = []
    for chrom_rows in by_chrom.values():
        chrom_rows.sort(key=lambda r: r["start"])
        run = []
        for row in chrom_rows:
            if run and row["start"] - run[-1]["end"] > 100000:
                run = []
            run.append(row)
            if len(run) >= n:
                groups.append(run[-n:])
    return groups


def empirical_threshold(groups, key, center, p=0.995):
    values = [abs(sum(r[key] for r in group) / len(group) - center) for group in groups]
    return max(quantile(values, p) or 0, 0.01 if key == "rd" else 0.002)


def power_from_observed(groups, key, center, threshold, shift, direction):
    values = [sum(r[key] for r in group) / len(group) - center for group in groups]
    if not values:
        return None
    if direction == "positive":
        return sum(v + shift > threshold for v in values) / len(values)
    return sum(v + shift < -threshold for v in values) / len(values)


def calibration(rows, rd_center, baf_center):
    result = {"method": "empirical_contiguous-bin_signal_injection_v1", "alpha_per_test": 0.01,
              "clinical_lod": False, "scenarios": []}
    for n in (3, 4, 5, 10, 20):
        groups = contiguous_groups(rows, n)
        rd_thr = empirical_threshold(groups, "rd", rd_center)
        baf_groups = [g for g in groups if all(r["bafdev"] is not None for r in g)]
        baf_thr = empirical_threshold(baf_groups, "bafdev", baf_center) if baf_groups else None
        for fraction in (0.05, 0.10, 0.20, 0.30, 0.50):
            loss_baf = fraction / (4 - 2 * fraction)
            gain_baf = fraction / (4 + 2 * fraction)
            result["scenarios"].extend([
                {"bins": n, "event": "mosaic_loss", "fraction": fraction,
                 "expected_rd_shift": -fraction / 2, "expected_mirrored_baf_shift": loss_baf,
                 "rd_detection_power": power_from_observed(groups, "rd", rd_center, rd_thr,
                                                            -fraction / 2, "negative")},
                {"bins": n, "event": "mosaic_gain", "fraction": fraction,
                 "expected_rd_shift": fraction / 2, "expected_mirrored_baf_shift": gain_baf,
                 "rd_detection_power": power_from_observed(groups, "rd", rd_center, rd_thr,
                                                            fraction / 2, "positive")},
                {"bins": n, "event": "copy_neutral_loh", "fraction": fraction,
                 "expected_rd_shift": 0, "expected_mirrored_baf_shift": fraction / 2,
                 "baf_detection_power": power_from_observed(baf_groups, "bafdev", baf_center,
                                                             baf_thr, fraction / 2, "positive")
                 if baf_thr is not None else None},
            ])
        result.setdefault("thresholds", []).append({"bins": n, "rd_abs_mean_shift": rd_thr,
                                                     "baf_abs_mean_shift": baf_thr,
                                                     "placements": len(groups)})
    result["limitations"] = [
        "Power is in-silico signal injection into this sample's retained 1 Mb bins, not a validated clinical LOD.",
        "Observed bins may contain biology as well as technical noise; adjacent bins are not independent.",
        "BAF injection assumes simple one-clone allelic states and is biased by genotype-selected heterozygotes.",
        "Distinct clones, tissue restriction, tumour admixture, treatment effects, and opposing gains/losses are not modeled.",
    ]
    return result


def screen_segments(rows, rd_center, baf_center, rd_noise, baf_noise):
    by_chrom = collections.defaultdict(list)
    for row in rows:
        if row["clean"]:
            by_chrom[row["chrom"]].append(row)
    candidates = []
    rd_cut = max(0.03, 3 * rd_noise / math.sqrt(5))
    baf_cut = max(0.005, 3 * baf_noise / math.sqrt(5))
    for chrom, chrom_rows in by_chrom.items():
        chrom_rows.sort(key=lambda r: r["start"])
        flags = []
        for i in range(len(chrom_rows)):
            local = chrom_rows[max(0, i - 2):min(len(chrom_rows), i + 3)]
            local = [r for r in local if abs(r["start"] - chrom_rows[i]["start"]) <= 2200000]
            rd_shift = median([r["rd"] for r in local]) - rd_center
            baf_values = [r["bafdev"] for r in local if r["bafdev"] is not None]
            baf_shift = median(baf_values) - baf_center if baf_values else 0
            if rd_shift < -rd_cut:
                flag = "loss_like"
            elif rd_shift > rd_cut:
                flag = "gain_like"
            elif baf_shift > baf_cut:
                flag = "copy_neutral_loh_like"
            else:
                flag = None
            flags.append(flag)
        start = 0
        while start < len(chrom_rows):
            flag = flags[start]
            end = start + 1
            while (end < len(chrom_rows) and flags[end] == flag and
                   chrom_rows[end]["start"] - chrom_rows[end - 1]["end"] <= 100000):
                end += 1
            segment = chrom_rows[start:end]
            if flag and len(segment) >= 3:
                candidates.append({
                    "chrom": chrom, "start": segment[0]["start"], "end": segment[-1]["end"],
                    "bins": len(segment), "screen_label": flag,
                    "median_rd": median([r["rd"] for r in segment]),
                    "median_bafdev": median([r["bafdev"] for r in segment if r["bafdev"] is not None]),
                    "het_count": sum(r["hets"] for r in segment),
                })
            start = end
    # Review every screened run against an empirical genome-wide threshold of
    # the same length. This separates discovery flags from stronger follow-up
    # evidence without converting this script into a CNV caller.
    threshold_cache = {}
    for segment in candidates:
        n = segment["bins"]
        if n not in threshold_cache:
            groups = contiguous_groups(rows, n)
            baf_groups = [g for g in groups if all(r["bafdev"] is not None for r in g)]
            threshold_cache[n] = (
                empirical_threshold(groups, "rd", rd_center),
                empirical_threshold(baf_groups, "bafdev", baf_center) if baf_groups else None,
            )
        rd_threshold, baf_threshold = threshold_cache[n]
        source = [r for r in by_chrom[segment["chrom"]]
                  if segment["start"] <= r["start"] and r["end"] <= segment["end"]]
        rd_mean = sum(r["rd"] for r in source) / len(source)
        baf_source = [r["bafdev"] for r in source if r["bafdev"] is not None]
        baf_mean = sum(baf_source) / len(baf_source) if baf_source else None
        segment.update({
            "mean_rd": rd_mean, "mean_bafdev": baf_mean,
            "empirical_rd_threshold": rd_threshold,
            "empirical_baf_threshold": baf_threshold,
            "passes_empirical_rd": abs(rd_mean - rd_center) > rd_threshold,
            "passes_empirical_baf": (baf_mean - baf_center > baf_threshold)
            if baf_mean is not None and baf_threshold is not None else False,
        })
        segment["joint_empirical_support"] = (segment["passes_empirical_rd"] and
                                               segment["passes_empirical_baf"])
    return candidates, {"rd_5bin_cut": rd_cut, "baf_5bin_cut": baf_cut,
                        "empirical_review_quantile_two_sided": 0.995}


def sex_coverage(regions_path, masks, autosomal_depth_center):
    gaps, centromeres = masks
    categories = collections.defaultdict(list)
    with gzip.open(regions_path, "rt") as handle:
        for line in handle:
            fields = line.split()
            chrom = fields[0]
            if chrom not in ("chrX", "chrY"):
                continue
            start, end, depth = int(fields[1]), int(fields[2]), float(fields[3])
            span = end - start
            if span < 500000:
                continue
            masked = overlap_bp(gaps.get(chrom), start, end) + overlap_bp(centromeres.get(chrom), start, end)
            if masked / span > 0.05:
                continue
            mid = (start + end) // 2
            in_par = any(a <= mid < b for a, b in PAR38[chrom])
            categories[(chrom, "PAR" if in_par else "nonPAR")].append(depth)
    ratios = {f"{chrom}_{part}": median(values) / autosomal_depth_center
              for (chrom, part), values in categories.items() if values}
    counts = {f"{chrom}_{part}": len(values) for (chrom, part), values in categories.items()}
    x = ratios.get("chrX_nonPAR")
    y = ratios.get("chrY_nonPAR")
    if x is not None and y is not None and 0.35 <= x <= 0.65 and y >= 0.20:
        pattern = "XY-like coverage"
    elif x is not None and y is not None and 0.80 <= x <= 1.20 and y < 0.10:
        pattern = "XX-like coverage"
    else:
        pattern = "unresolved coverage ploidy pattern"
    return {"ratios_to_clean_autosomal_depth": ratios, "retained_tile_counts": counts,
            "descriptive_pattern": pattern,
            "interpretation": "PAR and non-PAR were separated using GRCh38 boundaries; absent PAR estimates mean no qualifying 1 Mb tile was fully represented. This is coverage QC, not sex/gender metadata or a sex-chromosome mosaic call."}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior-windows", required=True)
    p.add_argument("--prior-summary", required=True)
    p.add_argument("--regions", required=True)
    p.add_argument("--umap-bed", required=True)
    p.add_argument("--superdups", required=True)
    p.add_argument("--gaps", required=True)
    p.add_argument("--centromeres", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    os.umask(0o077)
    out = Path(args.out)
    out.mkdir(mode=0o700, parents=True, exist_ok=True)
    summary_path = out / "summary.json"
    if summary_path.exists():
        raise SystemExit("Existing summary; inspect checkpoint, do not repeat")

    allowed = set(AUTOSOMES) | {"chrX", "chrY"}
    umap = load_intervals(args.umap_bed, 0, 1, 2, allowed)
    superdups = load_intervals(args.superdups, 1, 2, 3, allowed)
    gaps = load_intervals(args.gaps, 1, 2, 3, allowed)
    centromeres = load_intervals(args.centromeres, 1, 2, 3, allowed)
    rows = load_windows(args.prior_windows)
    annotate_masks(rows, umap, superdups, gaps, centromeres)
    raw_depth_center = median([r["depth"] for r in rows if r["clean"]])
    gc_correct(rows)
    clean = [r for r in rows if r["clean"]]
    rd_center = median([r["rd"] for r in clean])
    rd_noise = 1.4826 * mad([r["rd"] for r in clean], rd_center)
    baf_values = [r["bafdev"] for r in clean if r["bafdev"] is not None]
    baf_center = median(baf_values)
    baf_noise = 1.4826 * mad(baf_values, baf_center)

    low = [r for r in rows if r["old_rd"] < 0.8]
    low_reason_counts = collections.Counter(reason for r in low for reason in r["mask_reasons"])
    clean_low = [r for r in low if r["clean"]]
    segments, segment_thresholds = screen_segments(rows, rd_center, baf_center, rd_noise, baf_noise)
    whole = {}
    for chrom in AUTOSOMES:
        chrom_rows = [r for r in clean if r["chrom"] == chrom]
        whole[chrom] = {"clean_bins": len(chrom_rows),
                        "median_rd": median([r["rd"] for r in chrom_rows]),
                        "median_bafdev": median([r["bafdev"] for r in chrom_rows if r["bafdev"] is not None])}
    whole_rd_threshold = 0.03
    whole_candidates = [chrom for chrom, values in whole.items()
                        if abs(values["median_rd"] - rd_center) > whole_rd_threshold]
    max_whole_shift = max(abs(values["median_rd"] - rd_center) for values in whole.values())

    result = {
        "status": "complete",
        "method": "masked_completed_1mb_rd_baf_followup_v1",
        "scope": "Exploratory aggregate screen reusing completed MAPQ10 1 Mb depth and genotype-selected heterozygote BAF; no BAM reread.",
        "mask_policy": {"minimum_umap_k100_unique_fraction": 0.85,
                        "maximum_segmental_duplication_fraction": 0.10,
                        "maximum_gap_fraction": 0.01,
                        "exclude_any_centromere_overlap": True},
        "bins": {"prior_autosomal": len(rows), "retained_clean": len(clean),
                 "excluded": len(rows) - len(clean),
                 "old_rd_below_0_8": len(low), "old_low_rd_mask_flagged": len(low) - len(clean_low),
                 "old_low_rd_still_clean": len(clean_low),
                 "old_low_rd_mask_reasons_nonexclusive": dict(low_reason_counts)},
        "sample_noise": {"rd_center": rd_center, "rd_mad_sigma": rd_noise,
                         "bafdev_center": baf_center, "bafdev_mad_sigma": baf_noise},
        "whole_chromosome": whole,
        "whole_chromosome_screen": {
            "absolute_median_rd_threshold": whole_rd_threshold,
            "candidate_count": len(whole_candidates),
            "maximum_observed_absolute_median_rd_shift": max_whole_shift,
            "simple_single_clone_fraction_at_threshold": 2 * whole_rd_threshold,
            "model": "For a simple disomic-to-monosomic/trisomic mixture, expected RD shift is mosaic fraction / 2; this is an exploratory model, not a clinical LOD.",
        },
        "candidate_segments": {"count": len(segments), "screen_thresholds": segment_thresholds,
                               "labels": dict(collections.Counter(x["screen_label"] for x in segments)),
                               "passing_empirical_rd": sum(x["passes_empirical_rd"] for x in segments),
                               "passing_empirical_baf": sum(x["passes_empirical_baf"] for x in segments),
                               "joint_empirical_support": sum(x["joint_empirical_support"] for x in segments),
                               "max_bins": max((x["bins"] for x in segments), default=0)},
        "sex_chromosomes": sex_coverage(args.regions, (gaps, centromeres), raw_depth_center),
        "calibration": calibration(rows, rd_center, baf_center),
        "inputs": {name: {"path": value, "sha256": sha256(value)} for name, value in {
            "prior_windows": args.prior_windows, "prior_summary": args.prior_summary,
            "umap_k100_unique": args.umap_bed, "genomic_superdups": args.superdups,
            "assembly_gaps": args.gaps, "centromeres": args.centromeres}.items()},
        "resource_sources": {
            "umap_k100_unique": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Unique.Mappability.bb",
            "genomic_superdups": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz",
            "assembly_gaps": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz",
            "centromeres": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz",
        },
        "interpretation": {
            "low_window_review": "Mask overlap distinguishes likely reference/mapping artifacts from clean low-depth bins; isolated clean bins are not CNV calls.",
            "segments": "Candidate labels are robust 5-bin descriptive screening labels, not CNV genotypes or clinical calls.",
            "whole_chromosome": "Chromosome medians are normalized across this sample without a matched reference cohort.",
            "negative_result": "A flat bulk profile cannot exclude low-fraction, tissue-restricted, mixed, or opposing mosaic clones.",
        },
        "review_conclusion": {
            "old_low_rd_windows": f"Of {len(low)} input RD<0.8 bins, {len(low)-len(clean_low)} are mask-flagged and {len(clean_low)} remain clean; clean bins are not validated CNV calls.",
            "whole_chromosome": f"{len(whole_candidates)} autosomes exceed the absolute median-RD screen threshold {whole_rd_threshold:.4f}.",
            "subchromosomal": (f"{len(segments)} short 3-5-bin discovery flags remain; "
                                f"{sum(x['passes_empirical_rd'] for x in segments)} pass matched-length empirical RD review, "
                                f"{sum(x['passes_empirical_baf'] for x in segments)} pass empirical BAF review, and "
                                f"{sum(x['joint_empirical_support'] for x in segments)} have both. "
                                "Any joint flag requires orthogonal/local review and is not a CNV call."),
        },
        "tool_preflight": {
            "MoChA": "Not executed by this script; a separate MoChA workflow requires phased GT plus AD and reference-based phasing/common sites.",
            "CNVpytor": "Not executed by this script; provision separately when a fresh BAM-wide RD+BAF extraction is required.",
            "sources": ["https://github.com/freeseek/mocha", "https://github.com/abyzovlab/CNVpytor", "https://pmc.ncbi.nlm.nih.gov/articles/PMC8612020/"],
        },
        "limitations": [
            "No matched normal/reference cohort, clinical validation, absolute copy-number model, or constitutional/tumour tissue provenance.",
            "BAF comes from variant-only genotype-selected heterozygotes and remains ascertainment-biased despite mask filtering.",
            "1 Mb bins have limited breakpoint resolution and completed MAPQ10 depth is reused.",
            "This script does not establish DNA tissue or timing relative to treatment/transfusion.",
            "This script does not determine whether MoChA or CNVpytor was run elsewhere; review their separate stage records.",
        ],
    }

    with (out / "bins.masked.tsv").open("w") as handle:
        handle.write("chrom\tstart\tend\tunique_fraction\tsegdup_fraction\tgap_fraction\tcentromere_fraction\tclean\tmask_reasons\tcorrected_rd\thet_count\tmedian_abs_baf_minus_half\n")
        for r in rows:
            handle.write("\t".join(map(str, [r["chrom"], r["start"], r["end"], r["unique_fraction"],
                r["segdup_fraction"], r["gap_fraction"], r["centromere_fraction"], int(r["clean"]),
                ",".join(r["mask_reasons"]) or ".", r.get("rd", "."), r["hets"], r["bafdev"]])) + "\n")
    with (out / "segments.tsv").open("w") as handle:
        handle.write("chrom\tstart\tend\tbins\tscreen_label\tmedian_rd\tmedian_abs_baf_minus_half\thet_count\tmean_rd\tmean_bafdev\tempirical_rd_threshold\tempirical_baf_threshold\tpasses_empirical_rd\tpasses_empirical_baf\tjoint_empirical_support\n")
        for x in segments:
            handle.write("\t".join(str(x[k]) for k in ("chrom", "start", "end", "bins", "screen_label", "median_rd", "median_bafdev", "het_count", "mean_rd", "mean_bafdev", "empirical_rd_threshold", "empirical_baf_threshold", "passes_empirical_rd", "passes_empirical_baf", "joint_empirical_support")) + "\n")
    tmp = out / "summary.json.tmp"
    tmp.write_text(json.dumps(result, indent=2) + "\n")
    tmp.replace(summary_path)
    print(json.dumps({"status": "complete", "bins": result["bins"],
                      "candidate_segment_count": len(segments)}))


if __name__ == "__main__":
    main()
