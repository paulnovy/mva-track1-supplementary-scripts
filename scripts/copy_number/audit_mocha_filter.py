#!/usr/bin/env python3
"""Replay a retained MoChA README filter; no caller, phasing, network or WGS run.

Every one of the 32 gate ablations is checked row-for-row against the actual
retained AWK program, not merely against the uninformative full-filter zero.
The input files are never opened for writing. Output files are exclusive-create.
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import sys


# Exact substrings in the pinned README's "Filter callset" example. Fail closed
# if the supplied program differs; do not silently audit a different filter.
EXPRESSIONS = {
    "S_sample_qc": '!($(g["sample_id"]) in xcl)',
    "C_not_CNP": '$(g["type"])!~"^CNP"',
    "H_spacing_or_concordance": '( $(g["chrom"])~"X" && gender=="M" || bdev<0.1 || $(g["n50_hets"])<2e5 || lod_baf_conc!="nan" && lod_baf_conc>10.0 )',
    "P_precision_or_phase": '( $(g["bdev_se"])!="nan" || lod_baf_phase!="nan" && lod_baf_phase>10.0 )',
    "G_not_germline_like_gain": '( rel_cov<2.1 || bdev<0.05 || len>5e5 && bdev<0.1 && rel_cov<2.5 || len>5e6 && bdev<0.15 )',
}


def parse_table(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def evaluate(row, sample):
    """Independent numeric interpretation for retained finite/literal-nan data."""
    b, rc, length = float(row["bdev"]), float(row["rel_cov"]), int(row["length"])
    return dict(zip(EXPRESSIONS, (
        not (float(sample["call_rate"]) < .97 or float(sample["baf_auto"]) > .03),
        not row["type"].startswith("CNP"),
        ("X" in row["chrom"] and row["computed_gender"] == "M") or b < .1
        or float(row["n50_hets"]) < 2e5
        or (row["lod_baf_conc"] != "nan" and float(row["lod_baf_conc"]) > 10),
        row["bdev_se"] != "nan"
        or (row["lod_baf_phase"] != "nan" and float(row["lod_baf_phase"]) > 10),
        rc < 2.1 or b < .05 or (length > 5e5 and b < .1 and rc < 2.5)
        or (length > 5e6 and b < .15),
    )))


def snapshot(path):
    stat = path.stat()
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "mode": oct(stat.st_mode & 0o777), "sha256": digest}


def bcf_header(path):
    # BCF header only: gzip handles concatenated BGZF; no variant/GT scan.
    with gzip.open(path, "rb") as handle:
        if handle.read(5) != b"BCF\x02\x02":
            raise ValueError(f"not a BCF2.2 header: {path}")
        length = struct.unpack("<I", handle.read(4))[0]
        if length > 10_000_000:
            raise ValueError("unexpected BCF header size")
        return handle.read(length).rstrip(b"\0").decode()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    root, out = args.data_root.resolve(), args.output.resolve()
    if not out.is_relative_to(root / "results"):
        parser.error("output must be a new derived directory under DATA_ROOT/results")
    out.mkdir(mode=0o700, parents=True, exist_ok=True)
    if out.stat().st_mode & 0o077:
        parser.error("output directory must have private permissions")
    ops = root / "operations/mocha-interpretation-20260910"
    src = root / "results/mocha-wgs-20260909"
    old = root / "results/mocha-interpretation-20260910"
    priority = root / "operations/mocha-priority-20260910"
    paths = [src / "result/calls.tsv", src / "result/stats.tsv",
             src / "result/annotated.bcf", src / "result/annotated.bcf.csi",
             src / "work/all.unphased-gc.bcf", src / "work/all.unphased-gc.bcf.csi",
             src / "excluded-variants.bcf", src / "excluded-variants.bcf.csi",
             src / "work/excluded-regions.sorted.bed", src / "summary.json",
             root / "resources/mocha-20260909/cnps.bed",
             ops / "author-filter.awk", ops / "upstream-README.md", ops / "upstream-mocha.c",
             ops / "analyze.py", ops / "baseline.json", ops / "verification.json",
             old / "calls.author-example-filter.tsv", old / "calls.review.tsv", old / "summary.json",
             priority / "run_mocha_priority.sh", priority / "prepare_and_run_mocha_priority.sh",
             priority / "queue.log",
             root / "reports/mocha-review-20260910/MoChA_analiza_naukowa_2026-09-10.md",
             root / "reports/mva1-full-working-20260914/MVA1_WORKING_REPORT_EN.md"]
    # sudo is used only for reading root-owned historical evidence, never for
    # running this program or changing ownership of historical data.
    def read_bytes(path):
        try:
            return path.read_bytes()
        except PermissionError:
            return subprocess.run(["sudo", "-n", "cat", str(path)],
                                  check=True, capture_output=True).stdout

    def source_snapshot(path):
        try:
            return snapshot(path)
        except PermissionError:
            stat = path.stat()
            return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                    "mode": oct(stat.st_mode & 0o777),
                    "sha256": hashlib.sha256(read_bytes(path)).hexdigest()}

    before = [source_snapshot(path) for path in paths]
    calls_bytes, stats_bytes = read_bytes(paths[0]), read_bytes(paths[1])
    calls, stats = parse_table(calls_bytes.decode()), parse_table(stats_bytes.decode())
    samples = {row["sample_id"]: row for row in stats}
    if len(samples) != len(stats) or not calls:
        raise ValueError("empty calls or duplicate sample statistics")
    program = (ops / "author-filter.awk").read_text()
    readme = (ops / "upstream-README.md").read_text()
    block = readme.split("Filter callset\n==============", 1)[1].split("```", 2)[1]
    extracted = block.split('awk -F "\\t" \'', 1)[1].split("' \\\n", 1)[0] + "\n"
    if program != extracted:
        raise ValueError("retained AWK differs from pinned README extraction")
    for expr in EXPRESSIONS.values():
        if program.count(expr) != 1:
            raise ValueError(f"unexpected retained AWK expression: {expr}")
    for remote, local in [("public-pinned-README.md", "upstream-README.md"),
                          ("public-pinned-mocha.c", "upstream-mocha.c")]:
        if (out / remote).read_bytes() != (ops / local).read_bytes():
            raise ValueError(f"public pinned source differs: {local}")

    gates = [evaluate(row, samples[row["sample_id"]]) for row in calls]
    ids = [f"MC{i:02d}" for i in range(1, len(calls) + 1)]
    def replay(awk_program, stats_input=None):
        return subprocess.run(["awk", "-F", "\t", awk_program,
                               "-" if stats_input is not None else str(paths[1]), str(paths[0])],
                              input=stats_input, capture_output=True, text=True, check=True,
                              env={**os.environ, "LC_ALL": "C"}).stdout

    actual = replay(program)
    if actual.encode() != (old / "calls.author-example-filter.tsv").read_bytes():
        raise ValueError("AWK replay differs from retained filtered output")
    ablations = []
    for count in range(len(EXPRESSIONS) + 1):
        for disabled in itertools.combinations(EXPRESSIONS, count):
            active = [name for name in EXPRESSIONS if name not in disabled]
            altered = program
            for name in disabled:
                altered = altered.replace(EXPRESSIONS[name], "(1)")
            keep = [i for i, row_gates in enumerate(gates) if all(row_gates[k] for k in active)]
            if parse_table(replay(altered)) != [calls[i] for i in keep]:
                raise ValueError(f"AWK / independent Python mismatch: {disabled}")
            ablations.append({"disabled": ";".join(disabled) or "none", "retained": len(keep),
                              "retained_ids": ";".join(ids[i] for i in keep)})

    # Positive support under G ablation makes sample-gate checks non-vacuous.
    technical = program.replace(EXPRESSIONS["G_not_germline_like_gain"], "(1)")
    if not parse_table(replay(technical)):
        raise ValueError("no positive control for sample-gate verification")
    controls = []
    for field, value in [("call_rate", "0.96"), ("baf_auto", "0.04")]:
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=list(stats[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({**s, field: value} for s in stats)
        if parse_table(replay(technical, buffer.getvalue())):
            raise ValueError(f"sample gate did not remove positive controls: {field}")
        controls.append({"field": field, "synthetic_value": value, "retained": 0})

    def write(name, content):
        with (out / name).open("xb") as handle:
            handle.write(content.encode() if isinstance(content, str) else content)

    def table(name, rows):
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        write(name, buffer.getvalue())

    old_review = parse_table((old / "calls.review.tsv").read_text())
    historical_fields = ["sample_qc", "not_CNP", "marker_spacing_or_phase_concordance",
                         "bdev_precision_or_strong_phase", "not_germline_like_gain"]
    if len(old_review) != len(calls):
        raise ValueError("historical review call count mismatch")
    for i, row in enumerate(old_review):
        if any(row[k] != v for k, v in calls[i].items()):
            raise ValueError("historical review differs from raw call")
        expected = ";".join(old_name for name, old_name in zip(EXPRESSIONS, historical_fields)
                            if not gates[i][name])
        if expected != row["filter_failed_gates"]:
            raise ValueError("historical per-event failure attribution differs")

    events, waterfall, gate_counts = [], [], []
    for i, row in enumerate(calls):
        failed = [k for k, ok in gates[i].items() if not ok]
        events.append({"review_id": ids[i], **row, **{k: int(v) for k, v in gates[i].items()},
                       "all_failed_gates": ";".join(failed), "first_failed_gate": next(iter(failed), "none"),
                       "masked_fraction_historical": old_review[i]["span_masked_fraction"]})
    surviving = list(range(len(calls)))
    for name in EXPRESSIONS:
        removed = [i for i in surviving if not gates[i][name]]
        kept = [i for i in surviving if gates[i][name]]
        waterfall.append({"gate": name, "before": len(surviving), "excluded": len(removed),
                          "after": len(kept), "excluded_ids": ";".join(ids[i] for i in removed),
                          "retained_ids": ";".join(ids[i] for i in kept)})
        failing = [i for i in range(len(calls)) if not gates[i][name]]
        gate_counts.append({"gate": name, "failing_all_calls": len(failing),
                            "failing_non_CNP": sum(not calls[i]["type"].startswith("CNP") for i in failing),
                            "failing_ids": ";".join(ids[i] for i in failing)})
        surviving = kept
    sample_rows = [{**s, "fail_call_rate_lt_0.97": int(float(s["call_rate"]) < .97),
                    "fail_baf_auto_gt_0.03": int(float(s["baf_auto"]) > .03),
                    "baf_sd_present": int("baf_sd" in s), "baf_sd_filter_used": 0,
                    "baf_corr_filter_used": 0, "cov_sd_filter_used": 0} for s in stats]

    write("calls.unfiltered.tsv", calls_bytes)
    write("stats.unfiltered.tsv", stats_bytes)
    write("author-filter.awk", program)
    write("calls.exact-filter.tsv", actual)
    table("events.tsv", events)
    table("waterfall.tsv", waterfall)
    table("gate-counts.tsv", gate_counts)
    table("ablations.tsv", ablations)
    table("sample-qc.tsv", sample_rows)
    write("retained-queue.log", read_bytes(priority / "queue.log"))
    for name, path in [("input.header.vcf", paths[4]), ("annotated.header.vcf", paths[2])]:
        write(name, bcf_header(path))
    queue = read_bytes(priority / "queue.log").decode()
    # The log contains two passes; count only the pass associated with writing.
    pairs = re.findall(r"Using (\d+) out of (\d+) read variants from contig (\S+)\nWritten (\d+) variants for contig (\S+)", queue)
    if any(c != e or b != d for a, b, c, d, e in pairs):
        raise ValueError("inconsistent retained per-contig log")
    log_counts = [{"chrom": c, "modeled": int(a), "input": int(b), "written": int(d)}
                  for a, b, c, d, e in pairs]
    table("retained-log-counts.tsv", log_counts)
    historical_baseline = {x["path"]: x for x in json.loads((ops / "baseline.json").read_text())}
    matches = {x["path"]: x["sha256"] == historical_baseline[x["path"]]["sha256"]
               for x in before if x["path"] in historical_baseline}
    after = [source_snapshot(path) for path in paths]
    if before != after:
        raise ValueError("historical evidence changed during audit")
    summary = {"utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "command": shlex.join([sys.executable, *sys.argv]), "python": sys.version,
               "awk": subprocess.run(["awk", "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0],
               "script": snapshot(Path(__file__).resolve()), "raw_calls": len(calls),
               "type_counts": dict(collections.Counter(r["type"] for r in calls)),
               "exact_filter_survivors": len(parse_table(actual)),
               "waterfall": waterfall, "all_32_ablations_row_verified": len(ablations) == 32,
               "historical_per_event_attribution_agrees": True, "sample_gate_positive_controls": controls,
               "readme_awk_exact_match": True, "public_pinned_sources_exact_match": True,
               "public_commit": "95686b7b65f53a490513be76bb120e5fc20a8bcf",
               "historical_baseline_sha256_matches": matches,
               "log_modeled_sites": sum(x["modeled"] for x in log_counts),
               "log_written_loci": sum(x["written"] for x in log_counts),
               "historical_sources_unchanged": before == after,
               "mocha_or_phasing_or_wgs_rerun": False, "network_used_by_audit": False,
               "ablation_scope": "retained post-call TSV only; no upstream genotype/CNP/mask ablation"}
    write("summary.json", json.dumps(summary, indent=2) + "\n")
    write("source-manifest.json", json.dumps(before, indent=2) + "\n")
    write("artifact-manifest.json", json.dumps([snapshot(p) for p in sorted(out.iterdir()) if p.is_file()], indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ["raw_calls", "exact_filter_survivors", "all_32_ablations_row_verified", "historical_sources_unchanged"]}))


if __name__ == "__main__":
    main()
