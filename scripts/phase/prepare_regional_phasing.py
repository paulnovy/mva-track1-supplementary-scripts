#!/usr/bin/env python3
"""Prepare and read-phase a deliberately small, regional VCF with WhatsHap.

This tool never realigns reads or changes its source BAM/VCF. It writes only a
new regional candidate VCF and, when requested, a WhatsHap output VCF.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pysam


class ValidationError(RuntimeError):
    """Raised before any phasing work is started."""


@dataclass(frozen=True)
class Region:
    contig: str
    start: int  # one-based inclusive
    end: int  # one-based inclusive

    @property
    def label(self) -> str:
        return f"{self.contig}:{self.start}-{self.end}"


def die(message: str) -> None:
    raise ValidationError(message)


def parse_region(value: str) -> Region:
    try:
        contig, interval = value.rsplit(":", 1)
        start_text, end_text = interval.replace(",", "").split("-", 1)
        start, end = int(start_text), int(end_text)
    except (ValueError, AttributeError):
        die(f"Invalid region {value!r}; expected CONTIG:START-END (1-based inclusive).")
    if not contig or start < 1 or end < start:
        die(f"Invalid region {value!r}; expected positive START <= END.")
    return Region(contig, start, end)


def parse_target(value: str) -> tuple[str, int, str | None, str | None]:
    try:
        fields = value.split(":")
        if len(fields) not in (2, 4):
            raise ValueError
        contig, position = fields[:2]
        pos = int(position.replace(",", ""))
    except (ValueError, AttributeError):
        die(f"Invalid target {value!r}; expected CONTIG:POSITION (1-based).")
    if not contig or pos < 1:
        die(f"Invalid target {value!r}; expected a positive position.")
    if len(fields) == 2:
        return contig, pos, None, None
    ref, alt = fields[2:]
    if not ref or not alt or any(base not in "ACGT" for base in ref + alt):
        die(f"Invalid target {value!r}; REF and ALT must be A/C/G/T alleles.")
    return contig, pos, ref, alt


def existing_index(path: Path, kind: str) -> Path:
    if kind == "bam":
        choices = [Path(f"{path}.bai"), Path(f"{path}.csi"), path.with_suffix(".bai"), path.with_suffix(".csi")]
    else:
        choices = [Path(f"{path}.tbi"), Path(f"{path}.csi")]
    for candidate in choices:
        if candidate.is_file():
            return candidate
    die(f"Missing {kind.upper()} index for input; create an index before running this tool.")


def validate_bam(path: Path) -> pysam.AlignmentFile:
    if not path.is_file():
        die("BAM input is not a regular file.")
    index = existing_index(path, "bam")
    if index.stat().st_mtime < path.stat().st_mtime:
        die("BAM index is older than BAM; rebuild the index after BAM finalization.")
    try:
        errors = pysam.quickcheck(str(path))
    except Exception as exc:  # pysam's error varies by htslib version
        die(f"BAM quickcheck could not validate the BAM ({type(exc).__name__}).")
    if errors:
        die("BAM quickcheck failed; do not phase an incomplete or corrupt BAM.")
    try:
        bam = pysam.AlignmentFile(str(path), "rb")
        if bam.header.get("HD", {}).get("SO") != "coordinate":
            bam.close()
            die("BAM header is not coordinate-sorted (HD:SO must be coordinate).")
        if not bam.has_index():
            bam.close()
            die("BAM index cannot be opened by htslib.")
        return bam
    except ValidationError:
        raise
    except Exception as exc:
        die(f"BAM/index cannot be opened by htslib ({type(exc).__name__}).")


def validate_vcf(path: Path) -> pysam.VariantFile:
    if not path.is_file() or not (str(path).endswith(".vcf.gz") or str(path).endswith(".vcf.bgz")):
        die("VCF must be an existing bgzip-compressed .vcf.gz or .vcf.bgz file.")
    existing_index(path, "vcf")
    try:
        vcf = pysam.VariantFile(str(path))
        return vcf
    except Exception as exc:
        die(f"VCF/index cannot be opened by htslib ({type(exc).__name__}).")


def choose_sample(vcf: pysam.VariantFile, requested: str | None) -> str:
    samples = list(vcf.header.samples)
    if requested:
        if requested not in vcf.header.samples:
            die("Requested VCF sample is absent.")
        return requested
    if len(samples) != 1:
        die("VCF has multiple samples; select exactly one with --sample.")
    if not samples:
        die("VCF contains no sample column.")
    return samples[0]


def validate_read_group_sample(bam: pysam.AlignmentFile, sample: str, ignore_read_groups: bool) -> None:
    sample_names = {group.get("SM") for group in bam.header.get("RG", []) if group.get("SM")}
    if ignore_read_groups:
        if len(sample_names) != 1:
            die("--ignore-read-groups requires exactly one named BAM read-group sample.")
        return
    if sample not in sample_names:
        die("Selected VCF sample does not match a BAM read-group SM; use --ignore-read-groups only after confirming sample identity.")


def validate_reference(path: Path, bam: pysam.AlignmentFile, regions: Iterable[Region]) -> None:
    if not path.is_file() or not Path(f"{path}.fai").is_file():
        die("Reference FASTA and its .fai index are both required.")
    try:
        fasta = pysam.FastaFile(str(path))
    except Exception as exc:
        die(f"Reference FASTA/index cannot be opened ({type(exc).__name__}).")
    try:
        bam_lengths = dict(zip(bam.references, bam.lengths))
        fasta_lengths = dict(zip(fasta.references, fasta.lengths))
        for region in regions:
            if region.contig not in fasta_lengths:
                die("A requested region contig is absent from the reference FASTA.")
            if region.contig not in bam_lengths:
                die("A requested region contig is absent from the BAM header.")
            if fasta_lengths[region.contig] != bam_lengths[region.contig]:
                die("Reference FASTA and BAM disagree on a requested contig length.")
            if region.end > fasta_lengths[region.contig]:
                die("A requested region exceeds the reference contig length.")
    finally:
        fasta.close()


def depth_from_call(call, allow_ad: bool) -> tuple[int | None, str | None]:
    """Return FORMAT/DP, or explicitly opted-in sum of valid FORMAT/AD."""
    dp = call.get("DP")
    if dp is not None:
        return dp, "format_dp"
    ad = call.get("AD")
    if allow_ad and ad is not None and all(isinstance(value, int) and value >= 0 for value in ad):
        return sum(ad), "sum_format_ad"
    return None, None


def eligible_reason(record: pysam.VariantRecord, sample: str, min_qual: float | None, min_dp: int | None, min_gq: int | None, min_pl_delta: int | None, min_ad: int | None, dp_from_ad: bool) -> tuple[str | None, str | None]:
    gt = record.samples[sample].get("GT")
    if gt is None or tuple(gt) not in ((0, 1), (1, 0)):
        return "gt_not_called_diploid_ref_alt", None
    filters = set(record.filter.keys())
    if filters and filters != {"PASS"}:
        return "filter_not_pass", None
    if min_qual is not None and (record.qual is None or record.qual < min_qual):
        return "qual_below_minimum", None
    call = record.samples[sample]
    depth, depth_method = depth_from_call(call, dp_from_ad)
    if min_dp is not None and (depth is None or depth < min_dp):
        return "dp_below_minimum_or_missing", depth_method
    gq = call.get("GQ")
    pl = call.get("PL")
    # For a biallelic diploid 0/1 call, PL[1] must itself be the best
    # likelihood. A generic best-vs-second gap could instead support 0/0/1/1.
    pl_het_gap = None
    if pl is not None and len(pl) == 3 and all(value is not None for value in pl) and pl[1] == min(pl):
        pl_het_gap = min(pl[0], pl[2]) - pl[1]
    if min_gq is not None and (gq is None or gq < min_gq) and (min_pl_delta is None or pl_het_gap is None or pl_het_gap < min_pl_delta):
        return "gq_or_pl_below_minimum_or_missing", depth_method
    if min_ad is not None and call.get("AD") is not None:
        ad = call.get("AD")
        if len(ad) < 2 or ad[0] is None or ad[1] is None or ad[0] < min_ad or ad[1] < min_ad:
            return "ad_below_minimum", depth_method
    if len(record.alts or ()) != 1:
        return "not_biallelic", depth_method
    ref, alt = record.ref, record.alts[0]
    if not ref or not alt or any(base not in "ACGT" for base in ref + alt):
        return "alleles_not_acgt", depth_method
    if len(ref) == len(alt) == 1 and ref != alt:
        return None, depth_method
    if len(ref) != len(alt):
        return None, depth_method
    return "not_snv_or_indel", depth_method


def one_sample_header(source: pysam.VariantHeader, sample: str) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    for line in str(source).splitlines():
        if line.startswith("##"):
            header.add_line(line)
    header.add_sample(sample)
    return header


def copy_selected_record(record: pysam.VariantRecord, header: pysam.VariantHeader, sample: str) -> pysam.VariantRecord:
    copied = header.new_record(
        contig=record.contig,
        start=record.start,
        stop=record.stop,
        id=record.id,
        alleles=record.alleles,
        qual=record.qual,
    )
    for name in record.filter.keys():
        copied.filter.add(name)
    for key, value in record.info.items():
        copied.info[key] = value
    for key, value in record.samples[sample].items():
        copied.samples[sample][key] = value
    copied.samples[sample].phased = record.samples[sample].phased
    return copied


def new_output_directory(path: Path) -> None:
    if path.exists():
        die("Output directory already exists; choose a new run directory to avoid overwriting/repeating output.")
    if not path.parent.is_dir():
        die("Parent of output directory does not exist.")


def atomic_prepare(args: argparse.Namespace) -> None:
    out_dir = Path(args.out).resolve()
    new_output_directory(out_dir)
    bam_path, vcf_path, fasta_path = (Path(args.bam).resolve(), Path(args.vcf).resolve(), Path(args.fasta).resolve())
    regions = [parse_region(item) for item in args.region]
    targets = {parse_target(item) for item in args.target}
    bam = validate_bam(bam_path)
    vcf = validate_vcf(vcf_path)
    try:
        sample = choose_sample(vcf, args.sample)
        validate_read_group_sample(bam, sample, args.ignore_read_groups)
        validate_reference(fasta_path, bam, regions)
        tmp_dir = Path(tempfile.mkdtemp(prefix=f".{out_dir.name}.preparing-", dir=out_dir.parent))
        try:
            plain_vcf = tmp_dir / "candidates.vcf"
            header = one_sample_header(vcf.header, sample)
            counts: Counter[str] = Counter()
            depth_method_counts: Counter[str] = Counter()
            selected_records: list[pysam.VariantRecord] = []
            seen: set[tuple[str, int, str, tuple[str, ...] | None]] = set()
            with pysam.VariantFile(str(plain_vcf), "w", header=header) as output:
                for region in regions:
                    try:
                        records = vcf.fetch(region.contig, region.start - 1, region.end)
                        for record in records:
                            identity = (record.contig, record.pos, record.ref, record.alts)
                            if identity in seen:
                                continue
                            seen.add(identity)
                            reason, depth_method = eligible_reason(record, sample, args.min_qual, args.min_dp, args.min_gq, args.min_pl_delta, args.min_ad, args.dp_from_ad)
                            if reason:
                                counts[reason] += 1
                                continue
                            depth_method_counts[depth_method or "not_requested"] += 1
                            selected_records.append(copy_selected_record(record, header, sample))
                    except ValueError:
                        die("A requested region contig is absent from the indexed VCF.")
                contig_rank = {contig: rank for rank, contig in enumerate(vcf.header.contigs)}
                for record in sorted(selected_records, key=lambda item: (contig_rank[item.contig], item.start, item.stop, item.ref, item.alts)):
                    output.write(record)
            selected = len(selected_records)
            compressed = Path(pysam.tabix_index(str(plain_vcf), preset="vcf", force=True))
            manifest = {
                "format": "regional-phasing-prep-v1",
                "source_inputs": {"bam": str(bam_path), "vcf": str(vcf_path), "fasta": str(fasta_path)},
                "regions_1based_inclusive": [region.label for region in regions],
                "targets_1based": [f"{chrom}:{pos}" if ref is None else f"{chrom}:{pos}:{ref}:{alt}" for chrom, pos, ref, alt in sorted(targets)],
                "ignore_read_groups": args.ignore_read_groups,
                "min_qual": args.min_qual,
                "min_dp": args.min_dp, "min_gq": args.min_gq, "min_pl_delta": args.min_pl_delta, "min_ad": args.min_ad,
                "dp_from_ad": args.dp_from_ad,
                "selected_depth_method_counts": dict(sorted(depth_method_counts.items())),
                "candidate_count": selected,
                "excluded_counts": dict(sorted(counts.items())),
                "selection": "exact 0/1, biallelic A/C/G/T simple SNVs/indels with PASS or absent FILTER only; requested QUAL/DP/GQ-or-PL/AD thresholds are recorded above",
                "original_vcf_preserved": True,
            }
            (tmp_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
            (tmp_dir / "README-run-state.txt").write_text(
                "Prepared successfully. Run the phase subcommand once; it will refuse to overwrite a completed phase output.\n"
            )
            os.replace(tmp_dir, out_dir)
        except Exception:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            raise
    finally:
        vcf.close()
        bam.close()
    print(f"prepared candidates={selected} excluded={sum(counts.values())} regions={len(regions)}")


def load_manifest(out_dir: Path) -> dict:
    path = out_dir / "manifest.json"
    if not path.is_file() or not (out_dir / "candidates.vcf.gz").is_file() or not Path(f"{out_dir / 'candidates.vcf.gz'}.tbi").is_file():
        die("Prepared regional inputs are incomplete; run prepare into a new output directory first.")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        die("Prepared manifest is invalid.")


def summarize_phase(phased_vcf: Path, targets: set[tuple[str, int, str | None, str | None]]) -> dict:
    blocks: dict[tuple[str, str], int] = defaultdict(int)
    target_blocks: dict[tuple[str, int, str | None, str | None], tuple[tuple[str, str] | None, tuple[int | None, ...] | None]] = {}
    phased_count = 0
    total = 0
    with pysam.VariantFile(str(phased_vcf)) as vcf:
        samples = list(vcf.header.samples)
        if len(samples) != 1:
            die("Regional candidate VCF must contain exactly one sample.")
        sample = samples[0]
        for record in vcf:
            total += 1
            call = record.samples[sample]
            ps = call.get("PS")
            key: tuple[str, str] | None = None
            if call.phased and ps is not None:
                key = (record.contig, str(ps))
                blocks[key] += 1
                phased_count += 1
            for target in targets:
                contig, pos, ref, alt = target
                if record.contig == contig and record.pos == pos and (ref is None or (record.ref == ref and record.alts == (alt,))):
                    target_blocks[target] = (key, call.get("GT"))
    connected = sum(1 for key, _ in target_blocks.values() if key is not None and blocks[key] >= 2)
    pair_same_ps = False
    pair_orientation = None
    if len(targets) == 2 and len(target_blocks) == 2:
        first, second = (target_blocks[target] for target in sorted(targets))
        if first[0] is not None and first[0] == second[0]:
            pair_same_ps = True
            if first[1] in ((0, 1), (1, 0)) and second[1] in ((0, 1), (1, 0)):
                pair_orientation = "cis" if first[1] == second[1] else "trans"
    return {
        "candidate_records": total,
        "phased_records": phased_count,
        "unphased_records": total - phased_count,
        "phase_block_count": len(blocks),
        "phase_blocks": [
            {"contig": contig, "phase_set": ps, "records": count}
            for (contig, ps), count in sorted(blocks.items())
        ],
        "targets": {
            "requested": len(targets),
            "found_in_candidate_vcf": len(target_blocks),
            "connected_to_a_block_with_another_candidate": connected,
            "not_connected": len(target_blocks) - connected,
            "targets_same_phase_set": pair_same_ps,
            "target_pair_cis_or_trans": pair_orientation,
        },
        "interpretation": "A phase relationship is comparable only within the same contig and PS. This summary does not infer cis/trans across phase sets or confidence.",
    }


def phase(args: argparse.Namespace) -> None:
    out_dir = Path(args.out).resolve()
    if not out_dir.is_dir():
        die("Prepared output directory does not exist.")
    if (out_dir / "phase.complete.json").exists() or (out_dir / "phased.vcf.gz").exists():
        die("Completed or partial phased output already exists; refusing to overwrite or repeat it.")
    manifest = load_manifest(out_dir)
    bam_path = Path(manifest["source_inputs"]["bam"])
    fasta_path = Path(manifest["source_inputs"]["fasta"])
    regions = [parse_region(item) for item in manifest["regions_1based_inclusive"]]
    bam = validate_bam(bam_path)
    try:
        validate_reference(fasta_path, bam, regions)
        candidates = validate_vcf(out_dir / "candidates.vcf.gz")
        try:
            validate_read_group_sample(bam, choose_sample(candidates, None), manifest.get("ignore_read_groups", False))
        finally:
            candidates.close()
    finally:
        bam.close()
    targets = {parse_target(item) for item in manifest.get("targets_1based", [])}
    command = shlex.split(args.whatshap)
    if not command:
        die("--whatshap must name an executable or command.")
    with tempfile.TemporaryDirectory(prefix=".phasing-", dir=out_dir) as temp:
        temp_dir = Path(temp)
        temp_output = temp_dir / "phased.vcf.gz"
        command += [
            "phase", "--tag=PS", "--indels", "--reference", str(fasta_path), "--output", str(temp_output),
            str(out_dir / "candidates.vcf.gz"), str(bam_path),
        ]
        if manifest.get("ignore_read_groups", False):
            command.insert(command.index("phase") + 1, "--ignore-read-groups")
        with (temp_dir / "whatshap.stdout.log").open("w") as stdout, (temp_dir / "whatshap.stderr.log").open("w") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0 or not temp_output.is_file():
            for name in ("whatshap.stdout.log", "whatshap.stderr.log"):
                os.replace(temp_dir / name, out_dir / name)
            die("WhatsHap did not complete successfully; no phased output was published.")
        pysam.tabix_index(str(temp_output), preset="vcf", force=True)
        summary = summarize_phase(temp_output, targets)
        os.replace(temp_output, out_dir / "phased.vcf.gz")
        os.replace(Path(f"{temp_output}.tbi"), out_dir / "phased.vcf.gz.tbi")
        os.replace(temp_dir / "whatshap.stdout.log", out_dir / "whatshap.stdout.log")
        os.replace(temp_dir / "whatshap.stderr.log", out_dir / "whatshap.stderr.log")
    (out_dir / "phase-block-summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (out_dir / "phase.complete.json").write_text(json.dumps({"status": "complete", "phased_vcf": "phased.vcf.gz"}) + "\n")
    print(
        "phased_records={phased_records} phase_blocks={phase_block_count} targets_connected={connected}".format(
            connected=summary["targets"]["connected_to_a_block_with_another_candidate"], **summary
        )
    )


def add_prepare_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--bam", required=True, help="Finalized coordinate-sorted, indexed BAM")
    parser.add_argument("--vcf", required=True, help="Indexed bgzip VCF (.vcf.gz/.vcf.bgz)")
    parser.add_argument("--fasta", required=True, help="Exact reference FASTA with .fai")
    parser.add_argument("--region", action="append", required=True, help="CONTIG:START-END, repeatable")
    parser.add_argument("--target", action="append", default=[], help="Optional target CONTIG:POSITION, repeatable")
    parser.add_argument("--sample", help="Required only for a multi-sample VCF")
    parser.add_argument("--ignore-read-groups", action="store_true", help="Use only after confirming a single BAM sample is the selected VCF sample")
    parser.add_argument("--min-qual", type=float, help="Optional inclusive QUAL threshold for candidate heterozygotes")
    parser.add_argument("--min-dp", type=int, help="Optional inclusive FORMAT/DP threshold; missing DP is excluded")
    parser.add_argument("--min-gq", type=int, help="Optional inclusive FORMAT/GQ threshold; PL may satisfy --min-pl-delta instead")
    parser.add_argument("--min-pl-delta", type=int, help="Optional best-vs-second FORMAT/PL separation threshold")
    parser.add_argument("--min-ad", type=int, help="Optional inclusive ref and ALT FORMAT/AD threshold when AD is present")
    parser.add_argument("--dp-from-ad", action="store_true", help="When FORMAT/DP is absent, explicitly derive depth as sum of valid FORMAT/AD values")
    parser.add_argument("--out", required=True, help="New output directory")


def main(argv: list[str] | None = None) -> int:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="subcommand", required=True)
    prepare_parser = subparsers.add_parser("prepare", help="create the regional candidate VCF")
    add_prepare_arguments(prepare_parser)
    phase_parser = subparsers.add_parser("phase", help="run WhatsHap exactly once on prepared input")
    phase_parser.add_argument("--out", required=True, help="Existing prepared output directory")
    phase_parser.add_argument("--whatshap", default="whatshap", help="WhatsHap executable/command (default: whatshap)")
    run_parser = subparsers.add_parser("run", help="prepare then phase")
    add_prepare_arguments(run_parser)
    run_parser.add_argument("--whatshap", default="whatshap", help="WhatsHap executable/command (default: whatshap)")
    args = parser.parse_args(argv)
    try:
        if args.subcommand == "prepare":
            atomic_prepare(args)
        elif args.subcommand == "phase":
            phase(args)
        else:
            atomic_prepare(args)
            phase(args)
        return 0
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
