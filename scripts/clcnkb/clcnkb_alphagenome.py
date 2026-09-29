#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""Bounded, checkpointed AlphaGenome follow-up for CLCNKB INS79.

Prepare/analyse/verify run locally. Predict requires an explicitly configured,
credential-injecting runtime. The native PredictVariant insertion result is primary;
manually constructed insertion-containing sequences are sensitivity analyses.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

from alphagenome_mechanism import fetch, js, private, sha, tsv
from analyze_alphagenome_profiles import align_arrays, normalize_junction_coordinates
from run_alphagenome_api_stage import quiet_sdk

BASE = Path(os.environ.get("MVA_DATA_ROOT", ".")).resolve()
REPORT = BASE / "reports/clcnkb-functional-20260911"
RESULT = BASE / "results/clcnkb-functional-20260911"
OPS = BASE / "operations/clcnkb-functional-20260911"
SOURCE_REPORT = BASE / "reports/clcn-reconstruction-20260911/CLCNKA_CLCNKB_STRUCTURAL_RECONSTRUCTION.md"
SOURCE_DIR = BASE / "results/clcn-reconstruction-20260911/analysis"
SOURCE_FILES = [SOURCE_REPORT, SOURCE_DIR / "reconstructed-allele.vcf", SOURCE_DIR / "refined-allele.json", SOURCE_DIR / "transcript-context.json"]

CHROM = "chr1"
POS = 16050023
REF = "T"
# Public gnomAD v4 allele rs1553127751; see METHODS.md for provenance.
INSERT = "GCCATTATTTTTTCCTGCCCAGACAATGCCCATGCAGTGATCTGGGCCCCCAAGGACCCAGCTTCACCCCCACAGCACC"
ALT = REF + INSERT
SNP_POS = 16049963
SNP_REF = "T"
SNP_ALT = "C"
LENGTH = 2**20
MODES = ["RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"]
TERMS = ["CL:0002306", "CL:0002518", "CL:0002584", "CL:1000507", "UBERON:0001225", "UBERON:0001293", "UBERON:0002113"]
GENE = "CLCNKB"
GENE_START, GENE_END = 16043782, 16057326
FOCUS_START, FOCUS_END = 16049250, 16051200
EXONS = [(16049838, 16049916, 10), (16050516, 16050645, 11)]  # 1-based inclusive


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def interval_dict():
    # Match SDK Variant.reference_interval.resize for a one-base VCF anchor.
    start = POS - 1 - (LENGTH - 1) // 2
    return {"chrom": CHROM, "start": start, "end": start + LENGTH}


def build_sequences(reference: str, start0: int):
    """Return fixed-length reference, insertion-only and combined sequences.

    Insertion-bearing strings are right-trimmed by 79 bp. This exact matched
    construction lets combined-vs-insertion isolate the SNP in sensitivity
    analysis, but does not replace native PredictVariant insertion inference.
    """
    assert len(reference) == LENGTH
    anchor = POS - 1 - start0
    snp = SNP_POS - 1 - start0
    assert reference[anchor] == REF and reference[snp] == SNP_REF
    ins = (reference[: anchor + 1] + INSERT + reference[anchor + 1 :])[:LENGTH]
    combined_ref = reference[:snp] + SNP_ALT + reference[snp + 1 :]
    combined = (combined_ref[: anchor + 1] + INSERT + combined_ref[anchor + 1 :])[:LENGTH]
    assert len(ins) == len(combined) == LENGTH
    assert sum(a != b for a, b in zip(ins, combined)) == 1
    assert ins[snp] == SNP_REF and combined[snp] == SNP_ALT
    assert ins[anchor + 1 : anchor + 1 + len(INSERT)] == INSERT
    return {"reference": reference, "insertion_manual": ins, "combined_sensitivity": combined}


def prepare():
    os.umask(0o077)
    for p in (REPORT, RESULT, OPS):
        private(p)
    source_hashes = {str(p): sha(p) for p in SOURCE_FILES}
    allele = json.loads((SOURCE_DIR / "refined-allele.json").read_text())
    tx = json.loads((SOURCE_DIR / "transcript-context.json").read_text())
    assert allele["chrom"] == CHROM and allele["position1"] == POS
    assert allele["anchor_REF"] == REF and allele["inserted_sequence"] == INSERT and allele["insertion_length"] == 79
    assert allele["upstream_substitution_position1"] == SNP_POS and allele["upstream_substitution"] == "T>C"
    assert tx["gene"] == GENE and tx["intron_after_exon"] == 10
    assert tx["distance_from_upstream_exon_end"] == 107 and tx["distance_to_downstream_exon_start"] == 493
    vcf_rows = [x for x in (SOURCE_DIR / "reconstructed-allele.vcf").read_text().splitlines() if not x.startswith("#")]
    assert len(vcf_rows) == 1
    fields = vcf_rows[0].split("\t")
    assert fields[:5] == [CHROM, str(POS), "CLCNKB_INS79", REF, ALT]
    iv = interval_dict()
    reference = fetch(CHROM, iv["start"], iv["end"])
    sequences = build_sequences(reference, iv["start"])
    seq_hashes = {}
    for name, sequence in sequences.items():
        p = OPS / f"{name}.txt"
        p.write_text(sequence)
        p.chmod(0o600)
        seq_hashes[name] = sha(p)
    manifest = {
        "schema": "clcnkb_alphagenome_v1",
        "created_at": utc(),
        "assembly": "GRCh38",
        "sdk": "alphagenome==0.9.0",
        "model_context_bp": LENGTH,
        "interval_0based_halfopen": iv,
        "primary_variant": {"chrom": CHROM, "pos_1based": POS, "ref": REF, "alt": ALT, "inserted_bp": len(INSERT)},
        "nearby_variant": {"chrom": CHROM, "pos_1based": SNP_POS, "ref": SNP_REF, "alt": SNP_ALT, "distance_upstream_bp": 60},
        "requested_outputs": MODES,
        "ontology_terms": TERMS,
        "source_hashes": source_hashes,
        "sequence_hashes": seq_hashes,
        "api_tasks": ["score_INS79", "native_INS79", "native_SNP", "sequence_reference", "sequence_insertion_manual", "sequence_combined_sensitivity"],
        "interpretation_guardrails": [
            "Native PredictVariant INS79 is the primary indel-aware comparison.",
            "Manual insertion and combined sequences are sensitivity analyses only; their right-tail trim and alternate coordinate system can create positional artifacts.",
            "The nearby SNP and insertion are not treated as phased or as a diploid genotype.",
        ],
        "privacy": "Only selected alleles and public-reference-derived sequence windows; no reads, patient identifiers, phenotypes, BAM, or full VCF sent.",
    }
    js(OPS / "manifest.json", manifest)
    print(json.dumps({"status": "prepared", "api_tasks": len(manifest["api_tasks"]), "interval": iv}))


def serialize_output(output, path):
    from alphagenome.data import junction_data

    private(path)
    stats = {}
    for mode in MODES:
        data = getattr(output, mode.lower())
        if data is None:
            js(path / f"{mode}.metadata.json", {"available": False, "reason": "not_returned"})
            stats[mode] = {"available": False, "reason": "not_returned"}
            continue
        values = np.asarray(data.values)
        assert values.ndim == 2 and np.isfinite(values).all()
        np.savez_compressed(path / f"{mode}.npz", values=values)
        (path / f"{mode}.npz").chmod(0o600)
        tsv(path / f"{mode}.tracks.tsv", data.metadata)
        iv = None if data.interval is None else {"chrom": data.interval.chromosome, "start": data.interval.start, "end": data.interval.end, "strand": data.interval.strand}
        meta = {"available": True, "shape": list(values.shape), "interval": iv, "uns": str(data.uns)}
        if isinstance(data, junction_data.JunctionData):
            frame = pd.DataFrame([{"chrom": j.chromosome, "start": j.start, "end": j.end, "strand": j.strand} for j in data.junctions], columns=["chrom", "start", "end", "strand"])
            tsv(path / f"{mode}.junctions.tsv", frame)
            meta["kind"] = "junction"
        else:
            meta.update(kind="track", resolution=int(data.resolution))
        js(path / f"{mode}.metadata.json", meta)
        stats[mode] = meta
    return stats


def checkpoint(path: Path, task: str, call):
    path = private(path)
    done, started = path / "complete.json", path / "started.json"
    manifest_sha = sha(OPS / "manifest.json")
    if done.exists():
        state = json.loads(done.read_text())
        assert state["manifest_sha256"] == manifest_sha
        for rel, digest in state["hashes"].items():
            assert sha(path / rel) == digest
        return False
    if started.exists():
        raise RuntimeError(f"INCOMPLETE_REQUEST_REQUIRES_REVIEW:{task}")
    js(started, {"at": utc(), "task": task, "manifest_sha256": manifest_sha})
    began = time.monotonic()
    result = call()
    stats = result(path)
    hashes = {str(p.relative_to(path)): sha(p) for p in path.rglob("*") if p.is_file() and p.name != "complete.json"}
    js(done, {"status": "complete", "task": task, "at": utc(), "seconds": round(time.monotonic() - began, 2), "manifest_sha256": manifest_sha, "stats": stats, "hashes": hashes})
    print(json.dumps({"completed": task, "seconds": round(time.monotonic() - began, 1)}), flush=True)
    return True


def make_client():
    from alphagenome_native_client import create_client
    return create_client(timeout_seconds=900)

def predict():
    if os.environ.get("MVA_ALLOW_NETWORK") != "1":
        raise RuntimeError("set MVA_ALLOW_NETWORK=1 for opt-in prediction")
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    os.umask(0o077)
    manifest = json.loads((OPS / "manifest.json").read_text())
    assert manifest["primary_variant"]["alt"] == ALT and manifest["nearby_variant"]["pos_1based"] == SNP_POS
    for name, digest in manifest["sequence_hashes"].items():
        assert sha(OPS / f"{name}.txt") == digest
    ivd = manifest["interval_0based_halfopen"]
    interval = genome.Interval(CHROM, ivd["start"], ivd["end"])
    outputs = [dna_client.OutputType[x] for x in MODES]
    with quiet_sdk():
        client = make_client()
    completed = 0

    score_dir = RESULT / "scores/INS79"
    def score_call():
        insertion = genome.Variant(CHROM, POS, REF, ALT)
        scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[x] for x in MODES]
        with quiet_sdk():
            response = client.score_variant(interval=interval, variant=insertion, variant_scorers=scorers, organism=dna_client.Organism.HOMO_SAPIENS)
        def save(path):
            assert len(response) == len(scorers)
            summary = []
            for mode, scorer, adata in zip(MODES, scorers, response, strict=True):
                assert repr(adata.uns["variant_scorer"]) == repr(scorer)
                arrays = {"raw": np.asarray(adata.X)}
                arrays.update({str(k): np.asarray(v) for k, v in adata.layers.items()})
                assert all(np.isfinite(a).all() for a in arrays.values())
                np.savez_compressed(path / f"{mode}.npz", **arrays)
                (path / f"{mode}.npz").chmod(0o600)
                tsv(path / f"{mode}.obs.tsv", adata.obs.reset_index())
                tsv(path / f"{mode}.tracks.tsv", adata.var.reset_index())
                tidy = variant_scorers.tidy_anndata(adata, match_gene_strand=True)
                if tidy.empty:
                    tidy = pd.DataFrame(columns=["raw_score", "quantile_score", "gene_name", "output_type"])
                assert np.isfinite(tidy.raw_score.to_numpy(dtype=float)).all()
                tsv(path / f"{mode}.scores.tsv.gz", tidy)
                top = tidy.loc[tidy.raw_score.abs().sort_values(ascending=False).index].head(20)
                tsv(path / f"{mode}.top20.tsv", top)
                target = tidy[tidy.gene_name == GENE] if "gene_name" in tidy else tidy.iloc[0:0]
                summary.append({"mode": mode, "all_rows": len(tidy), "target_rows": len(target), "max_abs_all": float(tidy.raw_score.abs().max()) if len(tidy) else None, "max_abs_target": float(target.raw_score.abs().max()) if len(target) else None})
            js(path / "summary.json", summary)
            return {"modalities": summary}
        return save
    completed += checkpoint(score_dir, "score_INS79", score_call)

    def native_task(name, variant):
        def call():
            with quiet_sdk():
                prediction = client.predict_variant(interval, variant, requested_outputs=outputs, ontology_terms=TERMS)
            def save(path):
                return {"ref": serialize_output(prediction.reference, path / "ref"), "alt": serialize_output(prediction.alternate, path / "alt")}
            return save
        return call

    completed += checkpoint(RESULT / "profiles/native_INS79", "native_INS79", native_task("native_INS79", genome.Variant(CHROM, POS, REF, ALT)))
    completed += checkpoint(RESULT / "profiles/native_SNP", "native_SNP", native_task("native_SNP", genome.Variant(CHROM, SNP_POS, SNP_REF, SNP_ALT)))

    for sequence_name in ["reference", "insertion_manual", "combined_sensitivity"]:
        task = f"sequence_{sequence_name}"
        sequence = (OPS / f"{sequence_name}.txt").read_text()
        def sequence_call(sequence=sequence):
            with quiet_sdk():
                prediction = client.predict_sequence(sequence, interval=interval, requested_outputs=outputs, ontology_terms=TERMS)
            return lambda path: {"sequence": serialize_output(prediction, path / "sequence")}
        completed += checkpoint(RESULT / f"profiles/{task}", task, sequence_call)
    js(RESULT / "api.complete.json", {"status": "complete", "at": utc(), "completed_this_run": completed, "total_completed_requests": 6})
    print(json.dumps({"status": "complete", "completed_this_run": completed, "total_completed_requests": 6}))


def load_profile(path, mode):
    meta = json.loads((path / f"{mode}.metadata.json").read_text())
    if not meta.get("available"):
        return None, pd.DataFrame(), meta
    values = np.load(path / f"{mode}.npz", allow_pickle=False)["values"].astype(float)
    tracks = pd.read_csv(path / f"{mode}.tracks.tsv", sep="\t")
    return values, tracks, meta


def read_junctions(path, meta):
    try:
        frame = pd.read_csv(path / "SPLICE_JUNCTIONS.junctions.tsv", sep="\t")
    except pd.errors.EmptyDataError:
        frame = pd.DataFrame(columns=["chrom", "start", "end", "strand"])
    return normalize_junction_coordinates(frame, meta)


def compare_profiles(label, ref_path, alt_path):
    rows = []
    junctions = pd.DataFrame()
    for mode in MODES:
        ref, rtracks, rmeta = load_profile(ref_path, mode)
        alt, atracks, ameta = load_profile(alt_path, mode)
        if ref is None or alt is None:
            rows.append({"scenario": label, "mode": mode, "availability": "not_returned", "track_name": None, "biosample": None, "ontology_curie": None, "ref_sum": None, "alt_sum": None, "sum_ratio": None, "max_abs_delta": None})
            continue
        assert rtracks.name.tolist() == atracks.name.tolist()
        if mode == "SPLICE_JUNCTIONS":
            rj, aj = read_junctions(ref_path, rmeta), read_junctions(alt_path, ameta)
            junctions = align_arrays(rj, aj, ref, alt)
            for col in ["name", "biosample_name", "ontology_curie"]:
                junctions[col] = junctions.track.map(rtracks[col]) if col in rtracks else ""
            junctions["scenario"] = label
            junctions = junctions[(junctions.strand == "+") & (junctions.start < GENE_END) & (junctions.end > GENE_START)].copy()
            continue
        assert ref.shape == alt.shape and rmeta["resolution"] == ameta["resolution"]
        start, res = rmeta["interval"]["start"], rmeta["resolution"]
        coords = start + np.arange(ref.shape[0]) * res
        mask = (coords >= GENE_START) & (coords <= GENE_END)
        for idx, track in rtracks.iterrows():
            strand = str(track.get("strand", "."))
            if strand not in ["+", "."]:
                continue
            rv, av = ref[mask, idx], alt[mask, idx]
            denominator = float(rv.sum())
            rows.append({"scenario": label, "mode": mode, "availability": "returned", "track_name": track["name"], "biosample": str(track.get("biosample_name", "tissue agnostic")), "ontology_curie": str(track.get("ontology_curie", "")), "ref_sum": float(rv.sum()), "alt_sum": float(av.sum()), "sum_ratio": float(av.sum() / denominator) if denominator > 0 else None, "max_abs_delta": float(np.abs(av - rv).max()) if len(rv) else None})
    return pd.DataFrame(rows), junctions


def profile_comparisons():
    """Return only coordinate-compatible quantitative comparisons.

    PredictSequence receives fixed-length manual strings.  The INS79 string is
    shifted after the insertion and right-trimmed, so reference-to-INS79 bins
    are not homologous positions and must not be interpreted as an insertion
    effect.  The paired native API output is retained as the primary model
    comparison; the two manual insertion-bearing strings share the same
    constructed coordinate system, so their one-SNP difference is a bounded
    conditional sensitivity analysis.
    """
    return [
        ("INS79_native_primary", RESULT / "profiles/native_INS79/ref", RESULT / "profiles/native_INS79/alt"),
        ("SNP_native_separate", RESULT / "profiles/native_SNP/ref", RESULT / "profiles/native_SNP/alt"),
        ("manual_insertion_to_combined_conditional_SNP", RESULT / "profiles/sequence_insertion_manual/sequence", RESULT / "profiles/sequence_combined_sensitivity/sequence"),
    ]


def plot_native_ins(ref_path, alt_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = private(REPORT / "figures")
    made = []
    for term in TERMS:
        fig, axes = plt.subplots(4, 1, figsize=(13, 9), sharex=True, gridspec_kw={"height_ratios": [2, 1, 1, 1.6]})
        title = term
        for ax, mode in zip(axes[:3], MODES[:3]):
            ref, tracks, meta = load_profile(ref_path, mode)
            alt, _, _ = load_profile(alt_path, mode)
            if ref is None:
                ax.text(.03, .5, "Output not returned (unavailable, not zero)", transform=ax.transAxes)
                continue
            candidates = tracks.index[(tracks.get("ontology_curie", pd.Series("", index=tracks.index)) == term) & tracks.get("strand", pd.Series(".", index=tracks.index)).isin(["+", "."])].tolist()
            if mode == "SPLICE_SITES" and not candidates:
                candidates = tracks.index[tracks.get("strand", pd.Series(".", index=tracks.index)).isin(["+", "."])].tolist()
            if not candidates:
                ax.text(.03, .5, "No compatible returned track (unavailable, not zero)", transform=ax.transAxes)
                continue
            idx = candidates[0]
            coords = meta["interval"]["start"] + np.arange(ref.shape[0]) * meta["resolution"]
            mask = (coords >= FOCUS_START) & (coords <= FOCUS_END)
            ax.plot(coords[mask], ref[mask, idx], color="#2166ac", lw=.9, label="REF")
            ax.plot(coords[mask], alt[mask, idx], color="#b2182b", lw=.9, alpha=.8, label="INS79")
            title = str(tracks.loc[idx].get("biosample_name", term))
            ax.set_ylabel(mode, fontsize=8)
            ax.legend(fontsize=7, loc="upper right")
        ref, tracks, meta = load_profile(ref_path, "SPLICE_JUNCTIONS")
        alt, _, ameta = load_profile(alt_path, "SPLICE_JUNCTIONS")
        ax = axes[3]
        if ref is None:
            ax.text(.03, .5, "Junction output not returned (unavailable, not zero)", transform=ax.transAxes)
        else:
            jf = align_arrays(read_junctions(ref_path, meta), read_junctions(alt_path, ameta), ref, alt)
            jf["ontology_curie"] = jf.track.map(tracks["ontology_curie"]) if "ontology_curie" in tracks else ""
            jf = jf[(jf.ontology_curie == term) & (jf.start < FOCUS_END) & (jf.end > FOCUS_START)]
            both = jf[jf.returned_in_ref & jf.returned_in_alt].copy()
            if both.empty:
                ax.text(.03, .5, "No matched returned junctions in view (missing, not zero)", transform=ax.transAxes)
            else:
                both["peak"] = both[["ref", "alt"]].max(axis=1)
                both = both.nlargest(12, "peak")
                scale = max(float(both.peak.max()), 1.0)
                for _, row in both.iterrows():
                    x = np.linspace(row.start, row.end, 100)
                    phase = np.sin(np.linspace(0, np.pi, 100))
                    ax.plot(x, phase * (.15 + .75 * row.ref / scale), color="#2166ac", lw=.8)
                    ax.plot(x, -phase * (.15 + .75 * row.alt / scale), color="#b2182b", lw=.8)
            unmatched = jf[~(jf.returned_in_ref & jf.returned_in_alt)]
            if len(unmatched):
                ax.text(.01, .92, f"{len(unmatched)} coordinate/track rows returned on one side only; not imputed as zero", transform=ax.transAxes, fontsize=7)
        for start1, end1, exon in EXONS:
            ax.plot([start1 - 1, end1], [-.96, -.96], color="black", lw=5)
            ax.text((start1 + end1) / 2, -.83, f"E{exon}", ha="center", fontsize=7)
        ax.axhline(0, color="grey", lw=.5)
        ax.set_ylim(-1.1, 1.1)
        ax.set_ylabel("junctions\nREF↑ INS↓", fontsize=8)
        for panel in axes:
            panel.axvline(POS - 1, color="black", ls=":", lw=.9)
            panel.axvline(SNP_POS - 1, color="grey", ls="--", lw=.7)
            panel.set_xlim(FOCUS_START, FOCUS_END)
        axes[0].set_title("UNALIGNED TRACK OVERLAYS: REF genomic positions / ALT sequence positions (+79 downstream)", fontsize=8)
        axes[-1].set_title("Junction panel: matched server-returned coordinate labels; reference exon markers", fontsize=8)
        axes[-1].set_xlabel("Offset display positions — track panels are NOT indel-aligned GRCh38 comparisons")
        fig.suptitle(f"CLCNKB INS79 native PredictVariant — {title}\n1 Mb context; diagnostic overlays, not evidence of shifted splice sites")
        fig.tight_layout(rect=[0, 0, 1, .94])
        path = figures / f"INS79_{term.replace(':', '_')}.png"
        fig.savefig(path, dpi=150)
        path.chmod(0o600)
        plt.close(fig)
        made.append(str(path))
    return made


def analyze():
    os.umask(0o077)
    score = RESULT / "scores/INS79"
    all_scores, top_hits, target_scores = [], [], []
    score_summary = []
    for mode in MODES:
        frame = pd.read_csv(score / f"{mode}.scores.tsv.gz", sep="\t")
        frame["mode"] = mode
        all_scores.append(frame)
        ranked = frame.loc[frame.raw_score.abs().sort_values(ascending=False).index].head(20)
        top_hits.append(ranked)
        target = frame[frame.gene_name == GENE] if "gene_name" in frame else frame.iloc[0:0]
        target_scores.append(target)
        top_all = ranked.iloc[0] if len(ranked) else None
        top_target = target.loc[target.raw_score.abs().idxmax()] if len(target) else None
        score_summary.append({
            "mode": mode,
            "all_rows": len(frame),
            "target_rows": len(target),
            "target_status": "returned" if len(target) else "no_CLCNKB_rows_returned_not_zero",
            "max_abs_all": float(frame.raw_score.abs().max()) if len(frame) else None,
            "top_all_gene": None if top_all is None else str(top_all.get("gene_name", "")),
            "top_all_raw": None if top_all is None else float(top_all.raw_score),
            "max_abs_CLCNKB": None if top_target is None else float(abs(top_target.raw_score)),
            "top_CLCNKB_raw": None if top_target is None else float(top_target.raw_score),
            "top_CLCNKB_quantile": None if top_target is None else float(top_target.quantile_score),
            "top_CLCNKB_biosample": None if top_target is None else str(top_target.get("biosample_name", "tissue agnostic")),
        })
    tsv(REPORT / "discovery_top_hits.tsv", pd.concat(top_hits, ignore_index=True))
    tsv(REPORT / "CLCNKB_scores.tsv", pd.concat(target_scores, ignore_index=True) if any(len(x) for x in target_scores) else pd.DataFrame(columns=["mode", "raw_score", "quantile_score", "gene_name"]))
    js(REPORT / "score_summary.json", score_summary)

    comparisons = profile_comparisons()
    js(REPORT / "manual_reference_to_insertion_not_compared.json", {
        "status": "not_compared",
        "reason": "The manual INS79 sequence shifts downstream bases and right-trims 79 bases; its array positions are not coordinate-mapped to the reference sequence.",
        "valid_manual_comparison": "manual_insertion_to_combined_conditional_SNP",
    })
    tracks, junctions = [], []
    profile_summary = []
    for label, ref, alt in comparisons:
        rows, jf = compare_profiles(label, ref, alt)
        tracks.append(rows)
        junctions.append(jf)
        matched = jf[jf.returned_in_ref & jf.returned_in_alt] if len(jf) else jf
        profile_summary.append({
            "scenario": label,
            "track_rows": len(rows),
            "unavailable_track_rows": int((rows.availability != "returned").sum()) if len(rows) else 0,
            "junction_rows": len(jf),
            "matched_junction_rows": len(matched),
            "unmatched_junction_rows": int((~(jf.returned_in_ref & jf.returned_in_alt)).sum()) if len(jf) else 0,
            "max_abs_matched_junction_delta": float(matched.delta.abs().max()) if len(matched) else None,
        })
    track_frame = pd.concat(tracks, ignore_index=True)
    junction_frame = pd.concat(junctions, ignore_index=True) if any(len(x) for x in junctions) else pd.DataFrame()
    tsv(REPORT / "profile_track_changes.tsv", track_frame)
    tsv(REPORT / "profile_junction_changes.tsv", junction_frame)
    js(REPORT / "profile_summary.json", profile_summary)
    figures = plot_native_ins(RESULT / "profiles/native_INS79/ref", RESULT / "profiles/native_INS79/alt")

    # Cross-method baseline is an artifact diagnostic, never an effect estimate.
    baseline = []
    for mode in MODES:
        native, _, _ = load_profile(RESULT / "profiles/native_INS79/ref", mode)
        sequence, _, _ = load_profile(RESULT / "profiles/sequence_reference/sequence", mode)
        baseline.append({"mode": mode, "available_both": native is not None and sequence is not None, "shapes_equal": native is not None and sequence is not None and native.shape == sequence.shape, "max_abs_difference": float(np.max(np.abs(native - sequence))) if native is not None and sequence is not None and native.shape == sequence.shape else None})
    js(REPORT / "cross_method_reference_baseline.json", baseline)
    js(REPORT / "analysis_summary.json", {"score_modalities": score_summary, "profile_comparisons": profile_summary, "figures": figures, "missingness_rule": "not returned/unmatched is unavailable and is never replaced by zero", "primary_inference": "INS79 ScoreVariant indel-aware scalar scores; native paired profiles require coordinate caveats", "manual_reference_to_insertion": "not quantitatively compared because coordinates are shifted/right-trimmed", "combined_role": "sensitivity analysis only; phase and diploid state unresolved"})
    annotate_coordinate_limits()
    print(json.dumps({"status": "analysed", "figures": len(figures), "score_rows": sum(len(x) for x in all_scores), "track_comparisons": len(track_frame), "junction_rows": len(junction_frame)}))


def annotate_coordinate_limits():
    """Label existing derived metrics without changing predictions or numbers."""
    limits = {
        "INS79_native_primary": {
            "track_coordinate_system": "unaligned_REF_genomic_ALT_sequence_positions",
            "track_metric_validity": "diagnostic_only_not_indel_aligned_effect",
            "junction_coordinate_system": "matched_server_returned_0based_labels_no_global_liftover_claim",
        },
        "SNP_native_separate": {
            "track_coordinate_system": "same_length_reference_positions",
            "track_metric_validity": "within_scenario_same_position_model_comparison",
            "junction_coordinate_system": "matched_server_returned_0based_labels",
        },
        "manual_insertion_to_combined_conditional_SNP": {
            "track_coordinate_system": "constructed_insertion_sequence_positions_not_GRCh38",
            "track_metric_validity": "conditional_SNP_only_same_constructed_positions",
            "junction_coordinate_system": "constructed_sequence_labels_not_GRCh38_downstream_of_INS79",
        },
    }
    for filename, kind in [("profile_track_changes.tsv", "track"), ("profile_junction_changes.tsv", "junction")]:
        frame = pd.read_csv(REPORT / filename, sep="\t")
        frame["coordinate_system"] = frame.scenario.map(lambda scenario: limits[scenario][f"{kind}_coordinate_system"])
        frame["metric_validity"] = frame.scenario.map(lambda scenario: limits[scenario]["track_metric_validity"] if kind == "track" else "matched_labels_only_not_a_validated_isoform_or_genotype")
        tsv(REPORT / filename, frame)
    summary = json.loads((REPORT / "profile_summary.json").read_text())
    for row in summary:
        row.update(limits[row["scenario"]])
    js(REPORT / "profile_summary.json", summary)
    analysis = json.loads((REPORT / "analysis_summary.json").read_text())
    analysis.update(profile_comparisons=summary, primary_inference="INS79 ScoreVariant indel-aware scalar scores; raw native insertion tracks are unaligned diagnostic overlays", coordinate_limits="coordinate_interpretation.json")
    js(REPORT / "analysis_summary.json", analysis)
    js(REPORT / "coordinate_interpretation.json", {
        "scenarios": limits,
        "native_INS79_track_metrics": "Existing fixed-bin sums, ratios and max deltas retained for provenance only; they compare nonhomologous positions downstream and must not be used as biological effect estimates.",
        "manual_reference_to_insertion": "Not compared; downstream shift and right-tail trim are not mapped.",
        "manual_conditional_SNP": "Paired constructed strings differ only at the SNP, but downstream coordinate labels are not exact GRCh38 locations.",
        "native_junctions": "Internally matched server-returned labels; no documented global indel liftover or validated transcript mechanism claimed.",
        "score_variant_documentation": "https://www.alphagenomedocs.com/faqs.html",
    })


def finalize_labels():
    os.umask(0o077)
    annotate_coordinate_limits()
    figures = plot_native_ins(RESULT / "profiles/native_INS79/ref", RESULT / "profiles/native_INS79/alt")
    print(json.dumps({"status": "coordinate_labels_corrected", "figures": len(figures), "api_requests": 0}))


def verify():
    os.umask(0o077)
    manifest = json.loads((OPS / "manifest.json").read_text())
    source_checks = {path: {"expected": digest, "observed": sha(path)} for path, digest in manifest["source_hashes"].items()}
    assert all(x["expected"] == x["observed"] for x in source_checks.values())
    for name, digest in manifest["sequence_hashes"].items():
        assert sha(OPS / f"{name}.txt") == digest
    checkpoints = []
    for complete in sorted(RESULT.glob("**/complete.json")):
        state = json.loads(complete.read_text())
        assert state["status"] == "complete"
        for rel, digest in state["hashes"].items():
            assert sha(complete.parent / rel) == digest
        checkpoints.append({"path": str(complete), "files": len(state["hashes"]), "sha256": sha(complete)})
    assert len(checkpoints) == 6
    arrays = []
    for path in sorted(RESULT.rglob("*.npz")):
        with np.load(path, allow_pickle=False) as data:
            shapes = {name: list(data[name].shape) for name in data.files}
            assert all(np.isfinite(data[name]).all() for name in data.files)
            tracks = pd.read_csv(path.with_name(path.stem + ".tracks.tsv"), sep="\t")
            assert all(data[name].ndim == 2 and data[name].shape[1] == len(tracks) for name in data.files)
            metadata = path.with_name(path.stem + ".metadata.json")
            if metadata.exists():
                meta = json.loads(metadata.read_text())
                assert shapes["values"] == meta["shape"]
                if meta["kind"] == "junction":
                    assert shapes["values"][0] == len(pd.read_csv(path.with_name(path.stem + ".junctions.tsv"), sep="\t"))
            else:
                assert all(data[name].shape[0] == len(pd.read_csv(path.with_name(path.stem + ".obs.tsv"), sep="\t")) for name in data.files)
        arrays.append({"path": str(path), "shapes": shapes})
    assert arrays
    for root in (REPORT, RESULT, OPS):
        root.chmod(0o700)
        for path in root.rglob("*"):
            path.chmod(0o700 if path.is_dir() else 0o600)
    permission_checks = {"dirs_0700": all((p.stat().st_mode & 0o777) == 0o700 for root in (REPORT, RESULT, OPS) for p in [root, *[x for x in root.rglob("*") if x.is_dir()]]), "files_0600": all((p.stat().st_mode & 0o777) == 0o600 for root in (REPORT, RESULT, OPS) for p in root.rglob("*") if p.is_file())}
    assert all(permission_checks.values())
    excluded = {REPORT / "verification.json", OPS / "completion.json"}
    output_hashes = {str(path): sha(path) for root in (REPORT, RESULT, OPS) for path in root.rglob("*") if path.is_file() and path not in excluded}
    result = {"status": "pass", "at": utc(), "source_hashes_preserved": source_checks, "completed_api_requests": len(checkpoints), "checkpoints": checkpoints, "numeric_arrays": len(arrays), "numeric_arrays_finite": True, "permissions": permission_checks, "disk_free_gib": round(shutil.disk_usage(BASE).free / 1024**3, 2), "output_hashes": output_hashes, "limits": ["No genotype or phase inferred", "Combined sequence is sensitivity-only", "Manual insertion coordinate shift/right-tail trim can cause positional artifacts", "Missing/unreturned outputs are not zero", "Quantile ranks are not effect sizes"]}
    result["numeric_array_dimensions_match_metadata"] = True
    result["finalizer_code"] = {"path": str(Path(__file__).resolve()), "sha256": sha(__file__)}
    result["limits"].append("Native INS79 track sums and equal-bin deltas are unaligned diagnostics, not insertion-effect estimates")
    js(REPORT / "verification.json", result)
    js(OPS / "completion.json", {"status": "complete", "at": utc(), "api_requests": 6, "report": str(REPORT / "report.md"), "report_sha256": sha(REPORT / "report.md"), "verification": str(REPORT / "verification.json"), "verification_sha256": sha(REPORT / "verification.json"), "manifest_exclusions": [str(p) for p in sorted(excluded)], "limitations": result["limits"]})
    print(json.dumps({"status": "pass", "completed_api_requests": len(checkpoints), "numeric_arrays": len(arrays), "permissions": permission_checks}))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["prepare", "predict", "analyze", "finalize-labels", "verify"], required=True)
    args = parser.parse_args()
    try:
        {"prepare": prepare, "predict": predict, "analyze": analyze, "finalize-labels": finalize_labels, "verify": verify}[args.mode]()
    except Exception as error:
        if OPS.exists():
            js(OPS / f"{args.mode}-error.json", {"at": utc(), "mode": args.mode, "error_class": type(error).__name__, "frames": [{"file": Path(f.filename).name, "line": f.lineno} for f in traceback.extract_tb(error.__traceback__)]})
        print(json.dumps({"status": "failed", "mode": args.mode, "error_class": type(error).__name__}))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
