#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0"]
# ///
"""Bounded AlphaGenome scoring; requires a separately verified protected transport.

No credential value, phenotype, genotype or sample ID is read here. The trusted
client factory consumes only the gateway-injected opaque sentinel. Never run
with a plaintext API key. This script has no direct-network client fallback.
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import fcntl
import hashlib
import importlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import re
import shutil
import stat
import sys
from typing import Any

from prepare_alphagenome_api_stage import TARGET, WINDOW, _interval, _variants

SCORERS = ("RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS")
LIMITATIONS = [
    "Model predictions, not measured expression or RNA splicing.",
    "Scores do not establish pathogenicity, a diagnosis, or cis/trans phase.",
    "Single-variant predictions do not model combined haplotypes or aneuploidy.",
    "Raw and quantile scores are scorer-specific; not clinical probabilities.",
    "No signal does not exclude coding/protein effects or unmodelled tissues.",
]


class StageError(RuntimeError):
    """Static internal error codes only; no input values or remote messages."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path: Path, value: Any) -> None:
    # Exclusive create: existing evidence/checkpoints are never replaced.
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def read_json(path: Path) -> Any:
    return json.loads(path.read_bytes())


def check_private(path: Path, *, directory: bool = False) -> None:
    info = path.lstat()
    kind_ok = stat.S_ISDIR(info.st_mode) if directory else stat.S_ISREG(info.st_mode)
    if not kind_ok or info.st_mode & 0o077:
        raise StageError("UNPROTECTED_OR_UNEXPECTED_FILE_TYPE")


def validated_manifest(manifest_path: Path, source_path: Path) -> tuple[dict, dict]:
    for path in (manifest_path, source_path):
        check_private(path)
    manifest_hash = sha256(manifest_path)
    source_hash = sha256(source_path)
    manifest, source = read_json(manifest_path), read_json(source_path)
    if manifest.get("schema") != "alphagenome_api_stage_manifest_v1":
        raise StageError("MANIFEST_SCHEMA")
    if manifest.get("source", {}).get("sha256") != source_hash:
        raise StageError("SOURCE_HASH_MISMATCH")
    if manifest.get("target") != TARGET or manifest.get("reference_build") != "GRCh38":
        raise StageError("TARGET_OR_BUILD_MISMATCH")
    if source.get("status") != "complete" or source.get("reference_build", source.get("genome_build")) != "GRCh38":
        raise StageError("SOURCE_NOT_COMPLETE_GRCH38")
    if manifest.get("requested_scorers") != list(SCORERS):
        raise StageError("SCORER_SCOPE_MISMATCH")
    if manifest.get("requested_score_scope") != ["differential_splice", "differential_expression"]:
        raise StageError("SCORE_SCOPE_MISMATCH")
    if (manifest.get("scope") != "proposed AlphaGenome API differential scoring of two BUB1B SNVs"
            or source.get("scope") != "local exact-allele AlphaGenome Atlas annotation of two BUB1B SNVs"):
        raise StageError("SOURCE_SCOPE_MISMATCH")
    expected = []
    for variant in _variants(source):
        reference, query = _interval(variant["chrom"], variant["pos"])
        expected.append({**variant, "ref_interval": reference, "query_interval_1mb": query})
    if manifest.get("variants") != expected:
        raise StageError("VARIANTS_OR_COORDINATES_MISMATCH")
    binding = {"manifest_sha256": manifest_hash, "source_sha256": source_hash,
               "target": TARGET, "reference_build": "GRCh38", "scorers": list(SCORERS),
               "sequence_length": WINDOW, "variant_count": 2}
    return manifest, binding


@contextlib.contextmanager
def quiet_sdk():
    # SDK exceptions may contain request coordinates or transport metadata.
    # Suppress both Python and native writes. All errors are reported as codes.
    previous = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    sys.stdout.flush()
    sys.stderr.flush()
    saved_out, saved_err = os.dup(1), os.dup(2)
    try:
        with open(os.devnull, "w") as sink:
            os.dup2(sink.fileno(), 1)
            os.dup2(sink.fileno(), 2)
            with contextlib.redirect_stdout(sink), contextlib.redirect_stderr(sink):
                yield
    finally:
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        os.close(saved_out)
        os.close(saved_err)
        logging.disable(previous)


def make_client(spec: str, timeout: int):
    if not re.fullmatch(r"[A-Za-z_][\w.]*:[A-Za-z_]\w*", spec):
        raise StageError("INVALID_CLIENT_FACTORY")
    # Presence only; the selected factory validates native keys versus protected sentinels.
    if "ALPHAGENOME_API_KEY" not in os.environ:
        raise StageError("ALPHAGENOME_CREDENTIAL_NOT_PROVIDED")
    module, name = spec.split(":")
    return getattr(importlib.import_module(module), name)(target=TARGET, timeout_seconds=timeout)


def safe_failure(error: Exception) -> dict:
    result = {"status": "failed", "error": "STAGE_FAILURE", "at": utc()}
    if isinstance(error, StageError):
        result["error"] = str(error)
    # Whitelist only enum status names; never expose details/debug_error_string.
    try:
        code = error.code().name
        if code in {"UNAUTHENTICATED", "PERMISSION_DENIED", "UNAVAILABLE", "DEADLINE_EXCEEDED",
                    "RESOURCE_EXHAUSTED", "INVALID_ARGUMENT", "INTERNAL", "CANCELLED"}:
            result["grpc_status"] = code
    except Exception:
        pass
    return result


def sdk_version() -> str:
    return importlib.metadata.version("alphagenome")


def smoke(client, output: Path) -> dict:
    from alphagenome.models import dna_client
    metadata = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    counts = {}
    for name in SCORERS:
        frame = getattr(metadata, name.lower(), None)
        if frame is None or len(frame) == 0:
            raise StageError("METADATA_MODALITY_MISSING")
        counts[name] = len(frame)
    result = {"status": "metadata_smoke_complete", "at": utc(), "target": TARGET,
              "organism": "HOMO_SAPIENS", "sdk_version": sdk_version(),
              "modality_track_counts": counts, "patient_requests": 0}
    write_json(output / "metadata-smoke.json", result)
    return result


def score_one(client, variant: dict, output: Path) -> dict:
    import numpy as np
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    allele = genome.Variant(chromosome=variant["chrom"], position=variant["pos"],
                            reference_bases=variant["ref"], alternate_bases=variant["alt"])
    interval = allele.reference_interval.resize(WINDOW)
    planned = variant["query_interval_1mb"]
    if (interval.start, interval.end) != (planned["start"], planned["end"]):
        raise StageError("SDK_COORDINATE_MISMATCH")
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[name] for name in SCORERS]
    scores = client.score_variant(interval=interval, variant=allele, variant_scorers=scorers,
                                  organism=dna_client.Organism.HOMO_SAPIENS)
    if len(scores) != len(SCORERS):
        raise StageError("SCORER_RESULT_COUNT")
    stats = {}
    for name, scorer, matrix in zip(SCORERS, scorers, scores, strict=True):
        if repr(matrix.uns.get("variant_scorer")) != repr(scorer):
            raise StageError("SCORER_RESULT_ORDER")
        # Keep original numeric matrices and annotations independent of tidy filtering.
        arrays = {"raw_score": np.asarray(matrix.X)}
        # Newer AnnData exposes X again as layers[None]; X is already retained.
        arrays.update({key: np.asarray(value) for key, value in matrix.layers.items() if key is not None})
        if any(array.dtype.hasobject for array in arrays.values()):
            raise StageError("UNEXPECTED_SCORE_DTYPE")
        np.savez_compressed(output / f"{name}.npz", **arrays)
        matrix.obs.to_csv(output / f"{name}.obs.tsv.gz", sep="\t", index=True)
        matrix.var.to_csv(output / f"{name}.tracks.tsv.gz", sep="\t", index=True)
        frame = variant_scorers.tidy_anndata(matrix, match_gene_strand=True)
        if frame.empty:
            # No score is distinct from a score of zero; preserve the raw shape.
            import pandas as pd
            frame = pd.DataFrame(columns=["raw_score", "quantile_score", "output_type"])
        required = {"raw_score", "output_type"}
        if not required.issubset(frame.columns):
            raise StageError("TIDY_SCORE_SCHEMA")
        quantiles_available = "quantile_score" in frame.columns
        if not quantiles_available:
            frame["quantile_score"] = np.nan
        if len(frame) and not frame["output_type"].astype(str).str.contains(name, regex=False).all():
            raise StageError("OUTPUT_TYPE_MISMATCH")
        frame.to_csv(output / f"{name}.scores.tsv.gz", sep="\t", index=False)
        numeric = frame["raw_score"].to_numpy(dtype=float)
        quantile = frame["quantile_score"].to_numpy(dtype=float)
        finite = np.isfinite(numeric)
        stats[name] = {"tidy_rows": len(frame), "finite_raw_rows": int(finite.sum()),
                       "quantiles_available": quantiles_available,
                       "raw_matrix_shape": list(matrix.shape), "scorer": repr(scorer),
                       "max_absolute_raw_score": float(np.abs(numeric[finite]).max()) if finite.any() else None,
                       "quantile_abs_over_0995_rows": int((np.abs(quantile[np.isfinite(quantile)]) > .995).sum())}
        # Keep discovery hits per scorer, never rank raw values across modalities.
        ranked = frame.loc[finite].sort_values("raw_score", key=abs, ascending=False).head(20)
        ranked.to_csv(output / f"{name}.top20.tsv", sep="\t", index=False)
    for path in output.iterdir():
        path.chmod(0o600)
    return {"status": "complete", "at": utc(), "sdk_version": sdk_version(),
            "modalities": stats, "limitations": LIMITATIONS}


def verified_checkpoint(path: Path, binding: dict, ordinal: int) -> dict:
    check_private(path, directory=True)
    check_private(path / "complete.json")
    completion = read_json(path / "complete.json")
    if completion.get("binding") != binding or completion.get("ordinal") != ordinal:
        raise StageError("CHECKPOINT_BINDING_MISMATCH")
    if completion.get("status") != "complete":
        raise StageError("CHECKPOINT_NOT_COMPLETE")
    expected_files = {"request-started.json"}
    expected_files.update(f"{name}.{suffix}" for name in SCORERS
                          for suffix in ("npz", "obs.tsv.gz", "tracks.tsv.gz", "scores.tsv.gz", "top20.tsv"))
    if set(completion.get("files_sha256", {})) != expected_files:
        raise StageError("CHECKPOINT_FILE_SET_MISMATCH")
    for filename, expected_hash in completion["files_sha256"].items():
        if Path(filename).name != filename:
            raise StageError("CHECKPOINT_INVALID_FILENAME")
        check_private(path / filename)
        if sha256(path / filename) != expected_hash:
            raise StageError("CHECKPOINT_HASH_MISMATCH")
    return completion


def run(args, client_factory=make_client) -> dict:
    os.umask(0o077)
    if args.mode == "smoke":
        args.out.mkdir(mode=0o700)  # must be new, cannot overwrite prior smoke
        with quiet_sdk():
            return smoke(client_factory(args.client_factory, args.timeout_seconds), args.out)
    if not args.manifest or not args.atlas_summary:
        raise StageError("MANIFEST_AND_SOURCE_REQUIRED")
    manifest, binding = validated_manifest(args.manifest, args.atlas_summary)
    if args.mode == "validate":
        return {"status": "offline_validated", **binding, "patient_requests": 0}
    args.out.mkdir(mode=0o700, exist_ok=True)
    check_private(args.out, directory=True)
    with (args.out / "run.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise StageError("ANOTHER_RUN_ACTIVE") from None
        identity = args.out / "binding.json"
        if identity.exists():
            if read_json(identity) != binding:
                raise StageError("OUTPUT_BOUND_TO_OTHER_INPUT")
        else:
            write_json(identity, binding)
        if (args.out / "summary.json").exists():
            check_private(args.out / "summary.json")
            summary = read_json(args.out / "summary.json")
            if summary.get("binding") != binding or summary.get("status") != "complete":
                raise StageError("SUMMARY_BINDING_MISMATCH")
            for index in range(1, 3):
                verified_checkpoint(args.out / f"variant-{index:02d}", binding, index)
            return {"status": "already_complete", "completed_variants": 2, "patient_requests": 0}
        client = None
        completed, submitted = [], 0
        for index, variant in enumerate(manifest["variants"], 1):
            final = args.out / f"variant-{index:02d}"
            pending = args.out / f"variant-{index:02d}.pending"
            if final.exists():
                completed.append(verified_checkpoint(final, binding, index))
                continue
            if pending.exists():
                # An interrupted remote request may have executed. No silent resend.
                raise StageError("PARTIAL_ATTEMPT_REQUIRES_REVIEW")
            if shutil.disk_usage(args.out).free < 60 * 1024**3:
                raise StageError("DISK_FLOOR_60_GIB")
            _, current = validated_manifest(args.manifest, args.atlas_summary)
            if current != binding:
                raise StageError("SOURCE_CHANGED_BEFORE_REQUEST")
            pending.mkdir(mode=0o700)
            write_json(pending / "request-started.json", {"at": utc(), "ordinal": index, "binding": binding})
            try:
                with quiet_sdk():
                    if client is None:
                        client = client_factory(args.client_factory, args.timeout_seconds)
                    submitted += 1
                    result = score_one(client, variant, pending)
                _, current = validated_manifest(args.manifest, args.atlas_summary)
                if current != binding:
                    raise StageError("SOURCE_CHANGED_DURING_REQUEST")
                result.update(binding=binding, ordinal=index)
                result["files_sha256"] = {p.name: sha256(p) for p in pending.iterdir() if p.is_file()}
                write_json(pending / "complete.json", result)
                pending.rename(final)
                completed.append(result)
            except Exception as error:
                write_json(pending / "failure.json", safe_failure(error))
                raise
        summary = {"status": "complete", "at": utc(), "binding": binding,
                   "completed_variants": len(completed), "variant_results": completed,
                   "limitations": LIMITATIONS, "sources_unchanged": True}
        write_json(args.out / "summary.json", summary)
        return {"status": "complete", "completed_variants": 2, "patient_requests": submitted}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("validate", "smoke", "run"))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--atlas-summary", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--client-factory", default="")
    parser.add_argument("--timeout-seconds", type=int, default=600)
    args = parser.parse_args(argv)
    if args.mode in ("smoke", "run") and os.environ.get("MVA_ALLOW_NETWORK") != "1":
        parser.error("set MVA_ALLOW_NETWORK=1 for opt-in network execution")
    if not 1 <= args.timeout_seconds <= 1800:
        parser.error("timeout must be in 1..1800 seconds")
    try:
        result = run(args)
    except Exception as error:
        print(json.dumps(safe_failure(error), sort_keys=True))
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
