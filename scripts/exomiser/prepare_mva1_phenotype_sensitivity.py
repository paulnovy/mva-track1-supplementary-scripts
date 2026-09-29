#!/usr/bin/env python3
"""Build source-linked phenotype profiles from explicit private curation inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}

def load_config(path: Path) -> dict:
    """Load private curation and ordered profiles; no case metadata is built in."""
    config = json.loads(path.read_text(encoding="utf-8"))
    curation = config["curation"]
    profiles = config["profiles"]
    if not isinstance(curation, list) or not curation or not isinstance(profiles, dict) or not profiles:
        raise ValueError("nonempty_curation_and_profiles_required")
    by_id = {}
    for row in curation:
        for key in ("key", "hpo_id", "label", "feature_paragraphs", "onset", "subject", "justification"):
            if key not in row:
                raise ValueError(f"curation_field_missing:{key}")
        if row["hpo_id"] in by_id:
            raise ValueError("duplicate_curation_hpo_id")
        if not row["feature_paragraphs"] or any(
            type(number) is not int or number < 1 for number in row["feature_paragraphs"]
        ):
            raise ValueError("positive_source_paragraph_numbers_required")
        by_id[row["hpo_id"]] = row
    for profile, term_ids in profiles.items():
        if not re.fullmatch(r"[A-Za-z0-9_-]+", profile):
            raise ValueError("invalid_profile_directory_name")
        if not isinstance(term_ids, list) or not term_ids or len(set(term_ids)) != len(term_ids):
            raise ValueError("nonempty_unique_profile_terms_required")
        for term_id in term_ids:
            if term_id not in by_id or by_id[term_id]["subject"] != "proband":
                raise ValueError("profile_terms_must_be_curated_proband_features")
    return config


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def paragraphs(path: Path) -> dict[int, str]:
    with zipfile.ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))
    return {
        index: "".join(node.text or "" for node in paragraph.findall(".//w:t", NS)).strip()
        for index, paragraph in enumerate(root.findall(".//w:p", NS), 1)
    }


def ontology_labels(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".obo" or text.startswith("format-version:"):
        labels: dict[str, str] = {}
        ident = ""
        for line in text.splitlines():
            if line.startswith("id: HP:"):
                ident = line.removeprefix("id: ").strip()
            elif ident and line.startswith("name: "):
                labels[ident] = line.removeprefix("name: ").strip()
                ident = ""
        return labels
    graph = json.loads(text)["graphs"][0]
    labels: dict[str, str] = {}
    for node in graph["nodes"]:
        node_id = node.get("id", "")
        if node_id.startswith("http://purl.obolibrary.org/obo/HP_"):
            labels[node_id.rsplit("/", 1)[-1].replace("HP_", "HP:")] = node.get("lbl", "")
    return labels


def vcf_sample(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                if len(samples) != 1:
                    raise ValueError("expected_single_sample_vcf")
                return samples[0]
    raise ValueError("vcf_header_missing")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--docx", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True, help="Private JSON curation and ordered profiles")
    parser.add_argument("--hpo", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    curation = config["curation"]
    profile_terms = config["profiles"]
    args.output.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(args.output, 0o700)

    source = paragraphs(args.docx)
    labels = ontology_labels(args.hpo)
    sample = vcf_sample(args.vcf)
    for row in curation:
        if labels.get(row["hpo_id"]) != row["label"]:
            raise ValueError(
                f"pinned_hpo_label_mismatch:{row['hpo_id']}:{labels.get(row['hpo_id'])!r}:{row['label']!r}"
            )
        missing = [number for number in row["feature_paragraphs"] if not source.get(number)]
        if missing:
            raise ValueError(f"source_paragraph_missing:{row['key']}:{missing}")

    with (args.output / "curation.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = (
            "key", "hpo_id", "label", "status", "onset", "subject", "scored_profiles",
            "source_paragraphs", "source_text_sha256", "original_wording", "justification",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in curation:
            texts = [source[number] for number in row["feature_paragraphs"]]
            writer.writerow({
                "key": row["key"], "hpo_id": row["hpo_id"], "label": row["label"],
                "status": "present", "onset": row["onset"], "subject": row["subject"],
                "scored_profiles": ",".join(
                    profile for profile, terms in profile_terms.items() if row["hpo_id"] in terms
                ) or "none_interpretive_input",
                "source_paragraphs": ",".join(map(str, row["feature_paragraphs"])),
                "source_text_sha256": ",".join(hashlib.sha256(text.encode()).hexdigest() for text in texts),
                "original_wording": " | ".join(texts), "justification": row["justification"],
            })
    os.chmod(args.output / "curation.tsv", 0o600)

    common_metadata = {
        "created": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "createdBy": "phenotype-sensitivity-preparer",
        "resources": [{
            "id": "hp", "name": "Human Phenotype Ontology", "namespacePrefix": "HP",
            "url": "https://exomiser.readthedocs.io/en/15.1.0/advanced_analysis.html",
            "version": "Exomiser phenotype data 2512 (bundled HPO graph)",
            "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
        }],
    }
    # Preserve the input term order, including any historical baseline order.
    profile_hashes: dict[str, str] = {}
    for profile, term_ids in profile_terms.items():
        profile_dir = args.output / profile
        profile_dir.mkdir(mode=0o700, exist_ok=True)
        by_id = {row["hpo_id"]: row for row in curation}
        features = [
            {"type": {"id": term_id, "label": by_id[term_id]["label"]}, "excluded": False}
            for term_id in term_ids
        ]
        packet = {
            "id": f"{profile}-phenotype-sensitivity-v1",
            "subject": {"id": sample},
            "phenotypicFeatures": features,
            "metaData": common_metadata,
        }
        atomic_json(profile_dir / "phenopacket.json", packet)
        ped = profile_dir / "family.ped"
        ped.write_text(f"FAMILY\t{sample}\t0\t0\t0\t2\n", encoding="utf-8")
        os.chmod(ped, 0o600)
        profile_hashes[profile] = sha256(profile_dir / "phenopacket.json")

    if "family_history" in config:
        atomic_json(args.output / "family-history.json", config["family_history"])
    atomic_json(args.output / "provenance.json", {
        "schema": "mva1_phenotype_sensitivity_provenance_v1",
        "document": {"path": str(args.docx), "sha256": sha256(args.docx)},
        "curation_config": {"path": str(args.config), "sha256": sha256(args.config)},
        "vcf": {"path": str(args.vcf), "sha256": sha256(args.vcf), "single_sample": True},
        "hpo_graph": {"path": str(args.hpo), "sha256": sha256(args.hpo)},
        "profiles": {
            profile: {"term_ids": terms, "phenopacket_sha256": profile_hashes[profile]}
            for profile, terms in profile_terms.items()
        },
        "negative_features": "none inferred",
        "contains_private_source_data": True,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
