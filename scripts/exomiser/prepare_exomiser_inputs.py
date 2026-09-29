#!/usr/bin/env python3
"""Create protected, source-linked Exomiser inputs without emitting phenotype data."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
TERMS = (
    ("HP:0000121", re.compile(r"nephrocalcinos", re.I), "documented", "unknown"),
    ("HP:0001622", re.compile(r"prematur|preterm", re.I), "documented", "unknown"),
    ("HP:0002859", re.compile(r"rhabdomyosarcoma", re.I), "documented", "unknown"),
    ("HP:0001518", re.compile(r"small\s+for\s+gestational\s+age|\bSGA\b", re.I), "documented", "unknown"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.chmod(temporary, 0o600)
    temporary.replace(path)


def docx_paragraphs(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        root = ElementTree.fromstring(archive.read("word/document.xml"))
    paragraphs = []
    for paragraph in root.findall(".//w:p", NS):
        text = "".join(node.text or "" for node in paragraph.findall(".//w:t", NS)).strip()
        if text:
            paragraphs.append(text)
    return paragraphs


def ontology_labels(path: Path) -> dict[str, str]:
    graph = json.loads(path.read_text(encoding="utf-8"))["graphs"][0]
    labels = {}
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
    parser.add_argument("--hpo", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--review-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.chmod(args.output, 0o700)

    paragraphs = docx_paragraphs(args.docx)
    labels = ontology_labels(args.hpo)
    features, evidence = [], []
    for term_id, pattern, certainty, onset in TERMS:
        matches = [(index + 1, text) for index, text in enumerate(paragraphs) if pattern.search(text)]
        if not matches:
            continue
        label = labels.get(term_id)
        if not label:
            raise ValueError(f"HPO term unavailable: {term_id}")
        features.append({
            "type": {"id": term_id, "label": label},
            "excluded": False,
            "evidence": [{"evidenceCode": {"id": "ECO:0006017", "label": "author statement from published clinical study"}}],
        })
        evidence.append({
            "hpo_id": term_id,
            "certainty": certainty,
            "onset": onset,
            "source": str(args.docx),
            "paragraphs": [number for number, _ in matches],
            "source_text_sha256": [hashlib.sha256(text.encode()).hexdigest() for _, text in matches],
            "source_text": [text for _, text in matches],
        })
    if not features:
        raise ValueError("no_documented_phenotype_mappings")

    sample = vcf_sample(args.vcf)
    packet = {
        "id": "mva1-local-phenopacket-v1",
        "subject": {"id": "mva1-local-subject"},
        "phenotypicFeatures": features,
        "metaData": {
            "created": "2026-09-09T00:00:00Z",
            "createdBy": "local-phenotype-curation",
            "resources": [{
                "id": "hp",
                "name": "Human Phenotype Ontology",
                "namespacePrefix": "HP",
                "url": "https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01",
                "version": "2026-09-01",
                "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
            }],
        },
    }
    atomic_json(args.output / "phenopacket.json", packet)
    atomic_json(args.output / "source-evidence.json", {
        "schema": "mva1_hpo_source_evidence_v1",
        "document_sha256": sha256(args.docx),
        "review_json_sha256": sha256(args.review_json),
        "ontology_sha256": sha256(args.hpo),
        "features": evidence,
        "uncertainties": {
            "subject_origin": "unknown",
            "specimen_origin": "unknown",
            "collection_timing": "unknown",
            "negative_phenotypes": "not_curated_without_explicit_documentation",
            "onset": "unknown_unless_explicitly_encoded",
        },
    })
    ped = args.output / "family.ped"
    ped.write_text(f"MVA1\t{sample}\t0\t0\t0\t2\n", encoding="utf-8")
    os.chmod(ped, 0o600)
    atomic_json(args.output / "input-provenance.json", {
        "schema": "mva1_exomiser_input_provenance_v1",
        "vcf_path": str(args.vcf),
        "vcf_sha256": sha256(args.vcf),
        "reference_build_id": "GRCh38",
        "phenopacket_sha256": sha256(args.output / "phenopacket.json"),
        "ped_sha256": sha256(ped),
        "single_sample_vcf": True,
        "sample_identifier_redacted_from_reports": True,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
