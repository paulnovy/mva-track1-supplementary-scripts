#!/usr/bin/env python3
"""Query exact local SNVs from ZIP_STORED, tabix-indexed AlphaGenome tables.

The archive member is exposed read-only through a bounded loop device, avoiding a
full extraction. Candidate coordinates never appear in argv or stdout.
"""

import argparse
import csv
import json
import os
import stat
import struct
import subprocess
import zipfile
from pathlib import Path


IMAGE = "mva-mocha-prep:bcftools-1.20-samtools"


def protected_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    path.chmod(0o700)


def write_private(path: Path, data: str) -> None:
    path.write_text(data, encoding="utf-8")
    path.chmod(0o600)


def member_offset(archive: Path, info: zipfile.ZipInfo) -> int:
    with archive.open("rb") as fh:
        fh.seek(info.header_offset)
        header = fh.read(30)
    if len(header) != 30 or header[:4] != b"PK\x03\x04":
        raise RuntimeError(f"invalid local ZIP header in {archive.name}")
    filename_len, extra_len = struct.unpack_from("<HH", header, 26)
    return info.header_offset + 30 + filename_len + extra_len


def candidates(evidence: Path) -> list[dict]:
    doc = json.loads(evidence.read_text(encoding="utf-8"))
    out = []
    for item in doc["candidates"]:
        chrom, pos = item["coordinate"].rsplit(":", 1)
        alts = item["alt"]
        if len(alts) != 1 or len(item["ref"]) != 1 or len(alts[0]) != 1:
            raise RuntimeError("assigned AlphaGenome Atlas scope only supports exact SNVs")
        out.append({"chrom": chrom, "pos": int(pos), "ref": item["ref"], "alt": alts[0]})
    if len(out) != 2:
        raise RuntimeError("expected exactly two protected BUB1B candidates")
    return out


def run_table(name: str, archive: Path, work: Path, wanted: list[dict]) -> dict:
    with zipfile.ZipFile(archive) as zf:
        gz = next(i for i in zf.infolist() if i.filename.endswith(".tsv.gz"))
        tbi = next(i for i in zf.infolist() if i.filename.endswith(".tbi"))
        if gz.compress_type != zipfile.ZIP_STORED or tbi.compress_type != zipfile.ZIP_STORED:
            raise RuntimeError("archive members must be ZIP_STORED")
        offset = member_offset(archive, gz)
        index_path = work / f"{name}.tsv.gz.tbi"
        with zf.open(tbi) as src, index_path.open("wb") as dst:
            while block := src.read(1024 * 1024):
                dst.write(block)
        index_path.chmod(0o600)

    # Materialize only the selected compressed member. This is slower and needs
    # free disk equal to the BGZF member size, but avoids privileged loop devices
    # and host-specific containers.
    table_path = work / f"{name}.tsv.gz"
    with zipfile.ZipFile(archive) as zf, zf.open(gz) as src, table_path.open("wb") as dst:
        while block := src.read(8 * 1024 * 1024):
            dst.write(block)
    table_path.chmod(0o600)
    seqs = subprocess.check_output(["tabix", "-l", str(table_path)], text=True).splitlines()
    try:
        seqset = set(seqs)
        regions = []
        for item in wanted:
            chrom = item["chrom"]
            if chrom not in seqset and chrom.removeprefix("chr") in seqset:
                chrom = chrom.removeprefix("chr")
            elif chrom not in seqset and "chr" + chrom in seqset:
                chrom = "chr" + chrom
            if chrom not in seqset:
                raise RuntimeError("candidate contig is absent from the Atlas index")
            regions.append(f"{chrom}\t{item['pos']}\t{item['pos']}\n")
        region_path = work / f"{name}.regions.tsv"
        write_private(region_path, "".join(regions))
        raw_path = work / f"{name}.exact.tsv"
        with raw_path.open("wb") as out:
            subprocess.run(["tabix", "-R", str(region_path), str(table_path)], check=True, stdout=out)
        raw_path.chmod(0o600)
        header_path = work / f"{name}.header.txt"
        with header_path.open("wb") as out:
            subprocess.run(["tabix", "-H", str(table_path)], check=True, stdout=out)
        header_path.chmod(0o600)
    finally:
        table_path.unlink(missing_ok=True)

    rows = [line for line in raw_path.read_text(encoding="utf-8").splitlines() if line]
    return {
        "archive": archive.name,
        "member": gz.filename,
        "member_crc32": f"{gz.CRC:08x}",
        "member_size": gz.file_size,
        "index_member": tbi.filename,
        "indexed_contig_count": len(seqs),
        "queried_candidate_count": len(wanted),
        "regional_row_count": len(rows),
        "raw_result": raw_path.name,
        "header": header_path.name,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--evidence", type=Path, required=True)
    ap.add_argument("--splicing", type=Path, required=True)
    ap.add_argument("--avi", type=Path, required=True)
    ap.add_argument("--feature-importance", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    protected_dir(args.out)
    wanted = candidates(args.evidence)
    result = {
        "schema_version": 1,
        "scope": "exact protected BUB1B SNVs; local-only AlphaGenome Atlas lookup",
        "candidate_count": len(wanted),
        "tables": {},
    }
    result["tables"]["splicing"] = run_table("splicing", args.splicing, args.out, wanted)
    result["tables"]["avi"] = run_table("avi", args.avi, args.out, wanted)
    if args.feature_importance:
        result["tables"]["feature_importance"] = run_table(
            "feature_importance", args.feature_importance, args.out, wanted
        )
    write_private(args.out / "lookup-manifest.json", json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
