#!/usr/bin/env python3
"""Copy a VCF for PharmCAT, changing only non-finite QUAL values to missing."""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import os
from pathlib import Path


def sanitize(source: Path, output: Path, provenance: Path) -> dict:
    paths = [p.resolve() for p in (source, output, provenance)]
    if len(set(paths)) != 3:
        raise ValueError("input, output and provenance must be different files")
    if output.exists() or provenance.exists():
        raise FileExistsError("refusing to overwrite output or provenance")
    created = []
    counts = {"records": 0, "nonfinite_qual_replaced": 0,
              "finite_qual_retained": 0, "missing_qual_retained": 0}
    before = hashlib.sha256()
    after = hashlib.sha256()
    unchanged = hashlib.sha256()
    try:
        descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        created.append(output)
        with source.open("rb") as raw, os.fdopen(descriptor, "wb") as destination:
            compressed = raw.read(2) == b"\x1f\x8b"
            raw.seek(0)
            stream = gzip.GzipFile(fileobj=raw) if compressed else raw
            try:
                for line_number, line in enumerate(stream, 1):
                    before.update(line)
                    if not line.startswith(b"#"):
                        fields = line.split(b"\t")
                        if len(fields) < 8:
                            raise ValueError(f"invalid VCF field count at line {line_number}")
                        counts["records"] += 1
                        qual = fields[5]
                        if qual == b".":
                            counts["missing_qual_retained"] += 1
                        else:
                            try:
                                finite = Decimal(qual.decode("ascii")).is_finite()
                            except (InvalidOperation, UnicodeDecodeError):
                                raise ValueError(f"invalid QUAL at line {line_number}") from None
                            if finite:
                                counts["finite_qual_retained"] += 1
                            else:
                                fields[5] = b"."
                                counts["nonfinite_qual_replaced"] += 1
                        # This stream contains every non-QUAL byte, including
                        # sample fields, phasing, missingness and line endings.
                        unchanged.update(b"\t".join(fields[:5] + [b""] + fields[6:]))
                        line = b"\t".join(fields)
                    else:
                        unchanged.update(line)
                    destination.write(line)
                    after.update(line)
            finally:
                if compressed:
                    stream.close()
        result = {"status": "complete", "input": str(source.resolve()),
                  "output": str(output.resolve()), "input_compression": "gzip" if compressed else "none",
                  "counts": counts, "changed_field": "QUAL only; nonfinite to '.'",
                  "non_qual_fields_preserved": True, "genotypes_and_missingness_preserved": True,
                  "input_uncompressed_sha256": before.hexdigest(), "output_sha256": after.hexdigest(),
                  "non_qual_stream_sha256": unchanged.hexdigest(),
                  "interpretation": "Missing QUAL is not a new genotype or quality pass; upstream QC remains authoritative."}
        descriptor = os.open(provenance, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        created.append(provenance)
        with os.fdopen(descriptor, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        return result
    except BaseException:
        for path in created:
            path.unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    result = sanitize(args.input, args.output, args.provenance)
    print(json.dumps({"status": result["status"], "counts": result["counts"]}))


if __name__ == "__main__":
    main()
