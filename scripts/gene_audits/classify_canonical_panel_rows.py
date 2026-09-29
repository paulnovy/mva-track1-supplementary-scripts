#!/usr/bin/env python3
"""Classify saved VCF rows by canonical coding/intronic-splice geometry."""
import os
import argparse
import csv
import json
from pathlib import Path


def overlaps(span, interval) -> bool:
    return span[0] <= interval[1] and interval[0] <= span[1]


def called_alts(alt: str, gt: str):
    """Return called ALT alleles without changing the original VCF ALT/GT."""
    alts = alt.split(",")
    indices = []
    for allele in gt.replace("|", "/").split("/"):
        if allele.isdigit() and int(allele) > 0 and int(allele) <= len(alts):
            if int(allele) not in indices:
                indices.append(int(allele))
    return indices, [alts[index - 1] for index in indices]


def changed_span(pos: int, ref: str, alt: str):
    """Minimal changed reference span; insertions include their two flanks."""
    if not ref or not alt or any(base.upper() not in "ACGTN" for base in ref + alt):
        raise ValueError("Small-variant geometry requires sequence alleles; route symbolic/BND/star alleles to SV review")
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref_tail, alt_tail = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(ref_tail), len(alt_tail)) and ref_tail[-(suffix + 1)] == alt_tail[-(suffix + 1)]:
        suffix += 1
    if suffix:
        ref_tail, alt_tail = ref_tail[:-suffix], alt_tail[:-suffix]
    changed_start = pos + prefix
    if ref_tail:
        return changed_start, changed_start + len(ref_tail) - 1
    # In VCF an insertion is between bases; retain both adjacent coordinates.
    left = max(pos, changed_start - 1)
    return left, left + 1


def classify(spans, coding, splice):
    """Return coding first, then true intronic canonical donor/acceptor overlap."""
    if any(overlaps(span, interval) for span in spans for interval in coding):
        return "coding"
    roles = sorted({role for span in spans for role, interval in splice if overlaps(span, interval)})
    return "canonical_splice_intronic_" + "_".join(roles) if roles else "other_locus"


def load_geometry(coordinates: Path, structures: Path):
    geometry = {}
    with coordinates.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            gene = row[0]
            data = json.loads((structures / f"ensembl_{gene}_canonical_structure.json").read_text())
            translation = data["Translation"]
            cds_start, cds_end = sorted((translation["start"], translation["end"]))
            coding = []
            exons = sorted(((exon["start"], exon["end"]) for exon in data["Exon"]), key=lambda item: item[0])
            for start, end in exons:
                coding_start, coding_end = max(start, cds_start), min(end, cds_end)
                if coding_start <= coding_end:
                    coding.append((coding_start, coding_end))
            # Only true intronic flanks of internal exon junctions are splice windows.
            splice = []
            for (left_start, left_end), (right_start, right_end) in zip(exons, exons[1:]):
                if data["strand"] == 1:
                    splice.extend([("donor", (left_end + 1, left_end + 2)), ("acceptor", (right_start - 2, right_start - 1))])
                else:
                    splice.extend([("acceptor", (left_end + 1, left_end + 2)), ("donor", (right_start - 2, right_start - 1))])
            geometry[gene] = (coding, splice)
    return geometry


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--structures", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    os.umask(0o077)
    geometry = load_geometry(args.coordinates, args.structures)
    with args.records.open() as source, args.out.open("x", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=[*reader.fieldnames, "called_alt_indices", "called_alt", "allele_span", "canonical_region"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            pos = int(row["pos"])
            coding, splice = geometry[row["gene"]]
            indices, alts = called_alts(row["alt"], row["gt"])
            spans = [changed_span(pos, row["ref"], alt) for alt in alts]
            row["called_alt_indices"] = ",".join(map(str, indices))
            row["called_alt"] = ",".join(alts)
            row["allele_span"] = ",".join(f"{start}-{end}" for start, end in spans)
            row["canonical_region"] = classify(spans, coding, splice)
            writer.writerow(row)


if __name__ == "__main__":
    main()
