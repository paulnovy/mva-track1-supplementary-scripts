#!/usr/bin/env python3
"""Strict local normalization; quarantine incompatible records, account for each input."""
import argparse
from array import array
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import pysam
import pysam.bcftools as bcftools
import pysam.version


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def normalize(audit, aliases, output):
    os.umask(0o077)
    baseline = json.loads(audit.read_text())
    if baseline['status'] != 'complete':
        raise ValueError('Audit is incomplete')
    vcf, fasta = Path(baseline['input_vcf']), Path(baseline['reference_fasta'])
    initial_stat = vcf.stat()
    if (initial_stat.st_size, initial_stat.st_mtime_ns) != (
            baseline['input_bytes'], baseline['input_modified_ns']):
        raise ValueError('Source changed since audit')
    mapping = dict(line.split() for line in aliases.read_text().splitlines())
    if len(mapping.values()) != len(set(mapping.values())):
        raise ValueError('Contig aliases are not one-to-one')
    output.mkdir(parents=True, exist_ok=False)
    started = utc()
    counts = Counter()
    expected = array('I', [0])
    tag = 'MVA_SRC_ROW'
    supported = output / 'supported-original.bcf'
    quarantined = output / 'quarantined-original.bcf'
    quarantine_digest = hashlib.sha256()
    with pysam.FastaFile(str(fasta)) as reference, pysam.VariantFile(str(vcf)) as source:
        original_header = source.header.copy()
        if tag in source.header.info:
            raise ValueError('Provenance tag already exists')
        source.header.info.add(tag, 1, 'Integer', 'One-based record ordinal in source VCF')
        for old, new in mapping.items():
            length = reference.get_reference_length(new)
            if old not in source.header.contigs or source.header.contigs[old].length != length:
                raise ValueError('Alias length mismatch')
            if new not in source.header.contigs:
                source.header.contigs.add(new, length=length)
        with pysam.VariantFile(str(supported), 'wb', header=source.header) as accepted, \
                pysam.VariantFile(str(quarantined), 'wb', header=original_header) as rejected:
            for row, record in enumerate(source, 1):
                target = mapping.get(record.contig)
                if target is None:
                    reason = 'unsupported_contig'
                elif not record.ref or any(base not in 'ACGTN' for base in record.ref.upper()):
                    reason = 'invalid_REF'
                else:
                    observed = reference.fetch(target, record.start, record.start + len(record.ref)).upper()
                    reason = ('out_of_bounds' if len(observed) != len(record.ref) else
                              'match' if observed == record.ref.upper() else 'REF_mismatch')
                counts[reason] += 1
                expected.append(max(1, len(record.alts or ())) if reason == 'match' else 0)
                if reason == 'match':
                    # Preserve the original contig/genotype/alleles in this checkpoint.
                    record.info[tag] = row
                    accepted.write(record)
                else:
                    record.translate(original_header)
                    quarantine_digest.update(str(record).encode())
                    rejected.write(record)
                if row % 500000 == 0:
                    print(json.dumps({'stage': 'partition', 'records': row}), flush=True)
    if len(expected) - 1 != baseline['records'] or dict(counts) != baseline['counts']:
        raise ValueError('Partition disagrees with complete audit')
    # The accepted checkpoint retains original names. Alias renaming is explicit;
    # no liftover, REF swaps, duplicate removal or record filtering is performed.
    renamed = output / 'supported-renamed.bcf'
    bcftools.annotate('--rename-chrs', str(aliases), '-Ob', '-o', str(renamed),
                      str(supported), catch_stdout=False)
    normalized = output / 'normalized-unsorted.bcf'
    final = output / 'normalized.vcf.gz'
    print(json.dumps({'stage': 'bcftools_norm', 'input_records': counts['match']}), flush=True)
    bcftools.norm('-f', str(fasta), '--check-ref', 'e', '-m', '-any',
                  '--multi-overlaps', '.', '-Ob', '-o', str(normalized),
                  str(renamed), catch_stdout=False)
    bcftools.sort('-m', '512M', '-T', str(output / 'sort-tmp'), '-Oz',
                  '-o', str(final), str(normalized), catch_stdout=False)
    bcftools.index('--csi', str(final), catch_stdout=False)
    remaining = array('I', expected)
    total = 0
    with pysam.VariantFile(str(final)) as result:
        for record in result:
            row = record.info[tag]
            if not 0 < row < len(remaining) or remaining[row] == 0:
                raise ValueError('Unexpected or duplicate source lineage')
            remaining[row] -= 1
            total += 1
    if any(remaining) or total != sum(expected):
        raise ValueError('Normalized output lost source records or ALT alleles')
    actual_quarantine_digest = hashlib.sha256()
    quarantine_count = 0
    with pysam.VariantFile(str(quarantined)) as rejected:
        for record in rejected:
            actual_quarantine_digest.update(str(record).encode())
            quarantine_count += 1
    if (quarantine_count != baseline['records'] - counts['match'] or
            actual_quarantine_digest.digest() != quarantine_digest.digest()):
        raise ValueError('Quarantine conservation check failed')
    if (vcf.stat().st_size, vcf.stat().st_mtime_ns) != (
            initial_stat.st_size, initial_stat.st_mtime_ns):
        raise ValueError('Source changed during normalization')
    summary = {
        'status': 'complete', 'started_utc': started, 'finished_utc': utc(),
        'input_records': baseline['records'], 'partition_counts': dict(counts),
        'normalized_records': total, 'expected_normalized_records': sum(expected),
        'quarantined_records': quarantine_count, 'source_lineage_conserved': True,
        'quarantine_roundtrip_verified': True, 'source_metadata_unchanged': True,
        'pysam_version': pysam.__version__, 'htslib_version': pysam.__samtools_version__,
        'bcftools_version': pysam.version.__bcftools_version__,
        'audit_sha256': sha256(audit), 'aliases_sha256': sha256(aliases),
        'script_sha256': sha256(Path(__file__)), 'reference_fasta': str(fasta),
        'input_vcf': str(vcf),
        'outputs_sha256': {p.name: sha256(p) for p in
                           (supported, quarantined, final, Path(str(final) + '.csi'))},
        'policy': [
            'No liftover, REF swapping, deduplication, or genotype/QC filtering.',
            'All unsupported/invalid/out-of-bounds/mismatched records quarantined unchanged.',
            'Supported original BCF retains original contig names, GT, AD, PL and alleles.',
            'Normalized contigs use the explicit verified alias map; indels left-aligned.',
            'Multiallelic records split; other ALT alleles become missing GT, not reference.',
            'Split AD keeps REF and selected ALT depths; other ALT depths are not added to REF.',
            'DP and GQ remain source quantities; splitting does not re-genotype the sample.',
            'Normalization is not variant pathogenicity evidence or independent calling.',
        ],
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'status': 'complete', 'normalized_records': total,
                      'quarantined_records': quarantine_count}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('audit', 'aliases', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    normalize(**vars(parser.parse_args()))
