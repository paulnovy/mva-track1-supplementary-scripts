#!/usr/bin/env python3
"""Prepare a public, sites-only ClinVar VCF for exact GRCh38 allele matching.

No sample data, genotype interpretation, REF replacement or variant filtering.
Unsupported contigs are retained unchanged in a separate BCF.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import pysam
import pysam.bcftools as bcftools
import pysam.version


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def prepare(source, reference, output):
    os.umask(0o022)  # Public resource, also readable by the existing worker UID.
    source, reference, output = map(Path, (source, reference, output))
    initial = (source.stat().st_size, source.stat().st_mtime_ns)
    output.mkdir(parents=True, exist_ok=False)
    counts = Counter()
    aliases = {}
    accepted = output / 'supported-original.bcf'
    unsupported = output / 'unsupported-original.bcf'
    rejected_hash = hashlib.sha256()
    expected = 0
    with pysam.FastaFile(str(reference)) as fasta, pysam.VariantFile(str(source)) as vcf:
        if len(vcf.header.samples):
            raise ValueError('Expected public sites-only ClinVar, not a sample VCF')
        for old in vcf.header.contigs:
            candidate = 'chrM' if old in ('MT', 'M') else 'chr' + old
            target = old if old in fasta.references else candidate
            if target in fasta.references:
                declared = vcf.header.contigs[old].length
                if declared is not None and declared != fasta.get_reference_length(target):
                    raise ValueError('Reference length mismatch for contig alias')
                aliases[old] = target
        if not aliases or len(set(aliases.values())) != len(aliases):
            raise ValueError('Contig aliases must be nonempty and one-to-one')
        # Runtime probe of the existing tabix index, not merely its existence.
        probe_contig = next(iter(aliases))
        next(vcf.fetch(probe_contig, 0, min(1000000, fasta.get_reference_length(aliases[probe_contig]))), None)
    with pysam.VariantFile(str(source)) as vcf, \
            pysam.VariantFile(str(accepted), 'wb', header=vcf.header) as good, \
            pysam.VariantFile(str(unsupported), 'wb', header=vcf.header) as other:
        for record in vcf:
            counts['input_records'] += 1
            if record.contig not in aliases:
                counts['unsupported_records'] += 1
                other.write(record)
                rejected_hash.update(str(record).encode())
            else:
                counts['supported_records'] += 1
                expected += max(1, len(record.alts or ()))
                good.write(record)
        assert counts['input_records'] == counts['supported_records'] + counts['unsupported_records']
    alias_file = output / 'contig-aliases.tsv'
    alias_file.write_text(''.join(f'{old}\t{new}\n' for old, new in aliases.items()))
    renamed = output / 'supported-renamed.bcf'
    normalized = output / 'normalized-unsorted.bcf'
    final = output / 'clinvar.normalized.vcf.gz'
    bcftools.annotate('--rename-chrs', str(alias_file), '-Ob', '-o', str(renamed),
                      str(accepted), catch_stdout=False)
    # Stop on every REF mismatch. Never use --check-ref s (allele replacement).
    bcftools.norm('-f', str(reference), '--check-ref', 'e', '-m', '-any',
                  '--old-rec-tag', 'CLINVAR_ORIGINAL_REP', '-Ob', '-o', str(normalized),
                  str(renamed), catch_stdout=False)
    bcftools.sort('-m', '256M', '-T', str(output / 'sort-tmp'), '-Oz', '-o', str(final),
                  str(normalized), catch_stdout=False)
    bcftools.index('--csi', str(final), catch_stdout=False)
    with pysam.VariantFile(str(final)) as result:
        counts['normalized_records'] = sum(1 for _ in result)
        result.fetch(next(iter(aliases.values())), 0, 1)
    if counts['normalized_records'] != expected:
        raise ValueError('Split ALT count conservation failed')
    check_hash = hashlib.sha256()
    with pysam.VariantFile(str(unsupported)) as result:
        for record in result:
            check_hash.update(str(record).encode())
    if check_hash.digest() != rejected_hash.digest():
        raise ValueError('Unsupported record conservation failed')
    if (source.stat().st_size, source.stat().st_mtime_ns) != initial:
        raise ValueError('Source changed during preparation')
    summary = {
        'status': 'complete', 'completed_utc': datetime.now(timezone.utc).isoformat(),
        'counts': dict(counts), 'expected_normalized_records': expected,
        'source': str(source), 'source_sha256': sha256(source),
        'reference': str(reference), 'reference_fai_sha256': sha256(str(reference) + '.fai'),
        'aliases': aliases, 'source_index_runtime_probe': True,
        'ref_mismatch_policy': 'abort; no REF swapping',
        'unsupported_preserved': True, 'sample_data_processed': False,
        'versions': {'pysam': pysam.__version__, 'bcftools': pysam.version.__bcftools_version__},
        'script_sha256': sha256(__file__),
        'outputs_sha256': {p.name: sha256(p) for p in (final, Path(str(final) + '.csi'), unsupported)},
        'limitations': ['Exact normalized allele matches only; no pathogenicity inference.',
                       'ClinVar is not a comprehensive population or transcript annotation source.',
                       'Unsupported contigs are not annotated by this prepared resource.'],
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    # Remove only files created here, after the final resource and quarantine are verified.
    for tmp in (accepted, renamed, normalized):
        tmp.unlink()
    print(json.dumps({'status': 'complete', 'counts': dict(counts)}), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'reference', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    prepare(**vars(parser.parse_args()))
