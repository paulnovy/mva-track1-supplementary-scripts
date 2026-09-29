#!/usr/bin/env python3
"""Read-only REF compatibility audit; genomic exceptions stay in local output."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path

import pysam


def audit(vcf, fasta, aliases, output, expected_records):
    os.umask(0o077)
    output.mkdir(parents=True, exist_ok=False)
    mapping = dict(line.split() for line in aliases.read_text().splitlines())
    counts = Counter()
    started = datetime.now(timezone.utc).isoformat()
    with pysam.FastaFile(str(fasta)) as reference, pysam.VariantFile(str(vcf)) as variants, gzip.open(output / 'exceptions.tsv.gz', 'wt') as exceptions:
        exceptions.write('record_number\tcontig\tposition_1based\tREF\treference_contig\tobserved_REF\treason\n')
        for number, record in enumerate(variants, 1):
            target = mapping.get(record.contig)
            observed = ''
            if target is None:
                reason = 'unsupported_contig'
            elif not record.ref or any(base not in 'ACGTN' for base in record.ref.upper()):
                reason = 'invalid_REF'
            else:
                # VCF POS is 1-based; pysam start/fetch are 0-based. For a
                # symbolic SV compare the REF anchor, not the INFO/END span.
                observed = reference.fetch(target, record.start, record.start + len(record.ref)).upper()
                if len(observed) != len(record.ref):
                    reason = 'out_of_bounds'
                else:
                    reason = 'match' if observed == record.ref.upper() else 'REF_mismatch'
            counts[reason] += 1
            if reason != 'match':
                exceptions.write(f'{number}\t{record.contig}\t{record.pos}\t{record.ref}\t{target or ""}\t{observed}\t{reason}\n')
            if number % 500000 == 0:
                print(json.dumps({'records_checked': number}), flush=True)
    total = sum(counts.values())
    assert total == expected_records, 'Input record count differs from verified baseline'
    summary = {
        'status': 'complete', 'started_utc': started,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'records': total, 'expected_records': expected_records, 'counts': dict(counts),
        'all_records_REF_compatible': counts['match'] == total,
        'input_vcf': str(vcf), 'reference_fasta': str(fasta),
        'input_bytes': vcf.stat().st_size, 'input_modified_ns': vcf.stat().st_mtime_ns,
        'pysam_version': pysam.__version__,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'exceptions': 'exceptions.tsv.gz',
        'limitations': ['REF compatibility is not genotype validation or a biological conclusion.',
                        'Unsupported/mismatched records are retained in the original VCF; none are silently removed.'],
    }
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'status': 'complete', 'records': total, 'counts': dict(counts)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('vcf', 'fasta', 'aliases', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--expected-records', type=int, required=True)
    audit(**vars(parser.parse_args()))
