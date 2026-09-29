#!/usr/bin/env python3
"""Collect standard bcftools genotype statistics without filtering source records."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re

import pysam.bcftools as bcftools
import pysam.version


def aggregate(text):
    counts, genotypes, columns = {}, [], None
    for line in text.splitlines():
        fields = line.split('\t')
        if fields[0] == 'SN':
            counts[fields[2].rstrip(':')] = int(fields[3])
        elif fields[0] == '# PSC':
            columns = [re.sub(r'^\[\d+\]', '', field) for field in fields[1:]]
        elif fields[0] == 'PSC':
            if columns is None or len(columns) != len(fields) - 1:
                raise ValueError('Unrecognized bcftools sample statistics columns')
            genotypes.append({key: float(value) if key == 'average depth' else int(value)
                              for key, value in zip(columns, fields[1:])
                              if key not in ('id', 'sample')})
    if counts.get('number of samples') != 1 or len(genotypes) != 1:
        raise ValueError('Expected one sample with genotype statistics')
    return {'site_counts': counts, 'genotype_statistics': genotypes[0]}


def run(normalization, output):
    os.umask(0o077)
    source = json.loads(normalization.read_text())
    if source['status'] != 'complete':
        raise ValueError('Normalization incomplete')
    output.mkdir(parents=True, exist_ok=False)
    result = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'bcftools_version': pysam.version.__bcftools_version__, 'datasets': {}}
    for filename, expected in [('supported-original.bcf', source['partition_counts']['match']),
                               ('normalized.vcf.gz', source['normalized_records'])]:
        path = normalization.parent / filename
        before = path.stat()
        with path.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != source['outputs_sha256'][filename]:
                raise ValueError('Input digest differs from normalization manifest')
        stats = bcftools.stats('-s', '-', str(path))
        (output / (filename + '.stats.txt')).write_text(stats)
        observed = aggregate(stats)
        if observed['site_counts']['number of records'] != expected:
            raise ValueError('Stats record count disagrees with normalization')
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError('Input changed during statistics')
        result['datasets'][filename] = observed
        print(json.dumps({'dataset': filename, 'records': expected, 'status': 'counted'}), flush=True)
    result.update(status='complete', finished_utc=datetime.now(timezone.utc).isoformat(),
                  source_files_unchanged=True, notes=[
                      'Statistics include all FILTER states; no variant/genotype selection performed.',
                      'PSC ref/hom/het counters describe SNPs, not all variant classes.',
                      'Normalized missing GT includes other ALT alleles made missing during splitting.',
                      'Original full genotypes remain authoritative for inheritance assessment.',
                      'Single-sample ALT-only VCF statistics are not population AF, karyotype or calling evidence.'])
    (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'complete', 'datasets': len(result['datasets'])}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--normalization', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(**vars(parser.parse_args()))
