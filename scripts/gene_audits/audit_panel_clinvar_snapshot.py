#!/usr/bin/env python3
"""Exact called-allele matching of two normalized indexed VCFs in gene spans.

This is a snapshot lookup, not variant interpretation. Inputs must have been
normalized to the same assembly/reference and contig naming. Non-PASS calls
are retained; uncalled ALTs and reference genotypes do not constitute matches.
"""
import os
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def called_alts(record):
    return sorted({i for s in record.samples.values() for i in (s.get('GT') or ())
                   if i is not None and 0 < i <= len(record.alts or ())})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata', nargs='+', required=True)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--clinvar', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    os.umask(0o077)
    result = {'patient_vcf': args.vcf, 'clinvar_vcf': args.clinvar,
              'method': 'Exact CHROM/POS/REF/called ALT; same-reference normalized inputs; all FILTER states',
              'limitation': 'No match is not benign; not an exhaustive noncoding or structural-variant exclusion.',
              'genes': {}}
    with pysam.VariantFile(args.vcf) as patient, pysam.VariantFile(args.clinvar) as cv:
        if len(patient.header.samples) != 1:
            raise ValueError('Use a single-sample patient VCF')
        result['clinvar_header_provenance'] = [str(r).strip() for r in cv.header.records
            if str(r).startswith(('##fileDate=', '##source=', '##reference='))]
        for name in args.metadata:
            gene = json.loads(Path(name).read_text())
            if isinstance(gene, list):
                assert len(gene) == 1
                gene = gene[0]
            chrom = 'chr' + gene['seq_region_name'].removeprefix('chr')
            start, end = gene['start'], gene['end']
            clinical = defaultdict(list)
            classifications = Counter()
            region_records = list(cv.fetch(chrom, start-1, end))
            for r in region_records:
                if len(r.alts or ()) != 1:
                    raise ValueError('ClinVar input must be split biallelic')
                clinical[(r.pos, r.ref, r.alts[0])].append(r)
                classifications.update(r.info.get('CLNSIG', ()))
            matches = []
            called = 0
            for r in patient.fetch(chrom, start-1, end):
                for i in called_alts(r):
                    called += 1
                    alt = r.alts[i-1]
                    for c in clinical.get((r.pos, r.ref, alt), []):
                        matches.append({'variant': f'{chrom}:{r.pos}:{r.ref}>{alt}',
                            'filter': ';'.join(r.filter.keys()) or '.', 'qual': r.qual,
                            'samples': {s: {k: r.samples[s].get(k) for k in ('GT','DP','AD','GQ')}
                                        for s in r.samples},
                            'clinvar_id': c.id,
                            'clinvar': {k: c.info.get(k) for k in
                                ('ALLELEID','CLNSIG','CLNREVSTAT','CLNDN','GENEINFO','MC')}})
            result['genes'][gene['display_name']] = {
                'interval_1based_closed': [chrom,start,end],
                'clinvar_region_records': len(region_records),
                'clinvar_region_classifications': dict(classifications),
                'patient_called_alleles': called, 'matches': matches}
    with Path(args.output).open('x') as target:
        target.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({g: {'called': x['patient_called_alleles'],
                         'matched': len(x['matches'])} for g,x in result['genes'].items()}))


if __name__ == '__main__':
    main()
