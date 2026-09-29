#!/usr/bin/env python3
"""Cross-reference an existing partial-GT audit with retained follow-up artifacts."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import os
from pathlib import Path
import re

import pysam
from audit_partial_gt import digest, gt_text, semantic, state, called_alts, write_json


def row_keys(row):
    chrom = 'chr' + row['chrom'].removeprefix('chr')
    return [(chrom, int(row['pos']), row['ref'], a)
            for a in row.get('called_alt', row['alt']).split(',')]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('audit', 'source', 'core-original', 'clinvar', 'reports', 'exomiser-log', 'output'):
        ap.add_argument('--' + name, type=Path, required=True)
    ap.add_argument('--panel', nargs='+', required=True, help='label=retained TSV')
    a = ap.parse_args()
    os.umask(0o077)
    a.output.mkdir(mode=0o700, exist_ok=False)
    inputs = [a.audit / 'partial-lineage.jsonl.gz', a.audit / 'clinvar-partial-candidates.json',
              a.audit / 'core-partial-special.json', a.source, a.core_original, a.clinvar, a.exomiser_log]
    panel_keys = defaultdict(list)
    for spec in a.panel:
        name, path = spec.split('=', 1)
        path = Path(path)
        inputs.append(path)
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                for k in row_keys(row):
                    panel_keys[k].append({'panel': name, 'path': str(path), 'row': row})
    matched = json.loads((a.audit / 'clinvar-partial-candidates.json').read_text())
    special = {tuple(x['key']): x for x in matched if x['categories'] != ['other']}
    core_special = json.loads((a.audit / 'core-partial-special.json').read_text())
    with pysam.VariantFile(str(a.source)) as vcf, pysam.VariantFile(str(a.core_original)) as original:
        for x in core_special:
            k = tuple(x['key'])
            target = special.setdefault(k, x)
            chrom, pos, ref, alt = k
            found = [r for r in vcf.fetch(chrom, pos-1, pos)
                     if r.pos == pos and r.ref == ref and r.alts == (alt,)]
            target['source_exact_calls'] = [{'state': state(r), 'call': semantic(r)} for r in found]
            tag = x['normalization_original']
            tag = ','.join(tag) if isinstance(tag, list) else tag
            oc, op, oref, oalts, index = tag.split('|')
            if oc not in original.header.contigs:
                oc = oc.removeprefix('chr')
            records = [r for r in original.fetch(oc, int(op)-1, int(op))
                       if r.pos == int(op) and r.ref == oref and r.alts == tuple(oalts.split(','))]
            assert len(records) == 1, 'core original lineage is ambiguous'
            r = records[0]
            target['core_original'] = {'locus': [r.contig, r.pos, r.ref, list(r.alts)],
                                       'original_alt_index': int(index), **semantic(r)}
            assert int(index) in r.samples[0]['GT']
    log = a.exomiser_log.read_text()
    blacklist = set(re.search(r'GeneBlacklistFilter\{\[([^\]]+)\]', log)[1].split(', '))
    coding, exomiser, panel_hits = [], [], []
    panel_counts = Counter()
    geometry = Counter()
    with gzip.open(a.audit / 'partial-lineage.jsonl.gz', 'rt') as stream:
        for line in stream:
            x = json.loads(line)
            k = tuple(x['key'])
            oc, op, ref, alts = x['original']['locus']
            original_key = ('chr' + oc.removeprefix('chr'), op, ref, alts[x['original_alt_index']-1])
            evidence = list(panel_keys.get(k, []))
            if k != original_key:
                evidence.extend(panel_keys.get(original_key, []))
            if evidence:
                panel_hits.append({'key': k, 'source_row': x['source_row'], 'evidence': evidence})
                for panel in {e['panel'] for e in evidence}:
                    panel_counts[panel] += 1
                    if x['coding_hits']:
                        panel_counts[panel + '_coding'] += 1
            if k in special:
                special[k]['panel_evidence'] = evidence
            if x['exomiser']:
                exomiser.append(x)
            if x['coding_hits']:
                genes = {g for g, tx in x['coding_hits']}
                x['blacklisted_genes'] = sorted(genes & blacklist)
                x['all_coding_genes_blacklisted'] = genes <= blacklist
                coding.append(x)
                group = 'blacklisted' if genes <= blacklist else 'not_blacklisted'
                geometry[group] += 1
                geometry[group + '_exomiser_retained'] += bool(x['exomiser'])
                if not x['exomiser'] and not genes <= blacklist:
                    cv = set(x.get('clinvar', {}).get('CLNSIG') or ())
                    label = ('benign' if cv & {'Benign', 'Likely_benign', 'Benign/Likely_benign'}
                             else 'conflicting' if 'conflicting' in x.get('categories', []) else 'no_exact_ClinVar')
                    geometry['nonblacklisted_unranked_' + label] += 1
    # Retrieve raw snapshot evidence for exact identities, without network use.
    with pysam.VariantFile(str(a.clinvar)) as vcf:
        for k, x in special.items():
            chrom, pos, ref, alt = k
            rows = [r for r in vcf.fetch(chrom, pos-1, pos)
                    if r.pos == pos and r.ref == ref and r.alts == (alt,)]
            assert len(rows) == 1
            x['snapshot_record_id'] = rows[0].id
            x['snapshot_all_info'] = dict(rows[0].info)
            x['full_GT_preserved'] = x.get('original', x.get('core_original'))
    # Coordinate mention is evidence of documentation only, never inferred review.
    report_files = sorted(x for x in a.reports.rglob('*') if x.is_file() and x.suffix in {'.md', '.tsv', '.json'}
                          and 'screening-integrity-audit-20260914' not in x.parts and x.stat().st_size <= 5_000_000)
    report_hits = defaultdict(list)
    positions = {str(k[1]): k[1] for k in special}
    positions.update({f'{k[1]:,}': k[1] for k in special})
    pattern = re.compile(r'(?<!\d)(?:' + '|'.join(re.escape(p) for p in positions) + r')(?!\d)')
    for path in report_files:
        for number, line in enumerate(path.read_text(errors='replace').splitlines(), 1):
            for token in set(pattern.findall(line)):
                report_hits[positions[token]].append({'path': str(path), 'line': number,
                                                       'text': line[:1500]})
    summary = Counter()
    for k, x in special.items():
        x['historical_report_coordinate_mentions'] = report_hits[k[1]]
        complete = (x.get('core_recovery') == 'complete_called' or
                    any(r['state'] not in ('missing_or_partial', 'reference') for r in x.get('source_exact_calls', [])))
        x['complete_classification_route_recovered'] = complete
        for c in x['categories']:
            summary[c + '_unique_affected'] += 1
            summary[c + '_not_complete_in_either_route'] += not complete
        summary['unique_affected'] += 1
        summary['complete_route_recovered'] += complete
        summary['not_complete_in_either_route'] += not complete
        summary['with_report_coordinate_mentions'] += bool(report_hits[k[1]])
        summary['with_later_panel_evidence'] += bool(x.get('panel_evidence'))
    write_json(a.output / 'clinical-review-queue.json', list(special.values()))
    write_json(a.output / 'coding-review-ledger.json', coding)
    write_json(a.output / 'exomiser-retained-partial.json', exomiser)
    write_json(a.output / 'panel-recovery.json', panel_hits)
    write_json(a.output / 'summary.json', {'clinical': dict(summary), 'coding': dict(geometry),
               'panels': dict(panel_counts), 'report_files_scanned': len(report_files),
               'report_scan_limit': 'coordinate mentions in MD/TSV/JSON <=5 MB; mentions are not proof of allele review'})
    write_json(a.output / 'provenance.json', {'inputs': {str(x): digest(x) for x in inputs},
               'report_files_scanned': [str(x) for x in report_files], 'script_sha256': digest(Path(__file__))})
    print(json.dumps({'clinical': summary, 'coding': geometry, 'panels': panel_counts}))


if __name__ == '__main__':
    main()
