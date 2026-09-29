#!/usr/bin/env python3
"""Offline split-GT lineage audit. Never re-genotype, call, align, or query APIs.

Preserve original records and compare a representation-only bcftools replay to
the retained normalized and annotated rows. All outputs require a fresh private
directory. Coding means changed-span/CDS overlap, not predicted consequence.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import sys

import pysam
import pysam.bcftools as bcftools

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'clinvar'))
from clinvar_exact import terms
from clinvar_triage import genotype_state, PLP
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'gene_audits'))
from audit_panel_clinvar_snapshot import called_alts
from classify_canonical_panel_rows import changed_span


def key(r):
    return r.contig, r.pos, r.ref, r.alts[0]


def gt_text(gt, phased=False):
    return ('|' if phased else '/').join('.' if x is None else str(x) for x in gt or ())


def state(r):
    return genotype_state(next(iter(r.samples.values())).get('GT'))


def semantic(r):
    s = next(iter(r.samples.values()))
    return {'GT': s.get('GT'), 'phased': s.phased,
            **{x: s.get(x) for x in ('AD', 'PL', 'DP', 'GQ')}}


def categories(r):
    t = terms(r.info.get('CLNSIG'))
    result = []
    if t & PLP:
        result.append('PLP')
    if t & {'uncertain_significance', 'uncertain_significance/likely_benign',
            'uncertain_significance/likely_pathogenic'}:
        result.append('VUS')
    if any('conflict' in x for x in t):
        result.append('conflicting')
    return result or ['other']


def varint(handle):
    value = 0
    for shift in range(0, 70, 7):
        b = handle.read(1)
        if not b:
            if not shift:
                return None
            raise ValueError('truncated protobuf varint')
        value |= (b[0] & 127) << shift
        if b[0] < 128:
            return value
    raise ValueError('invalid protobuf varint')


def proto_fields(handle):
    """Read the wire-0/2 fields used by Exomiser 15.1.0 jannovar.proto."""
    while (tag := varint(handle)) is not None:
        field, wire = tag >> 3, tag & 7
        if wire == 0:
            value = varint(handle)
            if value is None:
                raise ValueError('missing protobuf value')
        elif wire == 2:
            size = varint(handle)
            value = handle.read(size)
            if len(value) != size:
                raise ValueError('truncated protobuf message')
        else:
            raise ValueError(f'unexpected wire type {wire}')
        yield field, value


def interval(raw, names, lengths):
    f = dict(proto_fields(io.BytesIO(raw)))
    # Jannovar stores strand-oriented zero-based half-open positions.
    chrom_id = f[3]
    start, end = f.get(4, 0), f.get(5, 0)
    if f.get(2, 0) == 1:
        start, end = lengths[chrom_id] - end, lengths[chrom_id] - start
    return 'chr' + names[chrom_id].removeprefix('chr'), start + 1, end


def coding_bins(path):
    bins = defaultdict(set)
    transcripts = 0
    names, lengths = {}, {}
    with path.open('rb') as raw:
        if raw.read(4) != b'JTPB':
            raise ValueError('expected Exomiser JTPB transcript snapshot')
        with gzip.GzipFile(fileobj=raw) as stream:
            for number, model in proto_fields(stream):
                if number == 1:
                    for n, entry in proto_fields(io.BytesIO(model)):
                        d = dict(proto_fields(io.BytesIO(entry)))
                        if n == 2:
                            names[d[1]] = d[2].decode()
                        elif n == 3:
                            lengths[d[1]] = d[2]
                    continue
                if number != 2:
                    continue
                fields = list(proto_fields(io.BytesIO(model)))
                f = dict(fields)
                chrom, cds_start, cds_end = interval(f[4], names, lengths)
                if cds_start > cds_end:
                    continue
                transcripts += 1
                gene, accession = f[2].decode(), f[1].decode()
                for n, exon in fields:
                    if n != 5:
                        continue
                    ec, a, b = interval(exon, names, lengths)
                    assert ec == chrom
                    a, b = max(a, cds_start), min(b, cds_end)
                    if a <= b:
                        for bin_id in range(a // 10000, b // 10000 + 1):
                            bins[(chrom, bin_id)].add((a, b, gene, accession))
    print(json.dumps({'stage': 'CDS_geometry', 'coding_transcripts': transcripts}), flush=True)
    return bins, transcripts


def coding_hits(k, bins):
    chrom, pos, ref, alt = k
    try:
        a, b = changed_span(pos, ref, alt)
    except ValueError:
        return None
    hits = set()
    for bin_id in range(a // 10000, b // 10000 + 1):
        for start, end, gene, tx in bins.get((chrom, bin_id), ()):
            if a <= end and start <= b:
                hits.add((gene, tx))
    return sorted(hits)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, obj):
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2)
        stream.write('\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('original', 'normalized', 'source', 'core', 'fasta', 'aliases', 'transcripts', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--exomiser', type=Path, nargs='+', required=True)
    args = p.parse_args()
    os.umask(0o077)
    args.output.mkdir(mode=0o700, exist_ok=False)
    out = args.output
    paths = [args.original, args.normalized, args.source, args.core, args.fasta,
             args.aliases, args.transcripts, *args.exomiser, Path(__file__).resolve()]
    before = {str(x): (x.stat().st_size, x.stat().st_mtime_ns) for x in paths}
    bins, transcript_count = coding_bins(args.transcripts)
    originals = {}
    counts = Counter()
    original_gt = Counter()
    with pysam.VariantFile(str(args.original)) as vcf:
        if len(vcf.header.samples) != 1:
            raise ValueError('single sample required')
        with pysam.VariantFile(str(out / 'original-multiallelic.bcf'), 'wb', header=vcf.header) as checkpoint:
            for r in vcf:
                counts['original_supported_rows'] += 1
                s = next(iter(r.samples.values()))
                counts['original_missing_or_partial'] += state(r) == 'missing_or_partial'
                if len(r.alts or ()) <= 1:
                    continue
                row = r.info['MVA_SRC_ROW']
                assert row not in originals
                originals[row] = {'locus': [r.contig, r.pos, r.ref, list(r.alts)], **semantic(r)}
                original_gt[gt_text(s.get('GT'), s.phased)] += 1
                checkpoint.write(r)
    print(json.dumps({'stage': 'originals', 'sites': len(originals)}), flush=True)
    bcftools.annotate('--rename-chrs', str(args.aliases), '-Ob', '-o', str(out / 'renamed.bcf'),
                      str(out / 'original-multiallelic.bcf'), catch_stdout=False)
    bcftools.norm('-f', str(args.fasta), '--check-ref', 'e', '-m', '-any', '--multi-overlaps', '.',
                  '--old-rec-tag', 'AUDIT_ORIGINAL', '-Ob', '-o', str(out / 'split-replay.bcf'),
                  str(out / 'renamed.bcf'), catch_stdout=False)
    replay = {}
    with pysam.VariantFile(str(out / 'split-replay.bcf')) as vcf:
        for r in vcf:
            src = r.info['MVA_SRC_ROW']
            tag = r.info['AUDIT_ORIGINAL']
            tag = ','.join(tag) if isinstance(tag, tuple) else tag
            index = int(tag.split('|')[-1])
            assert (src, key(r)) not in replay
            replay[(src, key(r))] = (index, semantic(r))
    partial = {}
    normalized_verified = set()
    with pysam.VariantFile(str(args.normalized)) as vcf:
        for r in vcf:
            counts['normalized_rows'] += 1
            src = r.info['MVA_SRC_ROW']
            if src in originals:
                identity = (src, key(r))
                index, expected = replay[identity]
                assert semantic(r) == expected
                assert identity not in normalized_verified
                normalized_verified.add(identity)
            if state(r) != 'missing_or_partial':
                continue
            counts['normalized_partial_rows'] += 1
            if not called_alts(r):
                counts['normalized_partial_without_called_alt'] += 1
                continue
            counts['normalized_partial_called_alt'] += 1
            assert key(r) not in partial, 'duplicate partial allele: needs multiset comparison'
            index, expected = replay[(src, key(r))]
            original = originals[src]
            assert index in original['GT'], 'selected ALT was not called in source'
            partial[key(r)] = {'source_row': src, 'key': key(r), 'original_alt_index': index,
                              'original': original, 'split': semantic(r),
                              'filter': ';'.join(r.filter.keys()) or '.',
                              'coding_hits': coding_hits(key(r), bins), 'core': [], 'exomiser': []}
    assert normalized_verified == set(replay)
    route_counts = {}
    core_partial_special = []
    for route, path in [('source', args.source), ('core', args.core)]:
        rc = Counter()
        for_row = set()
        with pysam.VariantFile(str(path)) as vcf:
            for r in vcf:
                rc['rows'] += 1
                is_partial = state(r) == 'missing_or_partial'
                present = bool(called_alts(r))
                rc['partial_rows'] += is_partial
                rc['partial_called_alt'] += is_partial and present
                k = key(r)
                if route == 'core' and k in partial:
                    partial[k]['core'].append({'state': state(r), 'called_alt': present,
                                               'call': semantic(r), 'filter': ';'.join(r.filter.keys()) or '.',
                                               'normalization_original': r.info.get('NORMALIZATION_ORIGINAL')})
                if route == 'source' and k in partial:
                    assert semantic(r) == partial[k]['split']
                    assert r.info['MVA_SRC_ROW'] == partial[k]['source_row']
                    for_row.add(k)
                if 'ALLELEID' not in r.info:
                    continue
                cats = categories(r)
                rc['exact_matches'] += 1
                group = 'partial_called' if is_partial and present else ('complete_called' if present and not is_partial else 'not_called')
                for cat in cats:
                    rc[f'{group}_{cat}'] += 1
                if not (is_partial and present):
                    continue
                clinical = {x: r.info.get(x) for x in ('ALLELEID', 'CLNSIG', 'CLNREVSTAT', 'GENEINFO', 'MC', 'CLNDN')}
                if route == 'source':
                    partial[k]['clinvar'] = clinical
                    partial[k]['categories'] = cats
                elif cats != ['other']:
                    core_partial_special.append({'key': k, 'categories': cats, 'clinvar': clinical,
                                                 'call': semantic(r), 'also_source_partial': k in partial,
                                                 'normalization_original': r.info.get('NORMALIZATION_ORIGINAL')})
        if route == 'source':
            assert for_row == set(partial)
        route_counts[route] = dict(rc)
        print(json.dumps({'stage': route, 'counts': rc}), flush=True)
    exomiser_summary = {}
    for path in args.exomiser:
        rows = 0
        hits = set()
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                rows += 1
                k = ('chr' + row['CONTIG'].removeprefix('chr'), int(row['START']), row['REF'], row['ALT'])
                if k in partial:
                    hits.add(k)
                    partial[k]['exomiser'].append({'path': str(path), **row})
        exomiser_summary[str(path)] = {'rows': rows, 'partial_alleles': len(hits)}
    aggregate = Counter()
    with gzip.open(out / 'partial-lineage.jsonl.gz', 'xt') as stream:
        for k, row in partial.items():
            aggregate['filter_' + row['filter']] += 1
            aggregate['coding_overlap'] += bool(row['coding_hits'])
            aggregate['unsupported_geometry'] += row['coding_hits'] is None
            aggregate['clinvar_match'] += 'clinvar' in row
            aggregate['exomiser_recovered'] += bool(row['exomiser'])
            recovery = ('complete_called' if any(x['called_alt'] and x['state'] != 'missing_or_partial' for x in row['core'])
                        else 'partial_called' if any(x['called_alt'] for x in row['core'])
                        else 'not_called' if row['core'] else 'absent')
            row['core_recovery'] = recovery
            aggregate['core_' + recovery] += 1
            if row['coding_hits']:
                aggregate['coding_core_' + recovery] += 1
                aggregate['coding_exomiser_recovered'] += bool(row['exomiser'])
                aggregate['coding_filter_' + row['filter']] += 1
            stream.write(json.dumps(row) + '\n')
    matched = [x for x in partial.values() if 'clinvar' in x]
    write_json(out / 'clinvar-partial-candidates.json', matched)
    write_json(out / 'core-partial-special.json', core_partial_special)
    with (out / 'coding-partial.tsv').open('x') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['chrom', 'pos', 'ref', 'alt', 'source_row', 'original_GT', 'split_GT',
                         'genes', 'transcripts', 'filter', 'clinvar_categories', 'core_recovery', 'exomiser_recovered'])
        for k, x in partial.items():
            if not x['coding_hits']:
                continue
            writer.writerow([*k, x['source_row'], gt_text(x['original']['GT'], x['original']['phased']),
                             gt_text(x['split']['GT'], x['split']['phased']),
                             ','.join(sorted({g for g, t in x['coding_hits']})),
                             ','.join(t for g, t in x['coding_hits']), x['filter'],
                             ','.join(x.get('categories', [])), x['core_recovery'], bool(x['exomiser'])])
    after = {str(x): (x.stat().st_size, x.stat().st_mtime_ns) for x in paths}
    assert before == after
    summary = {'status': 'complete', 'counts': dict(counts), 'original_multiallelic_GT': dict(original_gt),
               'coding_transcripts': transcript_count, 'routes': route_counts, 'partial': dict(aggregate),
               'exomiser': exomiser_summary, 'proof': {'replay_GT_AD_PL_DP_GQ_phasing_matches': len(replay),
               'annotation_partial_semantics_matches': len(partial), 'input_metadata_unchanged': True},
               'method': 'CDS changed-span overlap across saved Ensembl transcripts; no consequence/rarity inference; exact normalized allele route joins',
               'versions': {'pysam': pysam.__version__, 'htslib': pysam.__samtools_version__}}
    write_json(out / 'summary.json', summary)
    # FASTA is very large: stat and existing normalization manifest pin it; hash other inputs.
    write_json(out / 'input-provenance.json', {str(x): {'bytes': before[str(x)][0], 'mtime_ns': before[str(x)][1],
               'sha256': digest(x) if x != args.fasta else None} for x in paths})
    write_json(out / 'output-sha256.json', {x.name: digest(x) for x in out.iterdir() if x.is_file()})
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
