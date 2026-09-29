#!/usr/bin/env python3
"""Descriptive depth/BAF from completed artifacts; not a mosaic-CNV assay."""
import argparse
import bisect
import collections
import gzip
import json
import math
import os
from pathlib import Path
import statistics

import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--regions', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--vcf', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()
    os.umask(0o077)
    out = Path(args.out)
    out.mkdir(exist_ok=True, parents=True)
    if (out / 'summary.json').exists():
        raise SystemExit('Existing summary; inspect checkpoint, do not repeat')
    autosomes = {f'chr{i}' for i in range(1, 23)}
    windows = []
    excluded = collections.Counter()
    with pysam.FastaFile(args.fasta) as fa, gzip.open(args.regions, 'rt') as f:
        for line in f:
            chrom, start, end, depth, *_ = line.split()
            if chrom not in autosomes:
                continue
            start, end, depth = int(start), int(end), float(depth)
            seq = fa.fetch(chrom, start, end).upper()
            n = sum(seq.count(x) for x in 'ACGT')
            frac = n / len(seq) if seq else 0
            if len(seq) < 500000 or frac < .99:
                excluded['short_or_gap_rich_window'] += 1
                continue
            gc = (seq.count('G') + seq.count('C')) / n
            if not .25 <= gc <= .65:
                excluded['extreme_gc_window'] += 1
                continue
            windows.append({'chrom': chrom, 'start': start, 'end': end,
                            'depth': depth / frac, 'gc': gc, 'acgt_fraction': frac})
    assert len(windows) > 100, 'Insufficient autosomal coverage bins'
    for w in windows:
        # A chromosome is not used to estimate its own GC-depth expectation.
        peers = [q for q in windows if q['chrom'] != w['chrom']]
        near = sorted(peers, key=lambda q: abs(q['gc'] - w['gc']))[:200]
        expected = statistics.median(q['depth'] for q in near)
        assert expected > 0
        w['relative_depth'] = w['depth'] / expected
    center = statistics.median(w['relative_depth'] for w in windows)
    bychrom = collections.defaultdict(list)
    for w in windows:
        w['relative_depth'] /= center
        w['baf'] = []
        bychrom[w['chrom']].append(w)
    for rows in bychrom.values():
        rows.sort(key=lambda w: w['start'])
    starts = {c: [w['start'] for w in rows] for c, rows in bychrom.items()}
    reasons = collections.Counter()
    with pysam.VariantFile(args.vcf) as vcf:
        assert len(vcf.header.samples) == 1
        sample = next(iter(vcf.header.samples))
        for rec in vcf:
            if rec.contig not in autosomes:
                continue
            reasons['autosomal_records_scanned'] += 1
            if len(rec.ref) != 1 or len(rec.alts or ()) != 1 or len(rec.alts[0]) != 1:
                reasons['not_biallelic_snv'] += 1
                continue
            if set(rec.filter) - {'PASS', '.'} or (rec.qual or 0) < 30:
                reasons['filter_or_qual'] += 1
                continue
            call = rec.samples[sample]
            if call.get('GT') not in ((0, 1), (1, 0)):
                reasons['not_fully_called_het'] += 1
                continue
            ad = call.get('AD')
            if not ad or len(ad) != 2 or any(x is None for x in ad):
                reasons['missing_ad'] += 1
                continue
            depth = sum(ad)
            if not 20 <= depth <= 120 or min(ad) < 3:
                reasons['depth_or_allele_support'] += 1
                continue
            pl = call.get('PL')
            if not pl or len(pl) != 3 or any(x is None for x in pl) or min(pl[0], pl[2]) - pl[1] < 20:
                reasons['weak_or_missing_genotype_likelihood'] += 1
                continue
            rows = bychrom.get(rec.contig, [])
            idx = bisect.bisect_right(starts.get(rec.contig, []), rec.start) - 1
            if idx < 0 or rec.start >= rows[idx]['end']:
                reasons['outside_retained_windows'] += 1
                continue
            rows[idx]['baf'].append(ad[1] / depth)
            reasons['accepted_hets'] += 1
    summary = {'status': 'complete', 'method': 'descriptive_1mb_gc_depth_and_called_het_baf_v1',
               'inputs': vars(args), 'excluded_windows': dict(excluded), 'vcf_qc': dict(reasons),
               'coverage_source_mapq': 10, 'chromosomes': {},
               'limitations': ['Exploratory relative depth, not absolute copy number or calibrated mosaic fraction.',
                   'No mappability/segmental-duplication mask, control cohort or cytogenetic validation.',
                   'GC expectation uses median of 200 nearest-GC windows on other autosomes.',
                   'Variant-only, genotype-selected heterozygotes cause BAF ascertainment bias; not a common-SNP panel.',
                   'Depth is reused from the completed MAPQ10 mosdepth run; no new high-MAPQ BAM pass.',
                   'Window dispersion is descriptive, not a confidence interval or detection limit.',
                   'A flat profile does not exclude heterogeneous, tissue-specific or low-fraction MVA.']}
    with (out / 'windows.tsv').open('w') as f:
        f.write('chrom\tstart\tend\tgc\tacgt_fraction\tdepth_x\trelative_depth\thet_count\tmedian_abs_baf_minus_half\n')
        for w in windows:
            bafdev = statistics.median(abs(x - .5) for x in w['baf']) if w['baf'] else None
            f.write('\t'.join(map(str, [w['chrom'], w['start'], w['end'], w['gc'], w['acgt_fraction'],
                        w['depth'], w['relative_depth'], len(w['baf']), bafdev])) + '\n')
    for chrom, rows in sorted(bychrom.items(), key=lambda x: int(x[0][3:])):
        ratios = sorted(w['relative_depth'] for w in rows)
        baf = [x for w in rows for x in w['baf']]
        hist = [0] * 20
        for x in baf:
            hist[min(19, math.floor(x * 20))] += 1
        summary['chromosomes'][chrom] = {'windows': len(rows), 'median_relative_depth': statistics.median(ratios),
                'window_p10': ratios[int(.1 * (len(ratios) - 1))], 'window_p90': ratios[int(.9 * (len(ratios) - 1))],
                'het_count': len(baf), 'baf_histogram_20_equal_bins': hist,
                'median_abs_baf_minus_half': statistics.median(abs(x - .5) for x in baf) if baf else None}
    tmp = out / 'summary.json.tmp'
    tmp.write_text(json.dumps(summary, indent=2) + '\n')
    tmp.replace(out / 'summary.json')
    print(json.dumps({'status': 'complete', 'windows': len(windows), 'accepted_hets': reasons['accepted_hets']}))


if __name__ == '__main__':
    main()
