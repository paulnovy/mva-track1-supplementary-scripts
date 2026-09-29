#!/usr/bin/env python3
"""Audit regional SNP/indel evidence and read-pair phase without changing inputs.

All input/output coordinates are one-based VCF coordinates. Candidate SNPs are
screening hypotheses, not a replacement for a validated variant caller. Read
names remain in memory; outputs expose only aggregate fragment counts.
"""
import argparse
import itertools
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def base_calls(read):
    """Reference-oriented sequence, one-based positions, quality and end margin."""
    if read.query_sequence is None:
        return
    for q, r in read.get_aligned_pairs(matches_only=True):
        yield r + 1, read.query_sequence[q], read.query_qualities[q], min(q, read.query_length - q - 1)


def allele_at(read, pos, ref, alt):
    """Exact anchored simple SNP/indel allele; None for ambiguous/missing calls.

    CIGAR is intentionally not normalized across repeat-equivalent indels.
    Reference and insertion calls require an observed downstream flanking base.
    """
    if not read.query_sequence:
        return None
    pairs = read.get_aligned_pairs(matches_only=False)
    lookup = {r + 1: q for q, r in pairs if r is not None}
    q = lookup.get(pos)
    if q is None:
        return None
    if len(ref) == len(alt) == 1:
        if read.query_qualities[q] < 20 or min(q, read.query_length - q - 1) < 5:
            return None
        return ref if read.query_sequence[q] == ref else alt if read.query_sequence[q] == alt else None
    if not (ref.startswith(alt) or alt.startswith(ref)):
        return None
    end = pos + len(ref)
    qend = lookup.get(end)
    if qend is None or min(q, read.query_length - qend - 1) < 5:
        return None
    seq = read.query_sequence[q:qend]
    if not seq or min(read.query_qualities[q:qend + 1]) < 20:
        return None
    # For a deletion, absence of aligned query at deleted reference positions
    # plus a contiguous downstream flank distinguishes deletion from clipping.
    return ref if seq == ref else alt if seq == alt else None


def stat(path):
    s = Path(path).stat()
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def main():
    os.umask(0o077)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--vcf', required=True)
    p.add_argument('--region', required=True)
    p.add_argument('--target', action='append', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(mode=0o700, parents=True, exist_ok=False)
    chrom, span = a.region.split(':')
    start, end = map(int, span.split('-'))
    targets = {int(t.split(':')[0]): tuple(t.split(':')[1:]) for t in a.target}
    before = {x: stat(x) for x in [a.bam, a.bam+'.bai', a.vcf, a.vcf+'.tbi', a.fasta, a.fasta+'.fai']}
    bam = pysam.AlignmentFile(a.bam, 'rb')
    fa = pysam.FastaFile(a.fasta)
    reads = list(bam.fetch(chrom, start - 1, end))
    refseq = fa.fetch(chrom, start - 1, end)
    flags = Counter()
    tags = Counter()
    lengths, tlens = [], []
    # Each molecule has only one vote at a site. Discordant overlapping mates
    # are excluded, not adjudicated by picking whichever allele looks better.
    observations = {'strict': defaultdict(lambda: defaultdict(set)), 'exploratory': defaultdict(lambda: defaultdict(set))}
    strand = defaultdict(lambda: defaultdict(set))
    starts = defaultdict(lambda: defaultdict(set))
    primary = []
    for r in reads:
        for label, val in [('secondary', r.is_secondary), ('supplementary', r.is_supplementary), ('duplicate', r.is_duplicate), ('qcfail', r.is_qcfail), ('unmapped', r.is_unmapped), ('proper_pair', r.is_proper_pair)]:
            flags[label] += int(val)
        tags.update(k for k, v in r.get_tags())
        if r.query_length:
            lengths.append(r.query_length)
        if r.is_proper_pair and r.is_read1:
            tlens.append(abs(r.template_length))
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.is_unmapped:
            continue
        molecule = (r.get_tag('RG') if r.has_tag('RG') else '', r.query_name)
        if r.mapping_quality >= 20:
            primary.append((molecule, r))
        for pos, base, qual, margin in base_calls(r):
            if not start <= pos <= end or base not in 'ACGT':
                continue
            if r.mapping_quality >= 10 and qual >= 10:
                observations['exploratory'][pos][molecule].add(base)
            if r.mapping_quality >= 20 and qual >= 20 and margin >= 5:
                observations['strict'][pos][molecule].add(base)
                strand[pos][base].add(r.is_reverse)
                starts[pos][base].add((r.reference_start, r.is_reverse))
    counts = {mode: {pos: Counter(next(iter(v)) for v in mol.values() if len(v) == 1) for pos, mol in positions.items()} for mode, positions in observations.items()}
    vc = pysam.VariantFile(a.vcf)
    vcf_records = []
    candidates = {}
    for r in vc.fetch(chrom, start - 1, end):
        c = next(iter(r.samples.values()))
        gt = c.get('GT')
        vcf_records.append(dict(pos=r.pos, ref=r.ref, alt=list(r.alts or []), gt=gt, qual=r.qual, filters=list(r.filter), ad=c.get('AD'), pl=c.get('PL')))
        if gt and len(gt) == 2 and None not in gt and gt[0] != gt[1]:
            alleles = r.alleles
            # retain multiallelic 1/2 records as candidate pairs too
            pair = (alleles[gt[0]], alleles[gt[1]])
            candidates[r.pos] = dict(pos=r.pos, ref=pair[0], alt=pair[1], origin='core_vcf_het', vcf_ref=r.ref)
    raw_candidates = []
    for mode, positions in counts.items():
        for pos, c in positions.items():
            total = sum(c.values())
            ref = refseq[pos - start]
            for alt, n in c.items():
                if alt == ref or n < 3 or c[ref] < 3 or total < 10 or not .2 <= n / total <= .8:
                    continue
                row = dict(pos=pos, ref=ref, alt=alt, ref_count=c[ref], alt_count=n, total=total, mode=mode,
                           alt_strands=len(strand[pos][alt]), alt_distinct_starts=len(starts[pos][alt]))
                raw_candidates.append(row)
                if mode == 'strict' and len(strand[pos][alt]) == 2 and len(starts[pos][alt]) >= 3 and pos not in candidates:
                    candidates[pos] = dict(pos=pos, ref=ref, alt=alt, origin='raw_snp_screen_both_strands')
    for pos, (ref, alt) in targets.items():
        candidates.setdefault(pos, dict(pos=pos, ref=ref, alt=alt, origin='requested_target'))
    molecule_alleles = defaultdict(lambda: defaultdict(set))
    for mol, read in primary:
        for pos, candidate in candidates.items():
            if read.reference_start < pos <= read.reference_end:
                value = allele_at(read, pos, candidate['ref'], candidate['alt'])
                if value is not None:
                    molecule_alleles[mol][pos].add(value)
    candidate_counts = defaultdict(Counter)
    for m in molecule_alleles.values():
        for pos, values in m.items():
            if len(values) == 1:
                candidate_counts[pos].update(values)
    for pos, c in candidates.items():
        co = candidate_counts[pos]
        n = sum(co.values())
        c.update(fragment_counts=dict(co), accepted=(n >= 10 and co[c['ref']] >= 3 and co[c['alt']] >= 3 and .2 <= co[c['alt']] / n <= .8))
    edges = defaultdict(Counter)
    for m in molecule_alleles.values():
        valid = sorted(pos for pos, values in m.items() if len(values) == 1 and candidates[pos]['accepted'])
        for p1, p2 in itertools.combinations(valid, 2):
            edges[(p1,p2)][(next(iter(m[p1])), next(iter(m[p2])))] += 1
    edge_rows = [dict(pos1=p1,pos2=p2,distance=p2-p1,combination_counts={f'{a1}|{a2}':n for (a1,a2),n in c.items()},fragments=sum(c.values())) for (p1,p2),c in sorted(edges.items())]
    adjacency = defaultdict(set)
    for row in edge_rows:
        # Even one physically observed shared molecule is retained for the
        # connectivity audit; it is NOT asserted as a confident phase call.
        adjacency[row['pos1']].add(row['pos2'])
        adjacency[row['pos2']].add(row['pos1'])
    components = {}
    for pos in targets:
        seen, todo = {pos}, [pos]
        while todo:
            for nxt in adjacency[todo.pop()] - seen:
                seen.add(nxt)
                todo.append(nxt)
        components[pos] = sorted(seen)
    direct = defaultdict(Counter)
    for r in reads:
        for pos, (ref,alt) in targets.items():
            if r.reference_start < pos <= (r.reference_end or 0):
                for pp, base, qual, margin in base_calls(r):
                    if pp == pos:
                        category = 'strict' if not (r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) and r.mapping_quality >= 20 and qual >= 20 and margin >= 5 else 'filtered_or_low_quality'
                        direct[pos][category+':'+base] += 1
    def distribution(values):
        s = sorted(values)
        return {k: s[min(len(s)-1,int((len(s)-1)*v))] for k,v in [('min',0),('median',.5),('p99',.99),('max',1)]} if s else {}
    after = {x: stat(x) for x in before}
    assert before == after, 'Source metadata changed during audit'
    report = dict(region=a.region, coordinate_system='1-based inclusive', pysam_version=pysam.__version__, alignments=len(reads), flags=flags, tags=tags,
                  molecule_tags={k:tags[k] for k in ['BX','CB','RX','MI','HP','PS']}, read_length=distribution(lengths), proper_template_length=distribution(tlens),
                  strict_thresholds=dict(mapq=20,baseq=20,end_margin=5,ref_fragments=3,alt_fragments=3,depth=10,alt_balance=[.2,.8],dedup=True,qcfail_excluded=True,primary_only=True),
                  candidate_markers=list(sorted(candidates.values(),key=lambda c:c['pos'])),raw_snp_screen=raw_candidates,edges=edge_rows,target_components=components,
                  target_alignment_counts={p:dict(c) for p,c in direct.items()},target_fragment_base_counts={p:{mode:dict(counts[mode].get(p,{})) for mode in counts} for p in targets},
                  sources_unchanged=True,source_stat_before=before,source_stat_after=after,core_vcf_records=vcf_records)
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['alignments','molecule_tags','read_length','proper_template_length','target_components','target_fragment_base_counts','edges','sources_unchanged']},indent=2))


if __name__ == '__main__':
    main()
