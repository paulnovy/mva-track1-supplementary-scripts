#!/usr/bin/env python3
"""Summarise bounded short-read callability and existing DELLY evidence."""
from __future__ import annotations

import argparse, csv, json, os, re, statistics
from pathlib import Path

GENES = {
    "TRIM37": ("chr17", 58982633, 59106921),
    "CEP192": ("chr18", 12991283, 13125053),
}
def intervals_union(items):
    result=[]
    for start,end in sorted(items):
        if not result or start > result[-1][1]+1: result.append([start,end])
        else: result[-1][1]=max(result[-1][1],end)
    return result

def intersect(intervals, start, end):
    return [(max(a,start),min(b,end)) for a,b in intervals if a<=end and b>=start]

def transcript_intervals(doc, primary=False):
    selected=[]
    for tx in doc["Transcript"]:
        translation=tx.get("Translation")
        if tx.get("biotype") != "protein_coding" or not translation: continue
        mane=any(x.get("type") == "MANE_Select" for x in tx.get("MANE", []))
        if primary and not mane: continue
        cds=intersect([(x["start"],x["end"]) for x in tx["Exon"]], translation["start"], translation["end"])
        selected.extend(cds)
    return intervals_union(selected)

def load_depth(path):
    values={}
    with Path(path).open() as h:
        for row in csv.reader(h, delimiter="\t"):
            values[int(row[1])]=int(row[2])
    return values

def summarise(depths, intervals):
    pos=[p for a,b in intervals for p in range(a,b+1)]
    if not pos:
        raise ValueError("No coding intervals; check MANE metadata snapshot")
    out={"bases":len(pos)}
    for name, d in depths.items():
        x=[d.get(p,0) for p in pos]
        out[name]={"mean_depth":round(statistics.fmean(x),3), "median_depth":statistics.median(x),
                   "minimum_depth":min(x), "bases_ge_10":sum(v>=10 for v in x),
                   "bases_ge_20":sum(v>=20 for v in x),
                   "fraction_ge_10":round(sum(v>=10 for v in x)/len(x),6),
                   "fraction_ge_20":round(sum(v>=20 for v in x)/len(x),6)}
    q0=depths["mapq0"]; q20=depths["mapq20"]; q60=depths["mapq60"]
    out["mapping_sensitivity"]={"q20_lt_half_q0_bases":sum(q0.get(p,0)>=10 and q20.get(p,0)/q0[p]<.5 for p in pos),
                                "q60_lt_half_q20_bases":sum(q20.get(p,0)>=10 and q60.get(p,0)/q20[p]<.5 for p in pos)}
    return out

def overlaps(a,b,c,d): return a<=d and c<=b
def load_sv(path):
    results=[]
    endpoint=re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
    with Path(path).open() as h:
        fields="id chrom pos end svtype filter qual pe sr mapq srmapq gt gq ft dr dv rr rv alt consensus".split()
        for r in csv.DictReader(h, delimiter="\t", fieldnames=fields):
            chrom,pos,end=r["chrom"],int(r["pos"]),int(r["end"] if r["end"] not in ("", ".", None) else r["pos"])
            hits=[]
            for gene,(target_chrom,left,right) in GENES.items():
                if r["svtype"] == "BND":
                    mates=endpoint.findall(r["alt"])
                    if (chrom==target_chrom and left<=pos<=right) or any(c==target_chrom and left<=int(p)<=right for c,p in mates): hits.append(gene)
                elif chrom==target_chrom and overlaps(pos,end,left,right): hits.append(gene)
            if hits: results.append({**r,"target_genes":",".join(hits)})
    return results

def bin_summary(path):
    output={}
    with Path(path).open() as h:
        for row in csv.DictReader(h, delimiter="\t"):
            for gene,(chrom,start,end) in GENES.items():
                # local background: bins whose interval overlaps gene +/- 1 Mb
                if row["chrom"]==chrom and overlaps(int(row["start"]),int(row["end"]),max(1,start-1_000_000),end+1_000_000):
                    output.setdefault(gene,[]).append(row)
    return output

def span_bin_audit(path, hits):
    with Path(path).open() as handle:
        rows=list(csv.DictReader(handle, delimiter="\t"))
    answer=[]
    for hit in hits:
        gt = hit['gt']
        if not any(x.isdigit() and int(x) > 0 for x in gt.replace('|', '/').split('/')):
            continue
        kind = hit['svtype']
        if kind == 'BND':
            continue
        vid, gene, flt, chrom = hit['id'], hit['target_genes'], hit['filter'], hit['chrom']
        start, end = sorted((int(hit['pos']), int(hit['end'])))
        evidence = {k: hit[k] for k in ('pe','sr','mapq','srmapq','gq','ft','dr','dv','rr','rv')}

        clean=lambda r: r["clean"]=="1" and r["corrected_rd"]!="."
        inside=[r for r in rows if r["chrom"]==chrom and int(r["end"])>=start and int(r["start"])<=end and clean(r)]
        flanks=[r for r in rows if r["chrom"]==chrom and ((int(r["end"])<start and int(r["end"])>=start-2_000_000) or (int(r["start"])>end and int(r["start"])<=end+2_000_000)) and clean(r)]
        def stats(x):
            return {"clean_bins":len(x),"mean_corrected_rd":round(statistics.fmean(float(r["corrected_rd"]) for r in x),4),"median_abs_baf_minus_half":round(statistics.median(float(r["median_abs_baf_minus_half"]) for r in x if r["median_abs_baf_minus_half"] != "."),4) if any(r["median_abs_baf_minus_half"] != "." for r in x) else None} if x else None
        answer.append({"id":vid,"gene":gene,"type":kind,"filter":flt,"gt":gt,"coordinates":f"{chrom}:{start}-{end}","delly_evidence":evidence,"breakpoints_inside_gene":{g: [GENES[g][1] <= p <= GENES[g][2] for p in (start,end)] for g in gene.split(",")},"interpretation":"unvalidated balanced event; RD/BAF are not an inversion test" if kind == "INV" else "inspect normalized dosage and BAF; not an automatic CN classification","inside":stats(inside),"outside_clean_2mb":stats(flanks),"masked_or_uncallable_overlapping_bins":sum(1 for r in rows if r["chrom"]==chrom and int(r["end"])>=start and int(r["start"])<=end and not clean(r))})
    return answer

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trim-json', type=Path, required=True)
    p.add_argument('--cep-json', type=Path, required=True)
    p.add_argument('--raw', type=Path, required=True)
    p.add_argument('--bins', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--candidate-tsv', type=Path, help='Optional output of audit_transcript_panel_rows.py')
    a = p.parse_args()
    os.umask(0o077)
    # Fail before writing any files if this analysis was already produced.
    if (a.output / 'metrics.json').exists() or (a.output / 'sv_target_overlap.tsv').exists() or a.report.exists():
        p.error('Refusing to replace an existing report or analysis; select fresh output paths')
    docs = {'TRIM37': json.loads(a.trim_json.read_text()), 'CEP192': json.loads(a.cep_json.read_text())}
    for gene, (chrom, start, end) in GENES.items():
        doc = docs[gene]
        if doc.get('assembly_name') != 'GRCh38' or str(doc['seq_region_name']).removeprefix('chr') != chrom[3:] or (doc['start'], doc['end']) != (start, end):
            p.error(f'{gene}: metadata differs from the bounded GRCh38 audit; review scope explicitly')
    candidates = set()
    if a.candidate_tsv:
        with a.candidate_tsv.open() as handle:
            candidates = {int(r['pos']) for r in csv.DictReader(handle, delimiter='\t') if r['gene'] == 'CEP192' and ':coding' in r['primary_region']}
    payload = {'method': {'depth': 'samtools depth -aa -G UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -q 20 -Q {0,20,60}',
        'scope': 'Read-base depth; primary QC-pass nonduplicate nonsupplementary alignments. Overlapping mates are not collapsed for coverage; template counts are separate.'},
        'genes': {}, 'sv_hits': load_sv(a.raw / 'delly.all.tsv')}
    for gene, (chrom, start, end) in GENES.items():
        depths = {key: load_depth(a.raw / f'{gene}.gene.baseq20-{key}.tsv') for key in ('mapq0','mapq20','mapq60')}
        primary = transcript_intervals(docs[gene], True)
        all_cds = transcript_intervals(docs[gene])
        payload['genes'][gene] = {'coordinates': f'{chrom}:{start}-{end}', 'mane_cds_intervals': primary,
            'all_protein_coding_cds_union': all_cds, 'gene': summarise(depths, [(start,end)]),
            'mane_cds': summarise(depths, primary), 'all_coding_cds_union': summarise(depths, all_cds),
            'left_flank_50kb': summarise(depths, [(max(1,start-50000),start-1)]),
            'right_flank_50kb': summarise(depths, [(end+1,end+50000)])}
        if gene == 'CEP192':
            payload['cep192_candidate_base_depth'] = [{'position': pos, **{key: depths[key].get(pos,0) for key in depths}} for pos in sorted(candidates)]
    payload['masked_bins_local'] = bin_summary(a.bins)
    payload['alternate_genotype_span_bin_audit'] = span_bin_audit(a.bins, payload['sv_hits'])
    a.output.mkdir(parents=True, exist_ok=True)
    with (a.output / 'metrics.json').open('x') as handle:
        json.dump(payload, handle, indent=2)
        handle.write('\n')
    with (a.output / 'sv_target_overlap.tsv').open('x') as handle:
        fields = 'target_genes id chrom pos end svtype filter qual pe sr mapq srmapq gt gq ft dr dv rr rv alt consensus'.split()
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(payload['sv_hits'])
    lines = ['# TRIM37 / CEP192 callability and structural evidence', '',
        'Computed from the supplied files; no historical case conclusions are embedded in this producer.', '',
        '| target | bases | MAPQ20 mean | ≥20× |', '|---|---:|---:|---:|']
    for gene, values in payload['genes'].items():
        for label in ('gene','mane_cds','all_coding_cds_union'):
            q = values[label]
            lines.append(f"| {gene} {label} | {q['bases']} | {q['mapq20']['mean_depth']:.3f} | {q['mapq20']['fraction_ge_20']:.2%} |")
    lines += ['', '## All-filter structural calls', '',
        f"{len(payload['sv_hits'])} calls overlap a gene by span or either BND endpoint. All FILTER and GT states are retained in sv_target_overlap.tsv; reference/no-call genotypes are not alternate events.", '',
        '| event / type / FILTER / GT | clean inside bins, mean normalized RD, median abs(BAF−0.5) | outside 2-Mb bins |', '|---|---|---|']
    for row in payload['alternate_genotype_span_bin_audit']:
        def describe(x):
            return 'none' if x is None else f"{x['clean_bins']}, {x['mean_corrected_rd']}, {x['median_abs_baf_minus_half']}"
        lines.append(f"| {row['id']} / {row['type']} / {row['filter']} / {row['gt']} | {describe(row['inside'])} | {describe(row['outside_clean_2mb'])} |")
    lines += ['', '## Interpretation boundary', '',
        'Inspect normalized dosage and BAF together with masked bins and caller/read evidence. Near-baseline dosage can contradict a simple heterozygous deletion but cannot exclude every complex event. No automatic pathogenicity or copy-number classification is made.', '',
        'Inversions remain unvalidated by depth/BAF. The metrics explicitly record whether their reported breakpoints fall in a gene; an enclosing span is not direct gene-breakpoint disruption. A distant regulatory consequence is unknown.', '',
        'Published CEP192-site counts are in published_cep192_allele_support.json when that separate step was run. An absent alternate template is bounded negative evidence only. No reference result or no-call constitutes a gene exclusion.', '',
        'Depth is read-base callability, not calibrated independent-molecule copy number. Analysis of these saved inputs is not independent biological confirmation. Alternative transcripts, regulatory mechanisms, balanced/complex SVs, repeats and caller-missed alleles retain limitations.', '']
    a.report.parent.mkdir(parents=True, exist_ok=True)
    with a.report.open('x') as handle:
        handle.write('\n'.join(lines))

if __name__ == '__main__':
    main()
