#!/usr/bin/env python3
"""Count allele bases by independent template in two small regional SAM extracts."""
import argparse, json, os, re
from collections import defaultdict
from pathlib import Path

CIGAR=re.compile(r"(\d+)([MIDNSHP=X])")

def base_at(fields, site):
    ref=int(fields[3]); query=0
    for n,op in CIGAR.findall(fields[5]):
        n=int(n)
        if op in "M=X":
            if ref <= site < ref+n:
                i=query+site-ref
                return fields[9][i].upper(), ord(fields[10][i])-33
            ref += n; query += n
        elif op in "DN": ref += n
        elif op in "IS": query += n
    return None

def main(raw, out, sites):
    os.umask(0o077)
    out.mkdir(parents=True, exist_ok=True)
    result={}
    for site in sites:
        calls=defaultdict(list)
        for line in (raw/f"published-CEP192-{site}.sam").open():
            if line.startswith("@"): continue
            f=line.rstrip().split("\t")
            if len(f)<11 or f[9]=="*" or f[10]=="*": continue
            if int(f[1]) & 3844 or int(f[4]) < 20: continue
            hit=base_at(f,site)
            if hit and hit[1]>=20: calls[f[0]].append(hit[0])
        # A template is counted once only when all eligible reads agree; discordant templates retained separately.
        groups=defaultdict(int)
        for bases in calls.values():
            unique=set(bases)
            groups[next(iter(unique)) if len(unique)==1 else "discordant"] += 1
        ref,alt=sites[site]
        result[str(site)]={"ref":ref,"alt":alt,"eligible_templates":sum(groups.values()),"ref_templates":groups[ref],"alt_templates":groups[alt],"other_templates":sum(n for b,n in groups.items() if b not in (ref,alt,"discordant")),"discordant_templates":groups["discordant"],"all_template_base_counts":dict(sorted(groups.items()))}
    with (out/"published_cep192_allele_support.json").open("x") as target:
        target.write(json.dumps(result,indent=2)+"\n")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--mapped-alleles', required=True, type=Path)
    args = parser.parse_args()
    rows = json.loads(args.mapped_alleles.read_text())
    if any(r['chrom'] != 'chr18' or len(r['ref']) != 1 or len(r['alt']) != 1 for r in rows):
        parser.error('Expected mapped chr18 SNVs')
    main(args.raw, args.output, {r['pos']: (r['ref'], r['alt']) for r in rows})
