#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""Local-only LZTR1 profile figures with explicit unreturned junctions."""
import json
import analyze_alphagenome_profiles as m
m.GENE_START=20981296;m.GENE_END=21001032;m.BUB=[dict(pos=20999039)]
m.TERMS=['CL:0000121','UBERON:0000955','UBERON:0002113','CL:0002584','CL:0002552']
def main():
 c=json.loads((m.REPORT/'lztr1_transcript_coordinates.json').read_text());assert c['transcript']=='NM_006767.4' and c['chrom']=='chr22' and c['strand']=='+'
 p=m.RESULT/'noncoding_profiles/LZTR1_3prime_flank'
 rows,jf=m.analyze_one('LZTR1_3prime_flank',p/'ref',p/'alt',exons=c['exons'])
 m.js(m.REPORT/'lztr1_full_profile_summary.json',dict(track_rows=len(rows),junction_rows=len(jf),matched_junction_max_delta=float(jf.delta.abs().max()) if len(jf) else None,unmatched=int((~(jf.returned_in_ref&jf.returned_in_alt)).sum()) if len(jf) else 0,empty_junctions_means='not returned, not measured absence'))
if __name__=='__main__':main()
