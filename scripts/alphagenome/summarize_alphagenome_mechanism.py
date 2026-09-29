#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""Offline exon-masked RNA and returned-tissue inventory; no external calls."""
import json,os
import numpy as np,pandas as pd
from alphagenome_mechanism import REPORT,RESULT,MODES,js,tsv
from analyze_alphagenome_profiles import load

def summarize():
 os.umask(0o077)
 bub=json.loads((REPORT/'transcript_reconstruction.json').read_text())['exons']
 lz=json.loads((REPORT/'lztr1_transcript_coordinates.json').read_text())['exons']
 tasks=[('BUB1B_L737Ter',RESULT/'profiles/BUB1B_L737Ter/ref',RESULT/'profiles/BUB1B_L737Ter/alt',bub,'NM_001211.6'),('BUB1B_N1002K',RESULT/'profiles/BUB1B_N1002K/ref',RESULT/'profiles/BUB1B_N1002K/alt',bub,'NM_001211.6'),('hypothetical_cis',RESULT/'profiles/reference_sequence/sequence',RESULT/'profiles/cis_sequence/sequence',bub,'NM_001211.6'),('LZTR1_3prime_flank',RESULT/'noncoding_profiles/LZTR1_3prime_flank/ref',RESULT/'noncoding_profiles/LZTR1_3prime_flank/alt',lz,'NM_006767.4')]
 rows=[];inventory=[]
 for name,p,q,exons,tx in tasks:
  r,meta,rm=load(p,'RNA_SEQ');a,am,adm=load(q,'RNA_SEQ')
  assert rm['resolution']==1 and rm['interval']==adm['interval'] and r.shape==a.shape and meta.name.tolist()==am.name.tolist()
  mask=np.zeros(len(r),dtype=bool);start=rm['interval']['start']
  for e in exons:mask[e['genomic_start_1based']-1-start:e['genomic_end_1based']-start]=True
  for t,m in meta.iterrows():
   if m.get('strand') not in ['+','.']:continue
   ref=float(r[mask,t].mean());alt=float(a[mask,t].mean())
   rows.append(dict(scenario=name,transcript=tx,track_name=m['name'],biosample=m.get('biosample_name'),ontology_curie=m.get('ontology_curie'),exonic_bases=int(mask.sum()),ref_exonic_mean=ref,alt_exonic_mean=alt,alt_ref_ratio=alt/ref if ref else None,ln_ratio_with_0_001=float(np.log(alt+.001)-np.log(ref+.001)),interpretation='manual transcript exon mask; not identical to SDK GENCODE scorer; predicted not measured'))
  for mode in MODES:
   v,t,md=load(p,mode)
   for _,m in t.iterrows():inventory.append(dict(scenario=name,modality=mode,track_name=m['name'],ontology_curie=m.get('ontology_curie',''),biosample=m.get('biosample_name','tissue agnostic'),strand=m.get('strand','.'),values_returned=v.shape[0],resolution=md.get('resolution'),meaning='no returned junctions, not absence' if not len(v) else 'returned'))
 tsv(REPORT/'exon_masked_RNA_comparison.tsv',pd.DataFrame(rows));tsv(REPORT/'returned_tissue_inventory.tsv',pd.DataFrame(inventory))
 grouped=pd.DataFrame(rows).groupby('scenario').agg(min_alt_ref_ratio=('alt_ref_ratio','min'),max_alt_ref_ratio=('alt_ref_ratio','max'),track_count=('track_name','count')).reset_index()
 js(REPORT/'exon_masked_RNA_summary.json',grouped.to_dict(orient='records'))
 print(grouped.to_string(index=False))
if __name__=='__main__':summarize()
