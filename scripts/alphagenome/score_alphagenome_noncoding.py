#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5"]
# ///
"""Two explicitly exploratory noncoding variants; no AF/clinical claims."""
import os,json,time,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from alphagenome_mechanism import BASE,REPORT,RESULT,OPS,private,js,tsv,sha,utc,fetch
from run_alphagenome_api_stage import quiet_sdk as silence_dependencies
V=[dict(id='BUB1B_intronic_VUS',chrom='chr15',pos=40192892,ref='C',alt='T',gene='BUB1B'),dict(id='LZTR1_3prime_flank',chrom='chr22',pos=20999039,ref='G',alt='A',gene='LZTR1')]
M=['RNA_SEQ','SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS','DNASE','POLYADENYLATION']
def run():
 if os.environ.get('MVA_ALLOW_NETWORK') != '1': raise RuntimeError('set MVA_ALLOW_NETWORK=1 for opt-in prediction')
 from alphagenome.models import dna_client,variant_scorers
 from alphagenome.data import genome
 from alphagenome_native_client import create_client
 os.umask(0o077);out=private(RESULT/'noncoding_scores');summary=[]
 with silence_dependencies():client=create_client(timeout_seconds=600)
 for v in V:
  assert fetch(v['chrom'],v['pos']-1,v['pos'])==v['ref'];p=private(out/v['id']);done=p/'complete.json'
  if done.exists():
   c=json.loads(done.read_text())
   for n,h in c['hashes'].items():assert sha(p/n)==h
   summary+=c['summary'];continue
  if (p/'started.json').exists():raise RuntimeError('INCOMPLETE_REQUEST_REQUIRES_REVIEW')
  js(p/'started.json',dict(at=utc(),selected_allele=v,frequency='unknown_not_rare',purpose='exploratory_functional_not_pathogenicity'))
  variant=genome.Variant(v['chrom'],v['pos'],v['ref'],v['alt']);interval=variant.reference_interval.resize(2**20)
  scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS[x] for x in M]
  with silence_dependencies():r=client.score_variant(interval=interval,variant=variant,variant_scorers=scorers,organism=dna_client.Organism.HOMO_SAPIENS)
  rows=[]
  assert len(r)==len(M)
  for mode,scorer,a in zip(M,scorers,r,strict=True):
   assert repr(a.uns['variant_scorer'])==repr(scorer)
   arrays={'raw':np.asarray(a.X)};arrays.update({k:np.asarray(x) for k,x in a.layers.items() if k is not None});np.savez_compressed(p/f'{mode}.npz',**arrays)
   tsv(p/f'{mode}.obs.tsv',a.obs.reset_index());tsv(p/f'{mode}.tracks.tsv',a.var.reset_index())
   f=variant_scorers.tidy_anndata(a,match_gene_strand=True)
   if f.empty:f=pd.DataFrame(columns=['raw_score','quantile_score','gene_name','output_type'])
   assert np.isfinite(f.raw_score.to_numpy(dtype=float)).all();tsv(p/f'{mode}.scores.tsv.gz',f)
   t=f[f.gene_name==v['gene']] if 'gene_name' in f else f
   ranked=f.loc[f.raw_score.abs().sort_values(ascending=False).index].head(15);tsv(p/f'{mode}.top15.tsv',ranked)
   row=dict(allele=v['id'],mode=mode,all_rows=len(f),target_rows=len(t),max_abs_all=float(f.raw_score.abs().max()) if len(f) else None,max_abs_target=float(t.raw_score.abs().max()) if len(t) else None,high_quantile_all=int((f.quantile_score.abs()>.995).sum()) if 'quantile_score'in f else None)
   rows.append(row)
  hashes={str(q.relative_to(p)):sha(q) for q in p.iterdir() if q.is_file()};js(done,dict(status='complete',at=utc(),hashes=hashes,summary=rows));summary+=rows;print(json.dumps({'completed':v['id']}),flush=True)
 js(REPORT/'noncoding_score_summary.json',summary);tsv(REPORT/'noncoding_score_summary.tsv',pd.DataFrame(summary));print('{"status":"complete","exploratory_variants":2}')
if __name__=='__main__':
 try:run()
 except Exception as e:js(OPS/'noncoding-error.json',dict(at=utc(),error_class=type(e).__name__,frames=[dict(file=Path(f.filename).name,line=f.lineno) for f in traceback.extract_tb(e.__traceback__)]));print(json.dumps({'status':'failed','error_class':type(e).__name__}));raise SystemExit(1)
