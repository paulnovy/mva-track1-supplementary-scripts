#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""Bounded protected molecular profiles. No patient identifiers or raw reads.

Use --mode prepare locally; run --mode predict only in a trusted,
credential-injecting runtime. Each task is checkpointed; a started-but-incomplete task is not restarted automatically.
The native SDK can retry transient RPC failures within a task; see METHODS.md.
"""
from __future__ import annotations
import argparse, datetime as dt, hashlib, json, os, sys, time, traceback
from pathlib import Path
import numpy as np
import pandas as pd

BASE=Path(os.environ.get('MVA_DATA_ROOT','.')).resolve()
REPORT=BASE/'reports/alphagenome-mechanism-20260910'
RESULT=BASE/'results/alphagenome-mechanism-20260910'
OPS=BASE/'operations/alphagenome-mechanism-20260910'
OLD=BASE/'results/alphagenome-api-predictions-20260910'
FASTA=BASE/'references/GRCh38_standard/GRCh38_standard.fa'
MODES=['RNA_SEQ','SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS','DNASE']
TERMS=['CL:0000047','UBERON:0000955','UBERON:0002113','CL:0002584','CL:0002551','CL:2000001','EFO:0003042','CL:1001606']
BUB=[dict(id='BUB1B_L737Ter',chrom='chr15',pos=40209701,ref='T',alt='G'),dict(id='BUB1B_N1002K',chrom='chr15',pos=40220612,ref='T',alt='G')]
INTERVAL=dict(chrom='chr15',start=39690867,end=40739443)

def utc(): return dt.datetime.now(dt.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def private(p):p.mkdir(parents=True,exist_ok=True);p.chmod(0o700);return p
def js(p,x):
 p=Path(p);private(p.parent);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n');t.chmod(0o600);t.replace(p)
def tsv(p,x):
 p=Path(p);private(p.parent);x.to_csv(p,sep='\t',index=False);p.chmod(0o600)
def fetch(chrom,start,end):
 idx={s.split('\t')[0]:s.split('\t')[1:] for s in Path(str(FASTA)+'.fai').read_text().splitlines()};length,offset,bases,width=map(int,idx[chrom][:4]);assert 0<=start<end<=length
 with open(FASTA,'rb') as f:
  a=offset+start//bases*width+start%bases;b=offset+(end-1)//bases*width+(end-1)%bases+1;f.seek(a);s=f.read(b-a).replace(b'\n',b'').replace(b'\r',b'').decode().upper()
 assert len(s)==end-start and set(s)<=set('ACGTN');return s

def prepare():
 for p in (REPORT,RESULT,OPS):private(p)
 entries=[];discovery=[];relevant=[];summ=[];checks={}
 for i,v in enumerate(BUB,1):
  for mode in MODES[:4]:
   p=OLD/f'variant-{i:02d}'/f'{mode}.scores.tsv.gz';checks[str(p)]=sha(p);f=pd.read_csv(p,sep='\t');f['analysis_allele']=v['id']
   assert np.isfinite(f.raw_score).all();g=f[f.gene_name=='BUB1B'] if 'gene_name' in f else f
   if g.empty and mode=='SPLICE_SITES':g=f
   ranked=f.reindex(f.raw_score.abs().sort_values(ascending=False).index).head(15);discovery.append(ranked)
   rg=g[g.ontology_curie.isin(TERMS)] if 'ontology_curie' in g else g;relevant.append(rg)
   b=g.iloc[g.raw_score.abs().argmax()] if len(g) else None
   summ.append(dict(allele=v['id'],modality=mode,all_rows=len(f),target_rows=len(g),max_abs_all=float(f.raw_score.abs().max()),max_abs_target=float(g.raw_score.abs().max()) if len(g) else None,top_target_tissue=None if b is None else str(b.get('biosample_name','tissue agnostic')),high_quantile_rows=int((f.quantile_score.abs()>.995).sum()),target_high_quantile_rows=int((g.quantile_score.abs()>.995).sum())))
 tsv(REPORT/'discovery_top_hits.tsv',pd.concat(discovery));tsv(REPORT/'relevant_tissue_scores.tsv',pd.concat(relevant));tsv(REPORT/'existing_score_summary.tsv',pd.DataFrame(summ))
 # Same reference window for both alleles and cis scenario; no diploid average.
 seq=fetch(INTERVAL['chrom'],INTERVAL['start'],INTERVAL['end']);alt=list(seq)
 for v in BUB:assert seq[v['pos']-1-INTERVAL['start']]==v['ref'];alt[v['pos']-1-INTERVAL['start']]=v['alt']
 (OPS/'reference_window.txt').write_text(seq);(OPS/'cis_window.txt').write_text(''.join(alt))
 for p in (OPS/'reference_window.txt',OPS/'cis_window.txt'):p.chmod(0o600)
 manifest=dict(schema='alphagenome_mechanism_v1',created_at=utc(),reference='GRCh38',sdk='0.9.0',interval=INTERVAL,variants=BUB,requested_outputs=MODES,ontology_terms=TERMS,source_hashes=checks,sequence_hashes={n:sha(OPS/n) for n in ['reference_window.txt','cis_window.txt']},tasks=['BUB1B_L737Ter','BUB1B_N1002K','reference_sequence','cis_sequence'],privacy='Only selected alleles/windows or public-reference-derived synthetic haplotype; no patient IDs/phenotypes/VCF/BAM',model_limit='No diploid prediction or determination of actual phase')
 js(OPS/'manifest.json',manifest);js(REPORT/'existing_score_summary.json',summ);print(json.dumps({'prepared':True,'score_tables':8,'tasks':4}))

def serialize(out,path):
 from alphagenome.data import junction_data
 private(path);stats={}
 for mode in MODES:
  d=getattr(out,mode.lower());assert d is not None
  a=np.asarray(d.values);assert np.isfinite(a).all() and a.ndim==2
  np.savez_compressed(path/f'{mode}.npz',values=a)
  tsv(path/f'{mode}.tracks.tsv',d.metadata)
  met=dict(shape=list(a.shape),finite=True,interval=None if d.interval is None else dict(chrom=d.interval.chromosome,start=d.interval.start,end=d.interval.end,strand=d.interval.strand),uns=str(d.uns))
  if isinstance(d,junction_data.JunctionData):
   tsv(path/f'{mode}.junctions.tsv',pd.DataFrame([dict(chrom=j.chromosome,start=j.start,end=j.end,strand=j.strand) for j in d.junctions]));met['kind']='junction'
  else:met.update(kind='track',resolution=int(d.resolution))
  js(path/f'{mode}.metadata.json',met);stats[mode]=met
 for p in path.iterdir():p.chmod(0o600)
 return stats

def make_client():
 from alphagenome_native_client import create_client
 return create_client(timeout_seconds=900)

def predict():
 if os.environ.get('MVA_ALLOW_NETWORK') != '1': raise RuntimeError('set MVA_ALLOW_NETWORK=1 for opt-in prediction')
 from alphagenome.models import dna_client
 from alphagenome.data import genome
 from run_alphagenome_api_stage import quiet_sdk as silence_dependencies
 m=json.loads((OPS/'manifest.json').read_text());assert m['variants']==BUB and m['interval']==INTERVAL
 for n,h in m['sequence_hashes'].items():assert sha(OPS/n)==h
 # Do not print credentials or transport exception details.
 with silence_dependencies():client=make_client()
 if not (RESULT/'metadata.complete.json').exists():
  with silence_dependencies():metadata=client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
  for field in metadata.__dataclass_fields__:
   d=getattr(metadata,field)
   if isinstance(d,pd.DataFrame):tsv(RESULT/'metadata'/f'{field}.tsv',d)
  js(RESULT/'metadata.complete.json',dict(at=utc(),model_version=str(client._model_version)))
 outputs=[dna_client.OutputType[x] for x in MODES];iv=genome.Interval(INTERVAL['chrom'],INTERVAL['start'],INTERVAL['end'])
 tasks=[(v['id'],v) for v in BUB]+[('reference_sequence',None),('cis_sequence',None)]
 for name,v in tasks:
  p=private(RESULT/'profiles'/name);done=p/'complete.json';started=p/'started.json'
  if done.exists():
   c=json.loads(done.read_text());assert c['manifest_sha256']==sha(OPS/'manifest.json')
   for n,h in c['hashes'].items():assert sha(p/n)==h
   continue
  if started.exists():raise RuntimeError('INCOMPLETE_REQUEST_REQUIRES_REVIEW')
  js(started,dict(at=utc(),manifest_sha256=sha(OPS/'manifest.json'),task=name));t=time.monotonic()
  with silence_dependencies():
   if v:
    allele=genome.Variant(v['chrom'],v['pos'],v['ref'],v['alt']);o=client.predict_variant(iv,allele,requested_outputs=outputs,ontology_terms=TERMS)
   else:
    s=(OPS/('reference_window.txt' if name=='reference_sequence' else 'cis_window.txt')).read_text();o=client.predict_sequence(s,interval=iv,requested_outputs=outputs,ontology_terms=TERMS)
  if v:stats={'ref':serialize(o.reference,p/'ref'),'alt':serialize(o.alternate,p/'alt')}
  else:stats={'sequence':serialize(o,p/'sequence')}
  hashes={str(q.relative_to(p)):sha(q) for q in p.rglob('*') if q.is_file()}
  js(done,dict(status='complete',task=name,manifest_sha256=sha(OPS/'manifest.json'),at=utc(),seconds=round(time.monotonic()-t,2),stats=stats,hashes=hashes));print(json.dumps({'completed':name,'seconds':round(time.monotonic()-t,1)}),flush=True)
 js(RESULT/'profiles.complete.json',dict(status='complete',at=utc(),tasks=[x[0] for x in tasks]));print('{"status":"complete","tasks":4}')

def main():
 os.umask(0o077);a=argparse.ArgumentParser();a.add_argument('--mode',choices=['prepare','predict'],required=True);o=a.parse_args()
 try:prepare() if o.mode=='prepare' else predict()
 except Exception as e:
  # Static class is sufficient for operator, detailed inputs never in stdout.
  js(OPS/'last-error.json',dict(at=utc(),error_class=type(e).__name__,status='failed',frames=[dict(file=Path(f.filename).name,line=f.lineno) for f in traceback.extract_tb(e.__traceback__)]));print(json.dumps({'status':'failed','error_class':type(e).__name__}));raise SystemExit(1)
if __name__=='__main__':main()
