#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""One profile request to resolve empty gene-mask scores at LZTR1 3prime flank."""
import os,json,time,traceback
from pathlib import Path
from alphagenome_mechanism import RESULT,OPS,MODES,serialize,private,js,sha,utc,fetch
from run_alphagenome_api_stage import quiet_sdk

def run():
 if os.environ.get('MVA_ALLOW_NETWORK') != '1': raise RuntimeError('set MVA_ALLOW_NETWORK=1 for opt-in prediction')
 from alphagenome.models import dna_client
 from alphagenome.data import genome
 from alphagenome_native_client import create_client
 os.umask(0o077);p=private(RESULT/'noncoding_profiles/LZTR1_3prime_flank');done=p/'complete.json'
 if done.exists():
  c=json.loads(done.read_text())
  for n,h in c['hashes'].items():assert sha(p/n)==h
  print('{"status":"already_complete","requests":0}');return
 if (p/'started.json').exists():raise RuntimeError('INCOMPLETE_REQUEST_REQUIRES_REVIEW')
 assert fetch('chr22',20999038,20999039)=='G'
 v=genome.Variant('chr22',20999039,'G','A');iv=v.reference_interval.resize(2**20);terms=['CL:0000121','UBERON:0000955','UBERON:0002113','CL:0002584','CL:0002552']
 with quiet_sdk():client=create_client(timeout_seconds=900)
 js(p/'started.json',dict(at=utc(),variant='chr22:20999039:G>A',interval={'chrom':'chr22','start':iv.start,'end':iv.end},ontology_terms=terms,outputs=MODES,purpose='3prime context after empty splicing/polyA masks; exploratory not clinical'));t=time.monotonic()
 with quiet_sdk():r=client.predict_variant(iv,v,requested_outputs=[dna_client.OutputType[x] for x in MODES],ontology_terms=terms)
 stats=dict(ref=serialize(r.reference,p/'ref'),alt=serialize(r.alternate,p/'alt'));hashes={str(q.relative_to(p)):sha(q) for q in p.rglob('*') if q.is_file()}
 js(done,dict(status='complete',at=utc(),seconds=time.monotonic()-t,stats=stats,hashes=hashes));print('{"status":"complete","requests":1}')
if __name__=='__main__':
 try:run()
 except Exception as e:js(OPS/'lztr1-error.json',dict(at=utc(),error_class=type(e).__name__,frames=[dict(file=Path(f.filename).name,line=f.lineno) for f in traceback.extract_tb(e.__traceback__)]));print(json.dumps({'status':'failed','error_class':type(e).__name__}));raise SystemExit(1)
