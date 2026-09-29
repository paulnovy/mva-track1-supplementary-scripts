#!/usr/bin/env python3
"""Final bounded offline audit. Requires read access to protected child artifacts."""
import csv,json,os,stat,shutil
from pathlib import Path
import numpy as np
from alphagenome_mechanism import BASE,REPORT,RESULT,OPS,MODES,BUB,sha,js,utc

def rows(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
 os.umask(0o077)
 manifest=json.loads((OPS/'manifest.json').read_text())
 sources=dict(manifest['source_hashes'])
 tv=json.loads((OPS/'transcripts/verification.json').read_text())
 sources.update(tv['source_hashes_after'])
 source_checks={p:dict(expected=h,observed=sha(p)) for p,h in sources.items()}
 assert all(c['expected']==c['observed'] for c in source_checks.values())
 checkpoints=[]
 for directory in ['profiles','noncoding_scores','noncoding_profiles']:
  for p in sorted((RESULT/directory).glob('*/complete.json')):
   c=json.loads(p.read_text());assert c['status']=='complete'
   for rel,h in c['hashes'].items():assert sha(p.parent/rel)==h,(p,rel)
   checkpoints.append(dict(path=str(p),sha256=sha(p),checked_files=len(c['hashes'])))
 assert len(checkpoints)==7
 arrays=[]
 for p in sorted(RESULT.rglob('*.npz')):
  with np.load(p,allow_pickle=False) as z:
   shapes={}
   for k in z.files:
    a=z[k];assert np.isfinite(a).all(),(p,k);shapes[k]=list(a.shape)
  arrays.append(dict(path=str(p),arrays=shapes,finite=True))
 assert len(arrays)==52
 profile_metadata=[]
 for folder in ['profiles','noncoding_profiles']:
  for p in sorted((RESULT/folder).glob('*/*/*.metadata.json')):
   md=json.loads(p.read_text());mode=p.name.split('.')[0];assert mode in MODES
   assert md['shape'][1]==len(rows(p.with_name(mode+'.tracks.tsv')))
   if md['kind']=='track':assert md['shape'][0]*md['resolution']==2**20
   else:assert md['shape'][0]==len(rows(p.with_name(mode+'.junctions.tsv')))
   profile_metadata.append(str(p))
 assert len(profile_metadata)==40
 ts=json.loads((REPORT/'transcripts/summary.json').read_text())
 for n,h in ts['output_hashes'].items():assert sha(RESULT/'transcripts'/n)==h
 ref=(OPS/'reference_window.txt').read_text();cis=(OPS/'cis_window.txt').read_text();assert len(ref)==len(cis)==2**20
 changed=[i for i,(a,b) in enumerate(zip(ref,cis)) if a!=b]
 assert changed==[v['pos']-1-manifest['interval']['start'] for v in BUB]
 for v,i in zip(BUB,changed):assert ref[i]==v['ref'] and cis[i]==v['alt']
 table=rows(RESULT/'noncoding/noncoding_eligibility.tsv')
 technical=sum(r['technical_candidate']=='YES' for r in table);selected=sum(r['external_query_eligible']=='EXPLORATORY_SELECTED' for r in table)
 comparisons=len(rows(REPORT/'exon_masked_RNA_comparison.tsv'))
 figures=list((REPORT/'figures').glob('*.png'))
 skip=json.loads((REPORT/'exon17_skip_summary.json').read_text());assert isinstance(skip['tracks'],int) and skip['tracks']>=0
 assert json.loads((REPORT/'lztr1_full_profile_summary.json').read_text())['junction_rows']>=0
 assert shutil.disk_usage(BASE).free>=60*1024**3
 # Final outputs are data, including wrappers usable via bash; no executable bit is needed.
 for root in [REPORT,RESULT,OPS]:
  root.chmod(0o700)
  for p in root.rglob('*'):
   if p.is_dir():p.chmod(0o700)
   elif p.is_file():p.chmod(0o600)
 output_hashes={str(p):sha(p) for root in [REPORT,RESULT,OPS] for p in root.rglob('*') if p.is_file() and p.name not in ['verification.json']}
 v=dict(status='pass',at=utc(),sources_current_hash_match=source_checks,source_count_currently_rehashed=len(source_checks),reference_fasta_hash_inherited_not_rehashed=ts['input_hashes']['reference_fasta'],completed_api_checkpoints=checkpoints,numeric_arrays=arrays,profile_metadata_checked=len(profile_metadata),noncoding_rows=len(table),technical_candidates=technical,selected_exploratory_alleles=selected,exon_masked_comparisons=comparisons,figures=len(figures),transcript_products_checked=len(ts['output_hashes']),source_variants=BUB,reference_cis_exact_differences=2,permissions='all task directories0700/files0600',disk_free_gib=round(shutil.disk_usage(BASE).free/1024**3,2),source_or_pipeline_edits=False,output_hashes=output_hashes,limits=['LZTR1 zero returned junctions does not mean absent splicing','referenceSequence and PredictVariantREF not numerically identical: no crossmethod interaction claim','NMD and cis/trans remain conditional','reference full-file digest inherited; source VCF not newly whole-file rehashed in this bounded audit'])
 js(REPORT/'verification.json',v)
 print(json.dumps({k:v[k] for k in ['status','source_count_currently_rehashed','profile_metadata_checked','noncoding_rows','figures','disk_free_gib']}))
if __name__=='__main__':main()
