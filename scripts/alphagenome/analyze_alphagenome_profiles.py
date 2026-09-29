#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["alphagenome==0.9.0", "requests==2.32.5", "matplotlib>=3.9"]
# ///
"""Local, explicitly coordinate-aligned analysis of protected REF/ALT profiles."""
from pathlib import Path
import json,os
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from alphagenome_mechanism import REPORT,RESULT,OPS,MODES,TERMS,BUB,INTERVAL,private,js,tsv,sha
GENE_START=40160000;GENE_END=40225000

def load(p,mode):
 return (np.load(p/f'{mode}.npz')['values'].astype(float),pd.read_csv(p/f'{mode}.tracks.tsv',sep='\t'),json.loads((p/f'{mode}.metadata.json').read_text()))
def normalize_junction_coordinates(frame,metadata):
 frame=frame.copy();iv=metadata['interval']
 if frame.empty:return frame
 # SDK PredictSequence returns sequence-relative unlabelled junctions despite
 # attaching supplied genomic interval to TrackData. Native files stay intact.
 if frame.chrom.isna().all() or frame.chrom.fillna('').eq('').all():
  assert frame.start.min()>=0 and frame.end.max()<=iv['end']-iv['start']
  frame['chrom']=iv['chrom'];frame['start']+=iv['start'];frame['end']+=iv['start']
 return frame

def align_junctions(p,q):
 r,rt,rm=load(p,'SPLICE_JUNCTIONS');a,at,am=load(q,'SPLICE_JUNCTIONS');assert rt.name.tolist()==at.name.tolist()
 def read_junctions(path):
  try:return pd.read_csv(path,sep='\t')
  except pd.errors.EmptyDataError:return pd.DataFrame(columns=['chrom','start','end','strand'])
 rj=read_junctions(p/'SPLICE_JUNCTIONS.junctions.tsv');aj=read_junctions(q/'SPLICE_JUNCTIONS.junctions.tsv')
 rj=normalize_junction_coordinates(rj,rm);aj=normalize_junction_coordinates(aj,am)
 return align_arrays(rj,aj,r,a),rt

def align_arrays(rj,aj,r,a):
 key=['chrom','start','end','strand'];rr={tuple(x):i for i,x in enumerate(rj[key].itertuples(index=False,name=None))};aa={tuple(x):i for i,x in enumerate(aj[key].itertuples(index=False,name=None))};rows=[]
 assert len(rr)==len(rj) and len(aa)==len(aj)
 for k in sorted(rr.keys()|aa.keys()):
  for t in range(r.shape[1]):
   rv=float(r[rr[k],t]) if k in rr else None;av=float(a[aa[k],t]) if k in aa else None
   rows.append(dict(zip(key,k),track=t,ref=rv,alt=av,delta=None if rv is None or av is None else av-rv,returned_in_ref=k in rr,returned_in_alt=k in aa))
 return pd.DataFrame(rows,columns=key+['track','ref','alt','delta','returned_in_ref','returned_in_alt'])

def analyze_one(name,ref,alt,exons=None):
 rows=[];figdir=private(REPORT/'figures');curves={}
 for mode in MODES:
  if mode=='SPLICE_JUNCTIONS':continue
  r,rt,rm=load(ref,mode);a,at,am=load(alt,mode)
  assert rm['interval']==am['interval'] and rm['resolution']==am['resolution'] and rt.name.tolist()==at.name.tolist() and r.shape==a.shape
  start=rm['interval']['start'];res=rm['resolution'];coords=start+np.arange(len(r))*res
  mask=(coords>=GENE_START)&(coords<=GENE_END)
  # Track strand must match plus BUB1B or unstranded, never sum opposite strand.
  for t,track in rt.iterrows():
   if str(track.get('strand','.')) not in ['+','.']:continue
   rv=r[mask,t];av=a[mask,t];delta=av-rv;den=float(np.sum(rv));ratio=float(np.sum(av)/den) if den>0 else None
   rows.append(dict(scenario=name,modality=mode,track_name=track['name'],biosample=str(track.get('biosample_name','tissue agnostic')),ontology_curie=str(track.get('ontology_curie','')),resolution=res,ref_sum=float(np.sum(rv)),alt_sum=float(np.sum(av)),sum_ratio=ratio,max_abs_delta=float(np.abs(delta).max()) if len(delta) else 0,mean_abs_delta=float(np.abs(delta).mean()) if len(delta) else 0,n_bases_or_bins=len(rv)))
  curves[mode]=(coords,r,a,rt)
 jf,jt=align_junctions(ref,alt);jf=jf[(jf['strand']=='+')&(jf.start>=GENE_START)&(jf.end<=GENE_END)].copy()
 for col in ['name','biosample_name','ontology_curie']:
  jf[col]=jf.track.map(jt[col]) if col in jt else ''
 tsv(REPORT/f'{name}.junction_changes.tsv',jf)
 tsv(REPORT/f'{name}.track_changes.tsv',pd.DataFrame(rows))
 # Actual curves, same scale REF/ALT, plus separate deltas; one per selected context.
 present=curves['RNA_SEQ'][3]
 terms=[x for x in TERMS if x in set(present.ontology_curie)]
 for term in terms:
  fig,axes=plt.subplots(5,1,figsize=(13,10),sharex=True,gridspec_kw={'height_ratios':[2,1,1,1,1.7]})
  title=present.loc[present.ontology_curie==term,'biosample_name'].iloc[0]
  for ax,mode in zip(axes[:4],['RNA_SEQ','SPLICE_SITES','SPLICE_SITE_USAGE','DNASE']):
   x,r,a,meta=curves[mode];cand=meta.index[(meta.get('ontology_curie',pd.Series('',index=meta.index))==term)&meta.get('strand',pd.Series('.',index=meta.index)).isin(['+','.'])].tolist()
   if mode=='SPLICE_SITES':cand=meta.index[meta.get('strand',pd.Series('.',index=meta.index))=='+'].tolist()
   mask=(x>=GENE_START)&(x<=GENE_END)
   if not cand:ax.text(.05,.5,'Brak tego typu wyjścia dla tej tkanki',transform=ax.transAxes);continue
   # First compatible track to avoid silently averaging heterogeneous assays.
   t=cand[0];ax.plot(x[mask],r[mask,t],lw=.8,label='REF',color='#2b6cb0');ax.plot(x[mask],a[mask,t],lw=.8,label='ALT',color='#dd6b20',alpha=.8)
   ax.set_ylabel(mode,fontsize=8);ax.legend(fontsize=7,loc='upper right');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3))
  ax=axes[4];f=jf[jf.ontology_curie==term];f=f[f.returned_in_ref&f.returned_in_alt].copy();f['maxvalue']=f[['ref','alt']].max(axis=1);f=f.sort_values('maxvalue',ascending=False).head(12)
  if f.empty:ax.text(.03,.78,'Nie zwrócono połączeń — brak danych, nie wynik zerowy',transform=ax.transAxes,fontsize=8)
  for _,row in f.iterrows():
   w=row.end-row.start;maxk=max(f.maxvalue.max(),1)
   for field,color,sgn in [('ref','#2b6cb0',1),('alt','#dd6b20',-1)]:
    h=.15+.65*(row[field]/maxk);center=(row.start+row.end)/2;t=np.linspace(0,np.pi,100);xx=center-w/2*np.cos(t);yy=sgn*h*np.sin(t);ax.plot(xx,yy,color=color,lw=.8+1.5*row[field]/maxk,alpha=.7)

  exons=exons if exons is not None else json.loads((REPORT/'transcript_reconstruction.json').read_text())['exons']
  for exon in exons:
   ax.plot([exon['genomic_start_1based']-1,exon['genomic_end_1based']],[-.96,-.96],lw=4,color='black')
  ax.set_ylim(-1.08,1);ax.set_ylabel('Junctions\nREF↑ ALT↓',fontsize=8);ax.axhline(0,color='grey',lw=.5)
  for ax in axes:
   for v in BUB:ax.axvline(v['pos']-1,color='grey',ls=':',lw=.8)
   ax.set_xlim(GENE_START,GENE_END)
  axes[-1].set_xlabel(f"GRCh38 {rm['interval']['chrom']}, pozycja 0-based; bez uśredniania heterogenicznych torów")
  fig.suptitle(f'{name} — {title}\nPredykcje modelu, nie pomiar RNA; 1 Mb kontekstu, widok genu',fontsize=12)
  fig.tight_layout(rect=[0,0,1,.95]);p=figdir/f'{name}_{term.replace(":","_")}.png';fig.savefig(p,dpi=145);p.chmod(0o600);plt.close(fig)
 # Detail donor/acceptor and usage delta, original values not just scores.
 fig,axes=plt.subplots(3,1,figsize=(12,8),sharex=True)
 for ax,mode in zip(axes,['SPLICE_SITES','SPLICE_SITE_USAGE','DNASE']):
  x,r,a,meta=curves[mode];m=(x>=min(v['pos'] for v in BUB)-2500)&(x<=max(v['pos'] for v in BUB)+1400);tidx=meta.index[meta.get('strand',pd.Series('.',index=meta.index)).isin(['+','.'])]
  if len(tidx):
   d=(a-r)[:,tidx];ix=np.argmax(np.max(np.abs(d[m]),axis=0));t=tidx[ix];ax.plot(x[m],d[m,ix],color='#9b2c2c',lw=.7);ax.set_title(str(meta.loc[t,'name']),fontsize=8)
  ax.axhline(0,lw=.5,color='grey');ax.set_ylabel('ALT−REF\n'+mode,fontsize=8)
 for ax in axes:
  for v in BUB:ax.axvline(v['pos']-1,color='grey',ls=':')
 fig.suptitle(f'{name}: największa bezwzględna różnica w danej modalności\nSygnały z różnych skal — bez porównania wielkości między modalnościami');fig.tight_layout(rect=[0,0,1,.93]);p=figdir/f'{name}_local_deltas.png';fig.savefig(p,dpi=150);p.chmod(0o600);plt.close(fig)
 return rows,jf

def main():
 os.umask(0o077);allrows=[];checks={};jsum=[]
 for v in BUB:
  p=RESULT/'profiles'/v['id'];c=json.loads((p/'complete.json').read_text())
  for n,h in c['hashes'].items():assert sha(p/n)==h
  rows,jf=analyze_one(v['id'],p/'ref',p/'alt');allrows+=rows
  jsum.append(dict(scenario=v['id'],junction_rows=len(jf),max_abs_delta=float(jf.delta.abs().max()) if len(jf) else None,abs_delta_above_1=int((jf.delta.abs()>1).sum()),unmatched_junction_rows=int((~(jf.returned_in_ref&jf.returned_in_alt)).sum())))
 r=RESULT/'profiles/reference_sequence/sequence';a=RESULT/'profiles/cis_sequence/sequence';rows,jf=analyze_one('hypothetical_cis',r,a);allrows+=rows;jsum.append(dict(scenario='hypothetical_cis',junction_rows=len(jf),max_abs_delta=float(jf.delta.abs().max()) if len(jf) else None,abs_delta_above_1=int((jf.delta.abs()>1).sum()),unmatched_junction_rows=int((~(jf.returned_in_ref&jf.returned_in_alt)).sum())))
 tsv(REPORT/'profile_track_changes.tsv',pd.DataFrame(allrows));js(REPORT/'profile_junction_summary.json',jsum)
 # Validate sequence baseline equals reference output, or quantify discrepancies.
 baseline=[]
 for mode in MODES[:-2]+['DNASE']:
  x,_,_=load(RESULT/'profiles/BUB1B_L737Ter/ref',mode);y,_,_=load(r,mode);baseline.append(dict(modality=mode,shapes_equal=x.shape==y.shape,max_abs_difference=float(np.max(np.abs(x-y))) if x.shape==y.shape else None))
 js(REPORT/'reference_sequence_baseline_check.json',baseline);print(json.dumps({'complete':True,'track_comparisons':len(allrows),'figures':len(list((REPORT/'figures').glob('*.png'))),'junction_summary':jsum}))
if __name__=='__main__':main()
