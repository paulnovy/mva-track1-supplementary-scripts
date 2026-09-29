import argparse,os,resource,json,gzip,csv,hashlib,collections,time
from pathlib import Path
from datetime import datetime,timezone
parser=argparse.ArgumentParser(description="Annotate supplied SV calls using private phenotype terms")
parser.add_argument("--hpo-terms", type=Path, required=True, help="Private JSON object mapping HPO identifiers to labels")
args=parser.parse_args()
hponames=json.loads(args.hpo_terms.read_text())
if not isinstance(hponames,dict) or not hponames or not all(isinstance(k,str) and k.startswith("HP:") and isinstance(v,str) for k,v in hponames.items()):
 raise SystemExit("--hpo-terms must contain a nonempty HPO-to-label object")
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
os.umask(0o077)
base=Path(os.environ.get('MVA_DATA_ROOT','.')).resolve()
op=base/'operations/sv-gene-review-20260910'; out=base/'reports/sv-gene-review-20260910'
if out.exists(): raise SystemExit(f'refusing existing output directory: {out}')
paths={'bcf':base/'results/sv-screen-20260909/result/sv.delly.bcf','decoded':base/'operations/sv-review-20260910/delly.all.tsv','refseq':op/'ncbiRefSeqCurated.txt.gz','hpo':base/'results/hpo-exomiser-followup-20260909/resources/data/2512_phenotype/phenix/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt','exomiser':base/'results/hpo-exomiser-followup-20260909/work/rare-disease/native/MVA1-exomiser.genes.tsv'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(1024**2),b''):h.update(x)
 return h.hexdigest()
paths["private_hpo_terms"]=args.hpo_terms
before={str(p):sha(p) for p in paths.values()}
canonical={f'chr{i}' for i in range(1,23)}|{'chrX','chrY','chrM'}
# UCSC genePred coordinates are 0-based half-open.
genes={}; bins=collections.defaultdict(set)
for r in csv.reader(gzip.open(paths['refseq'],'rt'),delimiter='\t'):
 if r[2] not in canonical:continue
 key=(r[2],r[12]); g=genes.setdefault(key,{'gene':r[12],'chrom':r[2],'s':int(r[4]),'e':int(r[5]),'tx':[]})
 g['s']=min(g['s'],int(r[4]));g['e']=max(g['e'],int(r[5]))
 g['tx'].append({'id':r[1],'strand':r[3],'s':int(r[4]),'e':int(r[5]),'cds_s':int(r[6]),'cds_e':int(r[7]),'exons':list(zip(map(int,r[9].strip(',').split(',')),map(int,r[10].strip(',').split(','))))})
for k,g in genes.items():
 for i in range(g['s']//100000,(g['e']-1)//100000+1):bins[(g['chrom'],i)].add(k)
def query(c,s,e):
 if c not in canonical or e<=s:return set()
 hit=set()
 for i in range(max(0,s)//100000,(e-1)//100000+1):hit.update(bins.get((c,i),()))
 return {k for k in hit if genes[k]['s']<e and genes[k]['e']>s}
def merge(xs):
 res=[]
 for s,e in sorted(set(xs)):
  if s>=e:continue
  if res and s<=res[-1][1]:res[-1]=(res[-1][0],max(e,res[-1][1]))
  else:res.append((s,e))
 return res
hpo=collections.defaultdict(set)
for r in csv.reader(paths['hpo'].open(),delimiter='\t'):
 if len(r)==4 and r[3] in hponames:hpo[r[1]].add(r[3])
exorank={}
for r in csv.DictReader(paths['exomiser'].open(),delimiter='\t'):
 g=r['GENE_SYMBOL']; rank=int(r['#RANK']);exorank[g]=min(rank,exorank.get(g,99999))
mva={'BUB1B','BUB1','TRIP13','CEP57'}; renal_manual={'ATP6V1E1','PIK3C2A'}
focus=mva|renal_manual|{'LZTR1','FANCD2','GNRHR'}
target=focus|set(hpo)
for g in renal_manual:hpo[g].add('manually_selected_target')
fields=['chrom','pos','id','filter','svtype','end','chr2','pos2','svlen','pe','sr','mapq','srmapq','srq','precise','gt','gq','ft','rdcn','dr','dv','rr','rv']
rows=[]
for r in csv.reader(paths['decoded'].open(),delimiter='\t'):
 assert len(r)==len(fields)
 d=dict(zip(fields,r))
 for f in ['pos','end','pos2','svlen','pe','sr','mapq','srmapq','gq','dr','dv','rr','rv','rdcn']:
  d[f]=None if d[f]=='.' else int(d[f])
 d['eligible']=d['filter']=='PASS' and d['gt'] in ('0/1','1/1','1/0','0|1','1|1','1|0')
 rows.append(d)
eligible_count=sum(d['eligible'] for d in rows)
results=[]; candidates=[]; targetrows=[]; feature_counts=collections.Counter(); gene_counts=collections.Counter(); core=collections.defaultdict(collections.Counter)
for d in rows:
 c,p,t=d['chrom'],d['pos'],d['svtype']; end=d['end'] or p
 spans=[(c,max(0,p-1),max(p,end))] if t in ('DEL','DUP','INV') else []
 bp=[(c,max(0,p-1),p)]
 if t=='BND' and d['chr2']!='.' and d['pos2']:bp.append((d['chr2'],d['pos2']-1,d['pos2']))
 elif t in ('DEL','DUP','INV'):bp.append((c,max(0,end-1),end))
 ivkeys=set().union(*(query(*x) for x in spans)) if spans else set()
 bpkeys=set().union(*(query(*x) for x in bp))
 allkeys=ivkeys|bpkeys
 if not d['eligible']:allkeys={k for k in allkeys if k[1] in focus}
 effects=[]; cdsgenes=set();exongenes=set(); effgenes=set(); matched=set()
 for k in sorted(allkeys):
  g=genes[k]; intervals=spans if t in ('DEL','DUP') else bp
  ex=[];cds=[];txhits=[];exonlabels=[]
  for tx in g['tx']:
   txex=[];txcd=[]
   for j,(a,b) in enumerate(tx['exons']):
    exnum=j+1 if tx['strand']=='+' else len(tx['exons'])-j
    for cc,s,e in intervals:
     if cc!=g['chrom']:continue
     x,y=max(s,a),min(e,b)
     if x<y:txex.append((x,y));exonlabels.append(f"{tx['id']}:ex{exnum}")
     x,y=max(s,a,tx['cds_s']),min(e,b,tx['cds_e'])
     if x<y:txcd.append((x,y))
   if txex:txhits.append(tx['id'])
   ex.extend(txex);cds.extend(txcd)
  eb=sum(b-a for a,b in merge(ex));cb=sum(b-a for a,b in merge(cds))
  functional=k in bpkeys or t in ('DEL','DUP')
  if functional:effgenes.add(g['gene'])
  if eb:exongenes.add(g['gene'])
  if cb:cdsgenes.add(g['gene'])
  if functional and g['gene'] in target:matched.add(g['gene'])
  rel=('interval_dosage_hypothesis' if t in ('DEL','DUP') else 'breakpoint_in_gene' if k in bpkeys else 'inversion_interior_only')
  if g['gene'] in target:
   rec={**d,'gene':g['gene'],'gene_start_0':g['s'],'gene_end_0':g['e'],'relation':rel,'exon_overlap_bp':eb,'cds_overlap_bp':cb,'transcripts':'|'.join(sorted(set(txhits))),'exons':'|'.join(sorted(set(exonlabels))),'hpo':'|'.join(sorted(hpo[g['gene']])),'exomiser_rank':exorank.get(g['gene'],''),'MVA_gene':g['gene'] in mva}
   targetrows.append(rec)
   if d['eligible'] and functional:candidates.append(rec)
  if g['gene'] in focus:
   core[g['gene']]['all_records']+=1;core[g['gene']]['pass_nonref' if d['eligible'] else 'noneligible']+=1
   if d['eligible'] and cb:core[g['gene']]['pass_nonref_CDS']+=1
   if d['eligible'] and functional:core[g['gene']]['pass_nonref_effect_possible']+=1
 if d['eligible']:
  gs={k[1] for k in allkeys};size=abs(end-p) if t not in ('BND','INS') else abs(d['svlen'] or 0) if t=='INS' else 0
  rec={**d,'size_bp':size,'annotation_scope':('canonical' if c in canonical and (t!='BND' or d['chr2'] in canonical) else 'partial_canonical_primary_noncanonical_partner' if c in canonical else 'partial_noncanonical_primary_canonical_partner' if t=='BND' and d['chr2'] in canonical else 'noncanonical_not_annotated'),'genes_all':'|'.join(sorted(gs)),'effect_possible_genes':'|'.join(sorted(effgenes)),'coding_overlap_genes':'|'.join(sorted(cdsgenes)),'exon_overlap_genes':'|'.join(sorted(exongenes)),'phenotype_genes':'|'.join(sorted(matched)),'genes_n':len(gs),'coding_genes_n':len(cdsgenes)}
  results.append(rec)
  feature_counts['canonical' if c in canonical else 'noncanonical']+=1
  feature_counts['gene_overlap' if gs else 'no_gene_overlap']+=1
  if cdsgenes:feature_counts['coding_overlap']+=1
  if matched:feature_counts['phenotype_effect_possible']+=1
  if any(g in hpo and ('HP:0000121' in hpo[g] or 'manually_selected_target' in hpo[g]) for g in matched):feature_counts['renal_effect_possible']+=1
# Technical gate is exploratory; no pathogenicity assignment. PE/SR are different signal summaries, not independent molecules.
byid={r['id']:r for r in results}
for r in candidates:
 d=byid[r['id']];support=(r['pe'] or 0)+(r['sr'] or 0)
 good=support>=10 and max(r['mapq'] or 0,r['srmapq'] or 0)>=20 and (r['gq'] or 0)>=20 and r['ft']=='PASS'
 r['support_sum']=support;r['technical_gate']=good
 r['size_bp']=d['size_bp'];r['genes_n']=d['genes_n'];r['coding_genes_n']=d['coding_genes_n']
 r['large_span_flag']=r['svtype'] in ('DEL','DUP','INV') and d['size_bp']>=10000000
 r['copy_depth_inconsistent']=(r['svtype']=='DEL' and r['rdcn'] is not None and r['rdcn']>=2 or r['svtype']=='DUP' and r['rdcn'] is not None and r['rdcn']<=2) if r['chrom'] not in ('chrX','chrY','chrM') else False
 r['priority']=('MVA_core' if r['gene'] in mva else 'renal_direct' if 'HP:0000121' in hpo[r['gene']] or r['gene'] in renal_manual else 'rhabdomyosarcoma_direct' if 'HP:0002859' in hpo[r['gene']] else 'nonspecific_phenotype')
 r['review_class']=('broad_unvalidated_span' if r['large_span_flag'] else 'coding_or_breakpoint_candidate' if good and r['cds_overlap_bp'] else 'exonic_nonCDS_or_intragenic_candidate' if good else 'weak_technical_support')
# Keep all systematic candidates; short manual list prioritizes specific phenotype/coding, never substitutes the full table.
candidates.sort(key=lambda r:(r['priority']!='MVA_core',r['large_span_flag'],not r['technical_gate'],not bool(r['cds_overlap_bp']),r['priority']=='nonspecific_phenotype',r['copy_depth_inconsistent'],-r['support_sum'],r['id'],r['gene']))
manual=[r for r in candidates if r['technical_gate'] and r['cds_overlap_bp'] and not r['large_span_flag'] and r['priority']!='nonspecific_phenotype']
def write_tsv(path,rs,keys):
 with path.open('w') as f:
  w=csv.DictWriter(f,keys,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rs)
 os.chmod(path,0o600)
basefields=fields+['size_bp']
write_tsv(out/'pass_nonref.gene_annotation.tsv',results,basefields+['annotation_scope','genes_n','coding_genes_n','genes_all','effect_possible_genes','exon_overlap_genes','coding_overlap_genes','phenotype_genes'])
candfields=basefields+['gene','gene_start_0','gene_end_0','relation','exon_overlap_bp','cds_overlap_bp','transcripts','exons','hpo','exomiser_rank','MVA_gene','support_sum','technical_gate','genes_n','coding_genes_n','large_span_flag','copy_depth_inconsistent','priority','review_class']
write_tsv(out/'phenotype_candidates.tsv',candidates,candfields)
write_tsv(out/'specific_coding_candidates.tsv',manual,candfields)
write_tsv(out/'manual_gene_aware_shortlist.tsv',[r for r in candidates if r['review_class']=='coding_or_breakpoint_candidate'],candfields)
write_tsv(out/'core_gene_all_filters.tsv',[r for r in targetrows if r['gene'] in focus],fields+['gene','relation','exon_overlap_bp','cds_overlap_bp','transcripts','exons','hpo','MVA_gene'])
write_tsv(out/'gene_set_provenance.tsv',[{'gene':g,'hpo':'|'.join(sorted(hpo[g])),'MVA_core':g in mva,'manual_focus_target':g in focus,'annotated':any(k[1]==g for k in genes)} for g in sorted(target)],['gene','hpo','MVA_core','manual_focus_target','annotated'])
summary={'schema':'MVA1_SV_gene_aware_supplement_v1','status':'analysis_complete_pending_narrative','completed_utc':datetime.now(timezone.utc).isoformat(),'input_records':len(rows),'pass_nonreference':len(results),'annotation_gene_loci':len(genes),'annotation_unique_symbols':len(set(k[1] for k in genes)),'annotation_transcripts':sum(len(g['tx']) for g in genes.values()),'phenotype_target_genes':len(target),'phenotype_gene_counts':{h:sum(h in v for v in hpo.values()) for h in hponames},'counts':dict(feature_counts),'gene_variant_candidate_rows':len(candidates),'candidate_unique_variants':len(set(r['id'] for r in candidates)),'specific_coding_rows':len(manual),'specific_coding_unique_variants':len(set(r['id'] for r in manual)),'core_gene_counts':dict(core),'specific_coding_genes':sorted(set(r['gene'] for r in manual)),'source_hashes':before,'annotation_url':'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqCurated.txt.gz','annotation_server_last_modified':'2025-08-13T15:30:40Z','hpo_source_version':'Exomiser 2512 phenix ALL_SOURCES_ALL_FREQUENCIES; legacy snapshot; exact input HPO terms only','no_new_calling':True,'no_external_patient_query':True}
(out/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
after={str(p):sha(p) for p in paths.values()};assert before==after
ver={'checked_utc':datetime.now(timezone.utc).isoformat(),'source_hashes_before':before,'source_hashes_after':after,'unchanged':before==after,'input_records':len(rows),'pass_nonref_rows':len(results),'all_PASS_nonref_accounted':len(results)==eligible_count,'resources':{'affinity_cpu_max':2,'address_space_limit_gib':4,'single_thread':True},'notes':['Canonical contigs annotation only; noncanonical counted explicitly.','BND uses POS2/CHR2; no fictitious interval. INS breakpoint only.','INV interiors documented but not counted as functional overlap without breakpoint.','No individual genomic records sent to network.']}
(out/'verification.json').write_text(json.dumps(ver,indent=2)+'\n')
for p in list(op.iterdir())+list(out.iterdir()):
 if p.is_file():os.chmod(p,0o600)
print(json.dumps({k:summary[k] for k in ['counts','gene_variant_candidate_rows','candidate_unique_variants','specific_coding_rows','specific_coding_genes','core_gene_counts']},ensure_ascii=False))
