#!/usr/bin/env python3
"""Bounded exact-sequence reconstruction and molecule-level CLCN locus audit.

All coordinates in tables are 1-based unless labelled boundary0. No genotypes
or sample reads leave the local host. Pysam is provided by the existing worker.
"""
import argparse, csv, hashlib, json, os, statistics, math
from collections import Counter, defaultdict
from pathlib import Path
import pysam

ROOT=Path('/out')
REF='/ref/GRCh38_standard.fa'
L=16050023 # zero-based boundary after chr1:16050023
R=16059762 # zero-based boundary after chr1:16059762
REGIONS={'locus':('chr1',15989999,16080000),'distant_control':('chr1',16999999,17015000)}

def rc(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def canon(s): return min(s,rc(s))
def dump(p,x): p.write_text(json.dumps(x,indent=2)+'\n')
def tsv(p,rows,cols=None):
    if not rows:return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=cols or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fasta(p,items):
    with p.open('w') as f:
        for n,s in items.items(): f.write('>'+n+'\n'+'\n'.join(s[i:i+80] for i in range(0,len(s),80))+'\n')
def primary(x):return not (x.is_secondary or x.is_supplementary or x.is_duplicate or x.is_qcfail)
def read_key(x):return x.query_name,1 if x.is_read1 else 2
def phash(x):return hashlib.sha256(x.encode()).hexdigest()

def prepare():
    ref=pysam.FastaFile(REF)
    seqs={n:ref.fetch(*v).upper() for n,v in REGIONS.items()}
    fasta(ROOT/'raw/reference-windows.fa',seqs)
    # Equal context per boundary allele; REF windows and ALT all 2000 bp.
    hs={'DEL633_join':ref.fetch('chr1',L-1000,L)+ref.fetch('chr1',R,R+1000),
        'REF_left':ref.fetch('chr1',L-1000,L+1000),
        'REF_right':ref.fetch('chr1',R-1000,R+1000)}
    # Larger full-locus contigs used to retain every local paralog competitor.
    hs.update({'locus_reference':seqs['locus']})
    fasta(ROOT/'analysis/candidate-haplotypes.fa',hs)
    shifts=[]
    seed=hs['DEL633_join'][900:1100]
    for d in range(-200,201):
        s=ref.fetch('chr1',L+d-100,L+d)+ref.fetch('chr1',R+d,R+d+100)
        # Compare same fixed outer contexts, not translated sequences.
        fixed=ref.fetch('chr1',L-300,L+d)+ref.fetch('chr1',R+d,R+300)
        original=ref.fetch('chr1',L-300,L)+ref.fetch('chr1',R,R+300)
        if fixed==original:shifts.append(d)
    dump(ROOT/'analysis/breakpoint.json',{'left_boundary0':L,'right_boundary0':R,'removed_length_bp':R-L,
        'equivalent_shifts':shifts,'left_ambiguous_boundary0':[L+min(shifts),L+max(shifts)],
        'right_ambiguous_boundary0':[R+min(shifts),R+max(shifts)],'junction_200bp':seed})
    probes=[]
    for k in (31,51):
        for n in ('DEL633_join','REF_left','REF_right'):
            for offset in range(-(k-10),-9):
                s=hs[n][1000+offset:1000+offset+k].upper()
                probes.append({'id':f'{n}_k{k}_o{offset}','kind':n,'k':k,'chrom':'chr1','start1':0,'sequence':s,'canonical':canon(s),'gc':round(sum(b in 'GC' for b in s)/k,4),'window_gc':0})
        for n,(chrom,a,b) in REGIONS.items():
            for off in range(250,len(seqs[n])-250-k,25):
                s=seqs[n][off:off+k]; context=seqs[n][off-250:off+k+250]
                if set(s)<=set('ACGT'):
                    probes.append({'id':f'{n}_{a+off+1}_k{k}','kind':'dosage','k':k,'chrom':chrom,'start1':a+off+1,'sequence':s,'canonical':canon(s),'gc':round(sum(b in 'GC' for b in s)/k,4),'window_gc':round(sum(b in 'GC' for b in context)/len(context),4)})
    tsv(ROOT/'analysis/probes.tsv',probes)
    # Whole-reference exact scan input deduplicated by canonical sequence.
    with (ROOT/'analysis/scan-probes.tsv').open('w') as f:
        f.write('sequence\n'); f.write('\n'.join(sorted({p['canonical'] for p in probes}))+'\n')
    bam=pysam.AlignmentFile(ROOT/'raw/locus-with-mates.bam','rb')
    rows=list(bam); reads={read_key(x):x for x in rows if primary(x) and x.query_sequence}
    clip=Counter(); selected=set()
    with (ROOT/'raw/all-primary-reads.fa').open('w') as f:
        for key,x in reads.items():
            # SAM sequence is reference-strand oriented; search is strand symmetric.
            f.write('>'+phash(key[0])+f'/{key[1]}\n'+x.query_sequence+'\n')
            if x.reference_name=='chr1' and x.cigartuples and not x.is_unmapped:
                lc=x.cigartuples[0][1] if x.cigartuples[0][0]==4 else 0
                tc=x.cigartuples[-1][1] if x.cigartuples[-1][0]==4 else 0
                if lc>=15:clip[('leading',x.reference_start+1)]+=1
                if tc>=15:clip[('trailing',x.reference_end+1)]+=1
                if (tc>=10 and abs(x.reference_end-L)<=5) or (lc>=10 and abs(x.reference_start-R)<=5):selected.add(key[0])
    with (ROOT/'raw/prior-junction-templates.fa').open('w') as f:
        for key,x in reads.items():
            if key[0] in selected:f.write('>'+phash(key[0])+f'/{key[1]}\n'+x.query_sequence+'\n')
    tsv(ROOT/'analysis/clip-clusters.tsv',[{'side':s,'boundary1':p,'primary_read_ends':n} for (s,p),n in clip.most_common() if n>=3])
    dump(ROOT/'analysis/extraction.json',{'alignments':len(rows),'primary_qc_nonduplicate_read_ends':len(reads),'prior_clip_templates':len(selected),
        'regions_0based_halfopen':REGIONS,'pysam':pysam.__version__})

def eligible_kmers(x,k):
    s=x.query_sequence.upper();q=x.query_qualities
    if q is None:return
    bad=0
    for i in range(len(s)):
        bad+=int(q[i]<20 or s[i] not in 'ACGT')
        if i>=k:bad-=int(q[i-k]<20 or s[i-k] not in 'ACGT')
        if i>=k-1 and bad==0:yield canon(s[i-k+1:i+1])

def assemble():
    bam=pysam.AlignmentFile(ROOT/'raw/locus-with-mates.bam','rb')
    mol=defaultdict(list)
    for x in bam:
        if primary(x) and x.query_sequence and x.reference_name=='chr1' and 15990000<=x.reference_start<16080000:mol[x.query_name].append(x)
    count=Counter()
    for reads in mol.values():count.update(set(y for x in reads for y in eligible_kmers(x,51)))
    edges={k:v for k,v in count.items() if v>=3}
    outgoing=defaultdict(set);incoming=defaultdict(set)
    for k in edges:
        for s in (k,rc(k)):outgoing[s[:-1]].add(s[1:]);incoming[s[1:]].add(s[:-1])
    visited=set();contigs={};metadata=[]
    for u in sorted(outgoing):
        if len(incoming[u])==1 and len(outgoing[u])==1:continue
        for v in sorted(outgoing[u]):
            if (u,v) in visited:continue
            s=u+v[-1];a,b=u,v;support=[edges[canon(s)]]
            while True:
                visited.add((a,b));visited.add((rc(b),rc(a)))
                if len(outgoing[b])!=1 or len(incoming[b])!=1:break
                c=next(iter(outgoing[b]))
                if (b,c) in visited:break
                s+=c[-1];support.append(edges[canon(b+c[-1])]);a,b=b,c
            if len(s)>=100:
                n=f'unitig_{len(contigs)+1}';contigs[n]=s;metadata.append({'contig':n,'length':len(s),'minimum_kmer_templates':min(support),'median_kmer_templates':statistics.median(support)})
    fasta(ROOT/'analysis/unitigs-k51.fa',contigs)
    tsv(ROOT/'analysis/unitigs-k51.tsv',metadata)
    dump(ROOT/'analysis/assembly.json',{'k':51,'minimum_templates':3,'molecules':len(mol),'distinct_canonical_kmers':len(count),'retained_kmers':len(edges),'unitigs_ge100':len(contigs),'longest':max(map(len,contigs.values()),default=0),'note':'Nonbranching de Bruijn unitigs; no molecule-long phase implied; circular components omitted.'})

def count_probes():
    probes=list(csv.DictReader((ROOT/'analysis/probes.tsv').open(),delimiter='\t'))
    scan={r['canonical']:r for r in csv.DictReader((ROOT/'analysis/whole-reference-kmers.tsv').open(),delimiter='\t')}
    targets=set(scan);mol=defaultdict(list)
    for x in pysam.AlignmentFile(ROOT/'raw/locus-with-mates.bam','rb'):
        if primary(x) and x.query_sequence:mol[x.query_name].append(x)
    counts={q:Counter() for q in (0,20,60)};support=defaultdict(set);readrows=[]
    join={p['canonical'] for p in probes if p['kind']=='DEL633_join' and int(scan[p['canonical']]['reference_occurrences'])==0}
    for name,reads in mol.items():
        perq={q:set() for q in counts};nhit=set()
        for x in reads:
            found=set(s for k in (31,51) for s in eligible_kmers(x,k) if s in targets)
            if found&join:
                nhit|=found&join
                readrows.append({'template':phash(name),'read_end':1 if x.is_read1 else 2,'mapq':x.mapping_quality,'reference':x.reference_name,'start1':x.reference_start+1,'end1':x.reference_end,'cigar':x.cigarstring,'reverse':x.is_reverse,'mate_reference':x.next_reference_name,'mate_start1':x.next_reference_start+1,'template_length':x.template_length,'junction_k31':sum(len(s)==31 for s in found&join),'junction_k51':sum(len(s)==51 for s in found&join),'SA':x.get_tag('SA') if x.has_tag('SA') else ''})
            for q in perq:
                if x.mapping_quality>=q:perq[q]|=found
        for q in counts:counts[q].update(perq[q])
        for k in (31,51):
            if any(len(s)==k for s in nhit):support[k].add(name)
    for p in probes:
        z=scan[p['canonical']];p.update({'reference_occurrences':int(z['reference_occurrences']),'reference_positions':z['first_positions']})
        for q in counts:p[f'templates_mapq{q}']=counts[q][p['canonical']]
    tsv(ROOT/'analysis/probe-counts.tsv',probes)
    tsv(ROOT/'analysis/junction-template-support.tsv',readrows)
    dump(ROOT/'analysis/junction-support.json',{'templates_any_novel_k31':len(support[31]),'templates_any_novel_k51':len(support[51]),'qual_threshold_per_base':20,'overlapping_mates_collapsed':True,'primary_qc_nonduplicate_only':True,'minimum_anchor_bases':10})

def loadfasta(p):
    return {x.split('\n',1)[0]:''.join(x.splitlines()[1:]) for x in p.read_text().split('>')[1:]}

def refine():
    hs=loadfasta(ROOT/'analysis/candidate-haplotypes.fa');u=loadfasta(ROOT/'analysis/unitigs-k51.fa')
    left=u['unitig_121'];right=u['unitig_272']
    assert left[-25:]==right[:25]
    assembled=left+right[25:]
    # The two unbranched graph segments define a 79-bp insertion hypothesis;
    # independent raw read spanning of both junctions is evaluated in analyze.
    inserted=assembled[75:154]
    assert len(inserted)==79
    assert assembled[154:]==hs['REF_left'][1000:1124]
    assert sum(a!=b for a,b in zip(assembled[:75],hs['REF_left'][925:1000]))==1
    hs['INS79']=hs['REF_left'][:1000]+inserted+hs['REF_left'][1000:]
    h=list(hs['INS79']);h[939]='C';hs['INS79_upstream_T_C']=''.join(h)
    hs['REF_left_upstream_T_C']=hs['REF_left'][:939]+'C'+hs['REF_left'][940:]
    hs['DEL633_join_upstream_T_C']=hs['DEL633_join'][:939]+'C'+hs['DEL633_join'][940:]
    hs['locus_INS79']=hs['locus_reference'][:L-REGIONS['locus'][1]]+inserted+hs['locus_reference'][L-REGIONS['locus'][1]:]
    hs['locus_DEL9739']=hs['locus_reference'][:L-REGIONS['locus'][1]]+hs['locus_reference'][R-REGIONS['locus'][1]:]
    fasta(ROOT/'analysis/refined-haplotypes.fa',hs)
    fasta(ROOT/'analysis/reconstructed-junction-contigs.fa',{'left_junction_k51':left,'right_junction_k51':right,'merged_INS79_hypothesis':assembled,'inserted79':inserted})
    dump(ROOT/'analysis/refined-allele.json',{'chrom':'chr1','position1':L,'anchor_REF':hs['REF_left'][999],'inserted_sequence':inserted,'insertion_length':len(inserted),'merged_length':len(assembled),'assembly_overlap':25,'upstream_substitution_position1':L-61+1,'upstream_substitution':'T>C','note':'Insertion is a hypothesis until both junctions are observed on one read/template; do not infer deletion from only its left junction.'})
    old=list(csv.DictReader((ROOT/'analysis/probes.tsv').open(),delimiter='\t'))
    new=[]
    for n,boundary in [('INS79_left',1000),('INS79_right',1079)]:
        for k in (31,51):
            for off in range(-(k-10),-9):
                s=hs['INS79'][boundary+off:boundary+off+k]
                new.append({'id':f'{n}_k{k}_o{off}','kind':n,'k':k,'chrom':'chr1','start1':0,'sequence':s,'canonical':canon(s),'gc':round(sum(b in 'GC' for b in s)/k,4),'window_gc':0})
    tsv(ROOT/'analysis/refined-probes.tsv',new)
    with (ROOT/'analysis/refined-scan-probes.tsv').open('w') as f:f.write('sequence\n'+'\n'.join(sorted({r['canonical'] for r in new}))+'\n')

def analyze():
    probes=list(csv.DictReader((ROOT/'analysis/probes.tsv').open(),delimiter='\t'))+list(csv.DictReader((ROOT/'analysis/refined-probes.tsv').open(),delimiter='\t'))
    scan={}
    for fn in ['whole-reference-kmers.tsv','refined-whole-reference-kmers.tsv']:
        scan.update({r['canonical']:r for r in csv.DictReader((ROOT/'analysis'/fn).open(),delimiter='\t')})
    targets=set(scan);mol=defaultdict(list)
    for x in pysam.AlignmentFile(ROOT/'raw/locus-with-mates.bam','rb'):
        if primary(x) and x.query_sequence:mol[x.query_name].append(x)
    h=loadfasta(ROOT/'analysis/refined-haplotypes.fa')
    both=h['INS79'][990:1089]
    # Validate reference flanking 31-mers using exact exhaustive string search,
    # including both orientations and all contigs; only four queries are needed.
    flanks={'left31':h['REF_left'][969:1000],'right31':h['REF_left'][1000:1031],
        'left51':h['REF_left'][949:1000],'right51':h['REF_left'][1000:1051]}
    flankscan={n:{'sequence':s,'positions':[],'count':0} for n,s in flanks.items()}
    ref=pysam.FastaFile(REF)
    priorflanks=ROOT/'analysis/flank-whole-reference-uniqueness.json'
    if priorflanks.exists():
        previous=json.loads(priorflanks.read_text())
        assert all(previous[n]['sequence']==s for n,s in flanks.items())
        flankscan=previous
    else:
        for chrom in ref.references:
            s=ref.fetch(chrom).upper()
            for n,q in flanks.items():
                for strand,needle in [('+',q),('-',rc(q))]:
                    start=0
                    while (i:=s.find(needle,start))>=0:
                        z=flankscan[n];z['count']+=1
                        if len(z['positions'])<10:z['positions'].append([chrom,i+1,strand])
                        start=i+1
    dump(ROOT/'analysis/flank-whole-reference-uniqueness.json',flankscan)
    counts={q:Counter() for q in (0,20,60)};kindsets=defaultdict(set)
    for p in probes:
        if p['kind']!='dosage' and int(scan[p['canonical']]['reference_occurrences'])==0:
            kindsets[(p['kind'],int(p['k']))].add(p['canonical'])
        if p['kind']=='REF_left' and int(scan[p['canonical']]['reference_occurrences'])==1:
            kindsets[('REF_left_unique',int(p['k']))].add(p['canonical'])
    mol_support=defaultdict(set);readrows=[];full=set();full_unique=set();spanning=[]
    for name,reads in mol.items():
        byq={q:set() for q in counts};union=set()
        for x in reads:
            found=set(s for k in (31,51) for s in eligible_kmers(x,k) if s in targets)
            union|=found
            for q in byq:
                if x.mapping_quality>=q:byq[q]|=found
            seq=x.query_sequence.upper();q=x.query_qualities
            exact_both=False;matched_unique=[]
            for strand,needle in [('+',both),('-',rc(both))]:
                at=seq.find(needle)
                if at>=0 and q is not None and min(q[at:at+len(needle)])>=20:
                    exact_both=True
                    for n,z in flankscan.items():
                        for fl in [z['sequence'],rc(z['sequence'])]:
                            idx=seq.find(fl)
                            if idx>=0 and min(q[idx:idx+len(fl)])>=20 and z['count']==1:matched_unique.append(n)
                    full.add(name)
                    if matched_unique:full_unique.add(name)
                    spanning.append({'template':phash(name),'end':1 if x.is_read1 else 2,'strand_in_sam_sequence':strand,'mapq':x.mapping_quality,'start1':x.reference_start+1,'cigar':x.cigarstring,'unique_flanks':','.join(sorted(set(matched_unique))),'minimum_99bp_baseq':min(q[at:at+len(needle)]),'reverse_alignment':x.is_reverse})
            hitkind=[f'{n}_k{k}' for (n,k),ks in kindsets.items() if found&ks]
            if any(z.startswith('INS79') for z in hitkind):
                readrows.append({'template':phash(name),'end':1 if x.is_read1 else 2,'mapq':x.mapping_quality,'reference':x.reference_name,'start1':x.reference_start+1,'end1':x.reference_end,'cigar':x.cigarstring,'reverse':x.is_reverse,'both_junctions_exact99':exact_both,'junction_kmer_types':','.join(hitkind),'mate_reference':x.next_reference_name,'mate_start1':x.next_reference_start+1,'template_length':x.template_length})
        for q in counts:counts[q].update(byq[q])
        for nk,ks in kindsets.items():
            if union&ks:mol_support[nk].add(name)
    for p in probes:
        z=scan[p['canonical']];p.update({'reference_occurrences':int(z['reference_occurrences']),'reference_positions':z['first_positions']})
        for q in counts:p[f'templates_mapq{q}']=counts[q][p['canonical']]
    tsv(ROOT/'analysis/final-probe-counts.tsv',probes)
    tsv(ROOT/'analysis/INS79-template-support.tsv',readrows)
    tsv(ROOT/'analysis/INS79-both-junction-spanning.tsv',spanning)
    out={'primary_templates_examined':len(mol),'both_junctions_exact_99bp_bq20_templates':len(full),'both_junctions_exact_99bp_bq20_read_ends':len(spanning),'full_spanning_with_unique_flank_templates':len(full_unique),'support':{f'{n}_k{k}':len(v) for (n,k),v in mol_support.items()},'left_right_same_template_k31':len(mol_support[('INS79_left',31)]&mol_support[('INS79_right',31)]),'left_right_same_template_k51':len(mol_support[('INS79_left',51)]&mol_support[('INS79_right',51)]),'quality':'Each exact k-mer/99bp match has BQ>=20 at every base; primary nonduplicate QC-pass reads; mates deduplicated per template.'}
    dump(ROOT/'analysis/final-support.json',out)
    dosage=[]
    intervals={'CLCNKA_gene':(16022036,16034050),'CLCNKB_gene':(16043782,16057326),'DEL633_retained_interval':(L+1,R),'CLCNKB_before_insertion':(16043782,L-100),'CLCNKB_after_insertion':(L+100,16057326),'control_left':(15990000,16005000),'control_right':(16065000,16080000),'control_distant':(17000000,17015000)}
    for k in (31,51):
        pp=[p for p in probes if p['kind']=='dosage' and int(p['k'])==k]
        cc=[p for p in pp if int(p['reference_occurrences'])==1 and any(a<=int(p['start1'])<=b-k for n,(a,b) in intervals.items() if n.startswith('control'))]
        for q in (0,20,60):
            controlbygc=defaultdict(list)
            for p in cc:controlbygc[int(float(p['window_gc'])/0.05)].append(p[f'templates_mapq{q}'])
            for n,(a,b) in intervals.items():
                allp=[p for p in pp if a<=int(p['start1'])<=b-k]
                chosen=[p for p in allp if int(p['reference_occurrences'])==1]
                ratios=[]
                for p in chosen:
                    cm=controlbygc[int(float(p['window_gc'])/0.05)]
                    if len(cm)>=10 and statistics.median(cm)>0:ratios.append(p[f'templates_mapq{q}']/statistics.median(cm))
                vals=[p[f'templates_mapq{q}'] for p in chosen]
                dosage.append({'region':n,'start1':a,'end1':b,'k':k,'mapq':q,'probes_total':len(allp),'reference_unique_probes':len(chosen),'median_exact_kmer_templates':statistics.median(vals) if vals else None,'mean_exact_kmer_templates':statistics.fmean(vals) if vals else None,'pooled_control_median':statistics.median(p[f'templates_mapq{q}'] for p in cc),'gc_matched_probes':len(ratios),'median_gc_matched_ratio':statistics.median(ratios) if ratios else None,'mean_gc_matched_ratio':statistics.fmean(ratios) if ratios else None})
    tsv(ROOT/'analysis/paralog-discriminating-dosage.tsv',dosage)
    # Relative local-alignment score comparison, not a calibrated likelihood.
    scores=defaultdict(lambda:defaultdict(dict))
    models={'INS79':'insertion','INS79_upstream_T_C':'insertion','locus_INS79':'insertion','DEL633_join':'deletion','DEL633_join_upstream_T_C':'deletion','locus_DEL9739':'deletion','REF_left':'reference','REF_right':'reference','locus_reference':'reference','REF_left_upstream_T_C':'reference'}
    for line in (ROOT/'analysis/all-read-refined-models.paf').open():
        f=line.split();tags={z[:2]:z[5:] for z in f[12:]};n,end=f[0].rsplit('/',1);model=models[f[5]]
        scores[n][end][model]=max(scores[n][end].get(model,-999),int(tags.get('AS',-999)))
    rows=[]
    interesting={phash(n) for n in mol_support[('INS79_left',31)]|mol_support[('INS79_right',31)]}
    for n in sorted(interesting):
        ends=scores[n]
        complete=bool(ends) and all(all(m in v for m in ('insertion','deletion','reference')) for v in ends.values())
        totals={m:sum(v.get(m,0) for v in ends.values()) for m in ('insertion','deletion','reference')}
        best=sorted(totals,key=totals.get,reverse=True);delta=totals[best[0]]-totals[best[1]]
        rows.append({'template':n,'scored_read_ends':len(ends),'all_models_scored_for_each_end':complete,**totals,'best_model':best[0] if delta else 'tie','best_advantage_AS':delta})
    tsv(ROOT/'analysis/template-model-comparison.tsv',rows)
    dump(ROOT/'analysis/model-comparison-summary.json',{'templates':len(rows),'complete_templates':sum(r['all_models_scored_for_each_end'] for r in rows),'best_model_counts':dict(Counter(r['best_model'] for r in rows)),'advantage_ge20':dict(Counter(r['best_model'] for r in rows if r['best_advantage_AS']>=20 and r['all_models_scored_for_each_end'])),'interpretation':'AS is local alignment score, summed across observed ends, not a probability or formal diplotype fit. Full-locus reference competitor covers both paralogs.'})

def depth():
    import numpy as np
    probes=list(csv.DictReader((ROOT/'analysis/final-probe-counts.tsv').open(),delimiter='\t'))
    mol=defaultdict(list)
    for x in pysam.AlignmentFile(ROOT/'raw/locus-with-mates.bam','rb'):
        if primary(x) and x.query_sequence and not x.is_unmapped and x.reference_name=='chr1':mol[x.query_name].append(x)
    arrays={(n,q):np.zeros(b-a,dtype=np.int32) for n,(c,a,b) in REGIONS.items() for q in (0,20,60)}
    refspan=[]
    for name,reads in mol.items():
        covered={q:set() for q in (0,20,60)}
        for x in reads:
            quals=x.query_qualities
            if quals is None:continue
            aligned={r:q for q,r in x.get_aligned_pairs(matches_only=True) if quals[q]>=20}
            for mq in covered:
                if x.mapping_quality>=mq:covered[mq].update(aligned)
            # Consecutive reference-matching +/-10bp at the target anchor.
            positions=list(range(L-10,L+10))
            if all(r in aligned for r in positions):
                querypos=[aligned[r] for r in positions]
                if querypos==list(range(querypos[0],querypos[0]+20)):
                    refspan.append({'template':phash(name),'end':1 if x.is_read1 else 2,'mapq':x.mapping_quality,'start1':x.reference_start+1,'cigar':x.cigarstring,'observed_20bp':''.join(x.query_sequence[aligned[r]] for r in positions)})
        for n,(c,a,b) in REGIONS.items():
            for q in covered:
                ids=[r-a for r in covered[q] if a<=r<b]
                arrays[(n,q)][ids]+=1
    intervals={'CLCNKA_gene':(16022036,16034050),'CLCNKB_gene':(16043782,16057326),'DEL633_retained_interval':(L+1,R),'CLCNKB_before_insertion':(16043782,L-100),'CLCNKB_after_insertion':(L+100,16057326),'control_left':(15990000,16005000),'control_right':(16065000,16080000),'control_distant':(17000000,17015000)}
    for p in probes:
        if p['kind']!='dosage' or int(p['k'])!=51:continue
        center0=int(p['start1'])-1+25
        for rn,(c,a,b) in REGIONS.items():
            if a<=center0<b:
                for q in (0,20,60):p[f'base_depth_q{q}']=int(arrays[(rn,q)][center0-a])
    pp=[p for p in probes if 'base_depth_q0' in p and int(p['reference_occurrences'])==1]
    controls=[p for p in pp if any(a<=int(p['start1'])<=b-51 for n,(a,b) in intervals.items() if n.startswith('control'))]
    summary=[];windows=[]
    for q in (0,20,60):
        byg=defaultdict(list)
        for p in controls:byg[int(float(p['window_gc'])/.05)].append(p[f'base_depth_q{q}'])
        for n,(a,b) in intervals.items():
            selected=[p for p in pp if a<=int(p['start1'])<=b-51]
            vals=[p[f'base_depth_q{q}'] for p in selected];ratios=[]
            for p in selected:
                cm=byg[int(float(p['window_gc'])/.05)]
                if len(cm)>=10 and statistics.median(cm)>0:ratios.append(p[f'base_depth_q{q}']/statistics.median(cm))
            summary.append({'region':n,'mapq':q,'unique51_centers':len(vals),'mean_base_depth':statistics.fmean(vals),'median_base_depth':statistics.median(vals),'gc_matched_centers':len(ratios),'median_gc_matched_ratio':statistics.median(ratios),'mean_gc_matched_ratio':statistics.fmean(ratios)})
        for a in range(15990000,16080000,1000):
            vv=[p[f'base_depth_q{q}'] for p in pp if a<=int(p['start1'])<a+1000]
            if vv:windows.append({'start1':a,'end1':a+999,'mapq':q,'unique51_centers':len(vv),'median_base_depth':statistics.median(vv),'mean_base_depth':statistics.fmean(vv)})
    tsv(ROOT/'analysis/unique-segment-base-depth.tsv',summary)
    tsv(ROOT/'analysis/unique-segment-depth-windows.tsv',windows)
    tsv(ROOT/'analysis/unique-segment-probes.tsv',pp)
    tsv(ROOT/'analysis/REF-contiguous-spanners.tsv',refspan)
    dump(ROOT/'analysis/reference-spanning-status.json',{'templates_with_contiguous_20bp_target_alignment':len({r['template'] for r in refspan}),'read_ends':len(refspan),'note':'BQ20 all20bp; alignment-contiguous, not necessarily sequence-matching. Zero is not a calibrated genotype probability.'})

if __name__=='__main__':
    os.umask(0o077)
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','assemble','count','refine','analyze','depth']);a=parser.parse_args()
    {'prepare':prepare,'assemble':assemble,'count':count_probes,'refine':refine,'analyze':analyze,'depth':depth}[a.stage]()
