#!/usr/bin/env python3
import argparse, gzip, hashlib, json
from pathlib import Path

CODON = dict(zip(
    [a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],
    list('FFLLSSSSYY**CC*W')+list('LLLLPPPPHHQQRRRR')+
    list('IIIMTTTTNNKKSSRR')+list('VVVVAAAADDEEGGGG')
))

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def load_fai(fai):
    d={}
    with open(fai) as f:
        for line in f:
            n,l,o,b,w=line.rstrip().split('\t')[:5]
            d[n]=(int(l),int(o),int(b),int(w))
    return d

def fetch(fasta, fai, chrom, start0, end0):
    length, offset, bases, width=fai[chrom]
    assert 0 <= start0 <= end0 <= length
    out=[]
    with open(fasta,'rb') as f:
        p=start0
        while p < end0:
            line_i, in_line=divmod(p,bases)
            take=min(end0-p,bases-in_line)
            f.seek(offset+line_i*width+in_line)
            out.append(f.read(take).decode())
            p += take
    return ''.join(out).upper()

def fasta_record(name, seq, width=60):
    return '>'+name+'\n'+'\n'.join(seq[i:i+width] for i in range(0,len(seq),width))+'\n'

def translate(cds):
    assert len(cds)%3 == 0
    return ''.join(CODON[cds[i:i+3]] for i in range(0,len(cds),3))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--outdir', required=True)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    row=None
    with gzip.open(args.annotation,'rt') as f:
        for line in f:
            fields=line.rstrip().split('\t')
            if len(fields)>1 and fields[1]=='NM_001211.6': row=fields; break
    assert row is not None
    _,tx,chrom,strand,txs,txe,cdss,cdse,nex,starts,ends,_,gene,cdss_stat,cdse_stat,frames=row
    assert strand=='+' and gene=='BUB1B' and cdss_stat==cdse_stat=='cmpl'
    txs,txe,cdss,cdse,nex=map(int,(txs,txe,cdss,cdse,nex))
    starts=[int(x) for x in starts.rstrip(',').split(',')]
    ends=[int(x) for x in ends.rstrip(',').split(',')]
    assert len(starts)==len(ends)==nex
    fai=load_fai(args.reference+'.fai')
    exon_seqs=[fetch(args.reference,fai,chrom,s,e) for s,e in zip(starts,ends)]
    mrna=''.join(exon_seqs)
    genome_to_cdna={}
    k=0
    exon_cdna=[]
    for i,(s,e) in enumerate(zip(starts,ends),1):
        exon_cdna.append({'exon':i,'genomic_start_1based':s+1,'genomic_end_1based':e,
                          'cdna_start_1based':k+1,'cdna_end_1based':k+(e-s),'length':e-s})
        for g in range(s,e): genome_to_cdna[g+1]=k+(g-s)+1
        k += e-s
    cds_start_cdna=genome_to_cdna[cdss+1]
    cds_end_cdna=genome_to_cdna[cdse]
    cds=mrna[cds_start_cdna-1:cds_end_cdna]
    assert len(cds)%3==0
    prot_star=translate(cds)
    assert prot_star.endswith('*') and prot_star.count('*')==1
    variants=[
      {'id':'BUB1B_L737Ter','g_pos':40209701,'ref':'T','alt':'G','expected_c':2210,'expected_ref_codon':'TTA','expected_alt_codon':'TGA','expected_p':'p.Leu737Ter'},
      {'id':'BUB1B_N1002K','g_pos':40220612,'ref':'T','alt':'G','expected_c':3006,'expected_ref_codon':'AAT','expected_alt_codon':'AAG','expected_p':'p.Asn1002Lys'},
    ]
    seqs={'reference':mrna}
    details=[]
    for v in variants:
        cdna=genome_to_cdna[v['g_pos']]
        c=cdna-cds_start_cdna+1
        assert c==v['expected_c'] and mrna[cdna-1]==v['ref']
        alt=mrna[:cdna-1]+v['alt']+mrna[cdna:]
        alt_cds=alt[cds_start_cdna-1:cds_end_cdna]
        codon_start=((c-1)//3)*3
        ref_codon=cds[codon_start:codon_start+3]
        alt_codon=alt_cds[codon_start:codon_start+3]
        assert ref_codon==v['expected_ref_codon'] and alt_codon==v['expected_alt_codon']
        exon=next(x for x in exon_cdna if x['cdna_start_1based']<=cdna<=x['cdna_end_1based'])
        pstar=translate(alt_cds)
        first_stop=pstar.index('*')+1
        seqs[v['id']]=alt
        details.append({**v,'cdna_position':cdna,'cds_position':c,'codon_number':(c+2)//3,
                        'ref_codon':ref_codon,'alt_codon':alt_codon,'exon':exon['exon'],
                        'bases_from_ptc_end_to_exon_junction': exon['cdna_end_1based']-(cdna+2) if alt_codon in ('TAA','TAG','TGA') else None,
                        'first_stop_codon_number':first_stop,'translated_residues_before_stop':first_stop-1})
    cis=mrna
    for v in variants:
        cdna=genome_to_cdna[v['g_pos']]
        cis=cis[:cdna-1]+v['alt']+cis[cdna:]
    seqs['cis_both']=cis
    labels={
      'reference':'NM_001211.6_reference_GRCh38_local',
      'BUB1B_L737Ter':'NM_001211.6_c.2210T_G_p.Leu737Ter',
      'BUB1B_N1002K':'NM_001211.6_c.3006T_G_p.Asn1002Lys',
      'cis_both':'NM_001211.6_cis_c.2210T_G_and_c.3006T_G',
    }
    tx_fa=''; cds_fa=''; prot_fa=''
    proteins={}
    for key,s in seqs.items():
        cseq=s[cds_start_cdna-1:cds_end_cdna]
        p=translate(cseq)
        first=p.find('*')
        protein=p[:first] if first>=0 else p
        proteins[key]=protein
        tx_fa += fasta_record(labels[key]+' molecule=mRNA',s)
        cds_fa += fasta_record(labels[key]+' molecule=CDS_including_stop',cseq)
        prot_fa += fasta_record(labels[key]+' molecule=protein translated_to_first_stop',protein)
    # Explicit counterfactual: complete exon-17 skipping. This is not asserted to occur.
    e17=exon_cdna[16]
    skip_ref=mrna[:e17['cdna_start_1based']-1]+mrna[e17['cdna_end_1based']:]
    skip_n1002k=cis[:e17['cdna_start_1based']-1]+cis[e17['cdna_end_1based']:]
    for label,s in [('NM_001211.6_candidate_exon17_skipped',skip_ref),
                    ('NM_001211.6_candidate_exon17_skipped_with_genomic_N1002K',skip_n1002k)]:
        cseq=s[cds_start_cdna-1:cds_end_cdna-e17['length']]
        p=translate(cseq); assert p.endswith('*') and p.count('*')==1
        protein=p[:-1]
        tx_fa += fasta_record(label+' molecule=mRNA counterfactual_not_observed',s)
        cds_fa += fasta_record(label+' molecule=CDS_including_stop counterfactual_not_observed',cseq)
        prot_fa += fasta_record(label+' molecule=protein counterfactual_not_observed',protein)
    assert translate(skip_ref[cds_start_cdna-1:cds_end_cdna-e17['length']])[:-1] == proteins['reference'][:714]+proteins['reference'][761:]
    assert translate(skip_n1002k[cds_start_cdna-1:cds_end_cdna-e17['length']])[:-1][954]=='K'
    assert proteins['reference'].__len__()==1050
    assert proteins['BUB1B_L737Ter']==proteins['cis_both'] and len(proteins['BUB1B_L737Ter'])==736
    assert len(proteins['BUB1B_N1002K'])==1050 and proteins['reference'][1001]=='N' and proteins['BUB1B_N1002K'][1001]=='K'
    (out/'bub1b.transcripts.fa').write_text(tx_fa)
    (out/'bub1b.cds.fa').write_text(cds_fa)
    (out/'bub1b.proteins.fa').write_text(prot_fa)
    meta={
      'assembly':'GRCh38','chromosome':chrom,'strand':strand,'gene':gene,
      'transcript':tx,'protein':'NP_001202.5','mane_select':True,
      'annotation_source':'UCSC hg38 ncbiRefSeqCurated table (local snapshot)',
      'annotation_sha256':sha256(args.annotation),
      'reference_fai_sha256':sha256(args.reference+'.fai'),
      'tx_genomic_0based_halfopen':[txs,txe],'cds_genomic_0based_halfopen':[cdss,cdse],
      'transcript_length':len(mrna),'cds_length_including_stop':len(cds),
      'protein_length':len(proteins['reference']),'exon_count':nex,
      'cds_cdna_1based_inclusive':[cds_start_cdna,cds_end_cdna],
      'exons':exon_cdna,'variants':details,
      'last_exon_junction_after_cdna_base':exon_cdna[-2]['cdna_end_1based'],
      'ptc_to_last_exon_junction_nt':exon_cdna[-2]['cdna_end_1based']-(details[0]['cdna_position']+2),
      'scenarios':{
        'trans':'two hypothetical homologs represented by individual variant FASTA records; L737Ter product 736 aa if RNA escapes NMD, N1002K product 1050 aa',
        'cis':'both variants on one hypothetical homolog; translation terminates at codon 737, so N1002K is not reached on that translated product; second homolog sequence is unknown/reference only as a simplifying counterfactual',
        'candidate_exon17_skip':'counterfactual 141-nt in-frame exon omission, p.(Glu715_Ile761del); not established by this reconstruction',
        'candidate_exon17_skip_with_N1002K':'counterfactual 1003-aa product with the genomic N1002K codon shifted to N955K; not a demonstrated rescue',
      }
    }
    reconstruction=json.dumps(meta,indent=2)+'\n'
    (out/'reconstruction.json').write_text(reconstruction)
    # Compatibility name consumed by the downstream profile/summarization scripts.
    (out/'transcript_reconstruction.json').write_text(reconstruction)
    table=(
      'scenario\tallele_content\tRNA_prediction\tprotein_if_translated\tdomain_consequence\tinterpretation\n'
      'reference\tneither SNV\tnative transcript\t1050 aa\tintact\treference comparator\n'
      'ALT_L737Ter\tc.2210T>G alone\tpositional NMD candidate\t736 aa if escape\tdeletes aa737-1050 including kinase-like/pseudokinase core\tpredicted loss-of-function; RNA/protein unmeasured\n'
      'ALT_N1002K\tc.3006T>G alone\tno PTC; last-exon missense\t1050 aa with N1002K\tlocal C-terminal substitution\tfunctional effect unresolved; no direct activity evidence\n'
      'cis_both\tboth SNVs on one hypothetical homolog\tpositional NMD candidate driven by c.2210T>G\t736 aa if escape; identical to ALT_L737Ter protein\tN1002K codon is downstream and not translated\tcis is hypothetical; other homolog remains unspecified\n'
      'trans_pair\tone SNV per hypothetical homolog\tone NMD-candidate transcript plus one full-length missense transcript\t736 aa if stop transcript escapes, plus 1050 aa N1002K\tcomplementary truncating+missense hypothesis\ttrans is hypothetical; pathogenicity of N1002K unresolved\n'
      'candidate_exon17_skip\tcomplete exon17 omission; counterfactual\t3528-nt mRNA; in-frame 141-nt omission\t1003 aa, p.(Glu715_Ile761del)\t47-aa in-frame deletion; stop-variant site omitted\tnot observed or model-supported by this reconstruction\n'
      'candidate_exon17_skip_with_N1002K\texon17 omitted plus genomic N1002K allele\t3528-nt mRNA\t1003 aa, p.(Glu715_Ile761del) with N1002K becoming N955K in skipped product\t47-aa deletion plus C-terminal missense\tcounterfactual only; do not call rescue\n')
    (out/'consequences.tsv').write_text(table)
    report=f'''# BUB1B — rekonstrukcja transkryptu, CDS, białka i scenariuszy cis/trans

**Zakres:** lokalna, deterministyczna rekonstrukcja dwóch podanych SNV na GRCh38. Nie wykonano fazowania, analizy RNA ani funkcjonalnej. Scenariusze cis/trans są jawnie kontrfaktyczne i nie rozstrzygają rzeczywistej fazy.

## Wersja odniesienia i kontrola mapowania

- MANE Select: **NM_001211.6 / NP_001202.5**, BUB1B, nić `+`.
- Lokalna migawka UCSC `hg38/ncbiRefSeqCurated`: transkrypt chr15:`40161068-40221123`, CDS `40161220-40220759` (obie pary 0-based, half-open), 23 eksony.
- Złożony lokalnie mRNA ma {len(mrna)} nt, CDS z kodonem STOP {len(cds)} nt, białko referencyjne {len(proteins['reference'])} aa.
- FASTA odniesienia: NCBI GCA_000001405.15 GRCh38 no-alt analysis set, lokalny SHA-256 `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`; FAI SHA-256 `{meta['reference_fai_sha256']}`.
- Adnotacja UCSC lokalny SHA-256 `{meta['annotation_sha256']}`. Dokładne granice eksonów i współrzędne cDNA są w `reconstruction.json`.

## Odtworzone konsekwencje

| SNV GRCh38 | MANE HGVS | Ekson | Kodon | Produkt |
|---|---|---:|---|---|
| chr15:40209701:T>G | NM_001211.6:c.2210T>G; NP_001202.5:p.Leu737Ter | 17/23 | TTA→TGA | PTC w kodonie 737; 736 aa, jeśli transkrypt uniknie NMD |
| chr15:40220612:T>G | NM_001211.6:c.3006T>G; NP_001202.5:p.Asn1002Lys | 23/23 | AAT→AAG | pełna długość 1050 aa z N1002K |

Obserwowane **T>G** przy c.3006 jest odrębnym allelem DNA od publicznego c.3006T>A, mimo że oba dają p.Asn1002Lys. Nie przenoszono klasyfikacji między allelami.

## NMD — ocena warunkowa z pozycji

PTC c.2210 leży w eksonie 17. Koniec kodonu STOP jest **72 nt** przed złączem ekson 17→18 i **745 nt** przed ostatnim złączem ekson 22→23. Spełnia więc klasyczną przesłankę pozycyjną EJC/NMD (PTC co najmniej około 50–55 nt przed dalszym złączem), ale reguła ma wyjątki; bez RNA nie wiadomo, czy i w jakim stopniu ten konkretny transkrypt ulega NMD. [Kurosaki et al., Nat Rev Mol Cell Biol 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6855384/).

## Scenariusze cis/trans

- **Trans (hipotetyczny):** jeden homolog niesie L737Ter, drugi N1002K. Pierwszy transkrypt jest kandydatem do NMD; jeśli ucieknie, daje 736-aa produkt. Drugi daje pełne 1050-aa białko N1002K. Jest to model „truncating + missense”, nie dowód biallelicznego mechanizmu, bo faza i szkodliwość N1002K pozostają nierozstrzygnięte.
- **Cis (hipotetyczny):** oba SNV są na jednym homologu. Translacja zatrzymuje się w kodonie 737, więc kodon 1002 nie jest osiągany; białkowe FASTA `cis_both` jest identyczne z `ALT_L737Ter`. Drugi homolog jest nieznany; wariant referencyjny drugiego homologu można traktować tylko jako uproszczony kontrfaktyczny komparator, nie obserwację.

### Osobny kontrfaktyczny transkrypt z pominięciem eksonu 17

Pełne pominięcie eksonu 17 (141 nt) byłoby w ramce i usuwałoby 47 aa: **p.(Glu715_Ile761del)**. Usunęłoby też nukleotyd SNV L737Ter. Jeśli na tym samym hipotetycznym allelu byłby genomowy N1002K, kodon pozostałby osiągalny, lecz po delecji miałby pozycję **N955K** w 1003-aa produkcie. Ten kontrfaktyczny scenariusz wybrano wyłącznie dlatego, że wcześniejszy słaby scorer złącza wskazał współrzędne 40208770→40210109 (zgłoszone maksimum około 0,154); sama współrzędna ani taka wielkość nie dowodzi obecności izoformu. FASTA są oznaczone `candidate_exon17_skipped`; niniejsza rekonstrukcja nie nazywa wyniku „rescue”. Złącze to w zapisie 0-based ekson 16→18; kanoniczne złącza flankujące to 40208770→40209634 oraz 40209775→40210109.

## Konsekwencja domenowa — ostrożna interpretacja

L737Ter usuwa aa 737–1050 i tym samym cały lokalnie adnotowany rdzeń kinase-like/pseudokinase (InterPro w poprzednim przeglądzie: aa 791–984) oraz C-koniec zawierający N1002. Badania MVA wykazały dwie klasy defektów BUBR1: obniżenie obfitości białka oraz swoiste defekty checkpointu/przyłączenia mikrotubul; nie wolno zatem sprowadzać każdej zmiany do „utraty aktywności kinazy” ([Suijkerbuijk et al., Cancer Res 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2887387/)). Dla BUBR1 kataliza nie jest wymagana do poprawnej segregacji, natomiast fold pseudokinazowy wspiera stabilność ([Suijkerbuijk et al., Dev Cell 2012](https://pubmed.ncbi.nlm.nih.gov/22698286/)). Usunięcie/niestabilność C-końcowej pseudokinazy może też ograniczać fosforylację KARD i rekrutację PP2A-B56 ([Cordeiro et al., Cell Reports 2020](https://www.sciencedirect.com/science/article/pii/S2211124720313863)). Są to mechanizmy klasowe; **nie ma tu bezpośredniego testu N1002K**, jego aktywności, stabilności ani interakcji.

Truncacja po aa736 zachowuje sekwencję N-końcową 1–736, w tym regiony checkpointu/KEN i KARD; sama obecność tych fragmentów w przewidywanym polipeptydzie nie dowodzi ich ekspresji, lokalizacji ani funkcji. Jeżeli zadziała NMD, skrócony produkt może być bardzo ograniczony lub nieobecny.

## Pliki

- `bub1b.transcripts.fa`: referencja, każdy ALT osobno, cis-both oraz dwa jawnie kontrfaktyczne warianty exon17-skipped, mRNA.
- `bub1b.cds.fa`: odpowiadające CDS, z pełnym referencyjnym zakresem CDS i kodonami STOP.
- `bub1b.proteins.fa`: translacja do pierwszego kodonu STOP; cztery kanoniczne scenariusze 1050/736/1050/736 aa oraz dwa kontrfaktyczne produkty exon17-skipped po 1003 aa.
- `consequences.tsv`: skrót skutków i warunków.
- `reconstruction.json`: eksony, mapowanie i kontrole maszynowe.

**Wniosek:** rekonstrukcja potwierdza oba HGVS na MANE i pokazuje asymetrię scenariuszy. W trans N1002K może stanowić oddzielny pełnodługościowy allel, lecz jego szkodliwość nie jest dowiedziona. W cis L737Ter maskuje N1002K na poziomie translacji tego samego transkryptu. Rzeczywista faza pozostaje nierozstrzygnięta.
'''
    (out/'report.md').write_text(report)
    data_files=['bub1b.transcripts.fa','bub1b.cds.fa','bub1b.proteins.fa','consequences.tsv','reconstruction.json','transcript_reconstruction.json']
    summary={
      'status':'complete','scope':'local deterministic MANE transcript/CDS/protein reconstruction',
      'assembly':'GRCh38','transcript':'NM_001211.6','protein':'NP_001202.5','strand':'+',
      'variants':[
        {'variant':'chr15:40209701:T>G','hgvs_c':'c.2210T>G','hgvs_p':'p.Leu737Ter','exon':'17/23','protein_length_if_escape':736,'nmd':'positional candidate; RNA unmeasured'},
        {'variant':'chr15:40220612:T>G','hgvs_c':'c.3006T>G','hgvs_p':'p.Asn1002Lys','exon':'23/23','protein_length':1050,'effect':'unresolved'},
      ],
      'cis_trans':'unresolved; cis and trans FASTA/scenarios are hypothetical',
      'cis_key_result':'L737Ter stops translation before residue 1002; cis protein equals L737Ter-only protein',
      'trans_key_result':'separate 736-aa-if-escape truncation and full-length 1050-aa N1002K products',
      'candidate_exon17_skip':{'status':'counterfactual_not_observed','selection_reason':'coordinate matched a previously noted weak junction scorer (reported max ~0.154); not proof of isoform','junction_0based':'40208770->40210109','deleted_nt':141,'protein':'p.(Glu715_Ile761del)','length':1003,'with_genomic_N1002K':'N955K in skipped product'},
      'input_hashes':{'reference_fasta': '9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7',
                      'reference_fai':meta['reference_fai_sha256'],'annotation':meta['annotation_sha256']},
      'output_hashes':{name:sha256(out/name) for name in data_files},
      'limitations':['no phase evidence','no RNA/protein measurement','no direct functional assay of N1002K','scenario sequences are counterfactual, not inferred haplotypes']
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    for p in out.iterdir(): p.chmod(0o600)

if __name__=='__main__': main()
