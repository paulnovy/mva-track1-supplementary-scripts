#!/usr/bin/env python3
"""Bounded, allele-exact MVA1 annotation; never modifies historical sources.

The Ensembl skill wrapper performs VEP queries and prediction rendering. We add
documented transcript metadata flags that its CLI does not expose. Human GRCh38
is explicit; local FASTA/VCF checking precedes all public annotation queries.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import time

import pysam

ROOT = Path(os.environ.get('MVA_DATA_ROOT', '.')).resolve()
RUN = 'candidate-annotation-20260911'
RES, REP, OPS = (ROOT / kind / RUN for kind in ('results', 'reports', 'operations'))
API_DIR = Path(__file__).resolve().parents[1] / 'public_api'
WRAPPER = Path(os.environ.get('ENSEMBL_VEP_WRAPPER', API_DIR / 'ensembl_api.py'))
VEP_FLAGS = '&hgvs=1&mane=1&canonical=1&transcript_version=1&protein=1&numbers=1&variant_class=1'
FASTA = ROOT / 'references/GRCh38_standard/GRCh38_standard.fa'
VCF = ROOT / 'results/vcf-normalized-20260907/normalized.vcf.gz'
CLINVAR = ROOT / 'resources/clinvar/prepared-grch38/clinvar.normalized.vcf.gz'
EXROOT = ROOT / 'results/phenotype-sv-mei-followup-20260910/exomiser'
GENES = dict.fromkeys(('BUB1B', 'FANCD2', 'LZTR1', 'ATP6V1E1', 'GNRHR', 'PIK3C2A'),
                      'Selected target gene; no case interpretation encoded')
ABSENT_GENES = ['BUB1', 'CLCNKA', 'CLCNKB', 'SLC34A1', 'SLC34A3', 'CYP24A1', 'CLDN16', 'CLDN19', 'SLC12A1', 'CASR']
MANIFEST = {}
MANIFEST_PATH = None


def configure_manifest(path):
    """Read private inputs; no observations or inherited case conclusions ship here."""
    global MANIFEST, MANIFEST_PATH, GENES, ABSENT_GENES, FASTA, VCF, CLINVAR, EXROOT
    MANIFEST_PATH = path.resolve()
    MANIFEST = json.loads(path.read_text())
    if not isinstance(MANIFEST, dict):
        raise ValueError('manifest must be a JSON object')
    GENES = dict.fromkeys(MANIFEST.get('genes', list(GENES)), 'User-selected target gene')
    ABSENT_GENES = MANIFEST.get('audit_genes', ABSENT_GENES)
    paths = MANIFEST.get('input_paths', {})
    FASTA = private_path(paths['reference_fasta']) if 'reference_fasta' in paths else FASTA
    VCF = private_path(paths['normalized_vcf']) if 'normalized_vcf' in paths else VCF
    CLINVAR = private_path(paths['clinvar_vcf']) if 'clinvar_vcf' in paths else CLINVAR
    EXROOT = private_path(paths['exomiser_root']) if 'exomiser_root' in paths else EXROOT


def private_path(value):
    path = Path(value)
    return path if path.is_absolute() else MANIFEST_PATH.parent / path


def inherited_annotations():
    """Preserve explicitly supplied interpretation without presenting it as computed."""
    rows = []
    for annotation in MANIFEST.get('inheritance_annotations', []):
        row = dict(annotation)
        if 'gene' not in row:
            raise ValueError('each inheritance annotation requires gene')
        row['annotation_source'] = 'user_supplied_private_manifest; not inferred by this script'
        row['phase_rule_interpretation'] = phase_interpretation(
            row.get('phase', 'unresolved'), row.get('pathogenic_alleles_confirmed') is True)
        rows.append(row)
    return rows


def save(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_text(json.dumps(obj, indent=2, default=str) + '\n')
    path.chmod(0o600)


def tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(rows[0])
    with path.open('w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list, tuple)) else v for k, v in row.items()})
    path.chmod(0o600)


def fingerprint(path):
    s = path.stat()
    return {'path': str(path), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns,
            'inode': s.st_ino, 'owner_uid': s.st_uid, 'mode': oct(s.st_mode & 0o777)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def percent_to_fraction(value):
    """Exomiser 15.1.0 Frequency stores percentage, missing is not zero."""
    return None if value in (None, '') else str(Decimal(str(value)) / Decimal(100))


def phase_interpretation(phase, pathogenic_alleles_confirmed=False):
    if phase == 'cis_supported':
        return 'nominated_pair_not_biallelic; independent_single_allele_effects_retained'
    if phase == 'trans_supported' and pathogenic_alleles_confirmed:
        return 'biallelic_genotype_supported; clinical_fit_requires_separate_assessment'
    return 'biallelic_causality_not_established'


def minimal_alleles(pos, ref, alt):
    """Unpad VCF alleles for the VEP region endpoint, without shifting repeats."""
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos+1, ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref, alt


def select():
    rows_by_profile = {}
    selected = {}
    for profile in ('full', 'axis'):
        path = EXROOT / profile / f'work/rare-disease/native/MVA1-{profile}-exomiser.variants.tsv'
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        rows_by_profile[profile] = rows
        for r in rows:
            if r['GENE_SYMBOL'] not in GENES:
                continue
            key = '-'.join(r[x] for x in ('CONTIG', 'START', 'REF', 'ALT'))
            c = selected.setdefault(key, {'id': key, 'gene': r['GENE_SYMBOL'], 'chrom': 'chr' + r['CONTIG'],
                'pos': int(r['START']), 'ref': r['REF'], 'alt': r['ALT'], 'selection_reason': GENES[r['GENE_SYMBOL']], 'exomiser_rows': []})
            c['exomiser_rows'].append({'profile': profile, 'source': str(path), **r})
    for variant in MANIFEST.get('extra_variants', []):
        gene, chrom, pos, ref, alt = (variant[k] for k in ('gene', 'chrom', 'pos', 'ref', 'alt'))
        chrom = str(chrom).removeprefix('chr')
        pos = int(pos)
        reason = 'Additional target from private manifest'
        key = f'{chrom}-{pos}-{ref}-{alt}'
        selected.setdefault(key, dict(id=key, gene=gene, chrom='chr'+chrom, pos=pos, ref=ref, alt=alt, selection_reason=reason, exomiser_rows=[]))
    exclusions = []
    for gene in ABSENT_GENES:
        counts = {p: sum(r['GENE_SYMBOL'] == gene for r in rows) for p, rows in rows_by_profile.items()}
        exclusions.append(dict(gene=gene, **counts, status='no_candidate_in_current_ranked_small_variant_outputs' if not any(counts.values()) else 'present_in_ranked_small_variant_outputs',
            interpretation='Not a gene-wide absence claim; BUB1 is separate from BUB1B; structural CLCNKA/B task is separate.'))
    reference = pysam.FastaFile(str(FASTA))
    source = pysam.VariantFile(str(VCF))
    clinvar = pysam.VariantFile(str(CLINVAR))
    before = [fingerprint(p) for p in (FASTA, VCF, CLINVAR)]
    for c in selected.values():
        chrom, pos, ref, alt = (c[k] for k in ('chrom','pos','ref','alt'))
        observed_ref = reference.fetch(chrom, pos-1, pos-1+len(ref))
        assert observed_ref == ref, (c['id'], observed_ref, ref)
        found = [r for r in source.fetch(chrom, pos-1, pos) if r.pos == pos and r.ref == ref and r.alts == (alt,)]
        assert len(found) == 1, (c['id'], len(found))
        r = found[0]
        call = next(iter(r.samples.values()))
        c['source_call'] = dict(gt=call.get('GT'), phased=call.phased, dp=call.get('DP'), gq=call.get('GQ'), ad=call.get('AD'), filter=list(r.filter), qual=r.qual)
        c['reference_check'] = 'matched_exact_GRCh38_FASTA'
        matches = [r for r in clinvar.fetch(chrom, pos-1, pos) if r.pos == pos and r.ref == ref and r.alts == (alt,)]
        c['clinvar_snapshot'] = [dict(variation_id=r.id, chrom=r.chrom, pos=r.pos, ref=r.ref, alt=r.alts[0], info=dict(r.info)) for r in matches]
        c['clinvar_snapshot_date'] = MANIFEST.get('clinvar_snapshot_date')
        c['clinvar_assertion_evaluation_date'] = None  # VCF snapshot does not contain individual assertion dates.
        c['clinvar_match_status'] = 'exact_normalized_allele_match' if matches else 'no_exact_record_in_snapshot'
    if not selected:
        raise ValueError('No variants selected from input files and private manifest')
    save(RES / 'selection.json', list(selected.values()))
    save(OPS / 'source_before.json', before)
    tsv(RES / 'selection.tsv', [{k:v for k,v in c.items() if k not in ('exomiser_rows', 'clinvar_snapshot')} for c in selected.values()])
    tsv(RES / 'shortlist_exclusions.tsv', exclusions)
    save(RES / 'clinvar_exact_matches.json', {c['id']:c['clinvar_snapshot'] for c in selected.values()})
    save(RES / 'clinvar_header.json', {'header':str(clinvar.header)})


def query(archive=False):
    spec = importlib.util.spec_from_file_location('ensembl_skill_wrapper', WRAPPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    endpoint = 'https://sep2025.rest.ensembl.org' if archive else mod.BASE_URL
    mod._CLIENT_REGULAR = mod.http_client.HttpClient(endpoint, qps=2, timeout=45, max_retries=1,
        user_agent='MVA1-internal-annotation/1.0', default_headers={'Content-Type':'application/json'})
    mod.VEP_PLUGINS += VEP_FLAGS
    commands = []
    for c in json.loads((RES/'selection.json').read_text()):
        out = RES / ('vep_release115' if archive else 'vep') / (c['id'] + '.json')
        out.parent.mkdir(mode=0o700, exist_ok=True)
        arguments = argparse.Namespace(variant_str=c['id'].replace('-', ':'), species='human', assembly='GRCh38', output=str(out))
        # Split the ID manually: the IDs in this bounded selection contain no symbolic alleles.
        pos, ref, alt = minimal_alleles(c['pos'], c['ref'], c['alt'])
        assert ref  # This bounded shortlist contains only substitutions/deletions.
        arguments.variant_str = f"{c['chrom'][3:]}:{pos}:{ref}:{alt or '-'}"
        command = {'id':c['id'], 'endpoint':endpoint, 'call':'provided ensembl_api.py cmd_vep', 'arguments':vars(arguments),
                   'additional_documented_flags':VEP_FLAGS, 'wrapper_sha256':sha(WRAPPER)}
        if out.exists():
            command['status'] = 'existing_completed_response_reused'
        else:
            buffer = io.StringIO()
            start = time.time()
            try:
                with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
                    mod.cmd_vep(arguments)
                command['status'] = 'completed'
            except Exception as e:
                command['status'] = 'failed'
                command['error'] = repr(e)
            command['elapsed_seconds'] = time.time()-start
            out.with_suffix('.stdout.txt').write_text(buffer.getvalue())
            print(c['id'], command['status'], flush=True)
        commands.append(command)
        save(OPS/('vep_release115_commands.json' if archive else 'vep_commands.json'), commands)
        time.sleep(0.15)
    for name, route in [('data_release','/info/data'),('rest_version','/info/rest'),('software_version','/info/software')]:
        try:
            data = mod._get_client('GRCh38').fetch_json(route)
        except Exception as exc:
            data = {'status':'unavailable', 'error':repr(exc)}
        save(RES / ('vep_release115' if archive else 'vep') / f'ensembl_{name}.json', data)


def gnomad(qc=False):
    wrapper = Path(os.environ.get('GNOMAD_WRAPPER', API_DIR / 'get_variant_frequency.py'))
    spec = importlib.util.spec_from_file_location('gnomad_skill_wrapper', wrapper)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if qc:
        # The supplied wrapper does not request filters. Add these public
        # fields at its HTTP boundary; still execute its query/rate-limit code.
        original_fetch = mod.CLIENT.fetch_json
        def with_filters(url, **kwargs):
            body = kwargs['json_body']
            body['query'] = body['query'].replace('      rsids\n', '      rsids\n      flags\n').replace('      exome {\n', '      exome {\n        filters\n').replace('      genome {\n', '      genome {\n        filters\n')
            return original_fetch(url, **kwargs)
        mod.CLIENT.fetch_json = with_filters
    commands = []
    for c in json.loads((RES/'selection.json').read_text()):
        out = RES/('gnomad_qc' if qc else 'gnomad')/(c['id']+'.json')
        out.parent.mkdir(mode=0o700, exist_ok=True)
        command = {'id':c['id'], 'wrapper':str(wrapper), 'wrapper_sha256':sha(wrapper), 'dataset':'gnomad_r4', 'extra_fields':['variant.flags','exome.filters','genome.filters'] if qc else []}
        if out.exists():
            command['status'] = 'existing_completed_response_reused'
        else:
            try:
                mod.get_variant_frequency(c['id'], None, 'gnomad_r4', str(out))
                command['status'] = 'response_saved'
            except Exception as exc:
                command.update(status='failed', error=repr(exc))
            print(c['id'], command['status'], flush=True)
        commands.append(command)
        save(OPS/('gnomad_qc_commands.json' if qc else 'gnomad_commands.json'), commands)
        time.sleep(6.1)  # honour 10 requests/minute even across resume boundaries


def assemble():
    candidates = json.loads((RES/'selection.json').read_text())
    frequencies, transcripts, summaries = [], [], []
    appendix = ['# Complete MANE/canonical VEP prediction rows', '',
        'Human GRCh38; official Ensembl September 2025 REST archive, release 115. '
        'The supplied ensembl-database wrapper renders these rows. Other transcripts '
        'and all unprinted response fields remain in the complete JSON and stdout files.', '']
    for c in candidates:
        cid = c['id']
        raw = RES/'vep_release115'/(cid+'.json')
        vep = json.loads(raw.read_text()) if raw.exists() else []
        primary = []
        if vep:
            assert len(vep) == 1
            top = vep[0]
            assert top['assembly_name'] == 'GRCh38'
            assert top['seq_region_name'] == c['chrom'][3:]
            vr, va = top['allele_string'].split('/')
            vkey = minimal_alleles(top['start'], vr.replace('-',''), va.replace('-',''))
            assert vkey == minimal_alleles(c['pos'],c['ref'],c['alt']), (cid, vkey)
            tc = [t for t in top.get('transcript_consequences',[]) if t.get('gene_symbol') == c['gene']]
            primary = [t for t in tc if t.get('mane_select')] or [t for t in tc if t.get('canonical')]
            assert primary, (cid, 'No nominated gene primary transcript')
            for t in primary:
                transcripts.append({'id':cid,'gene':c['gene'],'selection':'MANE_Select' if t.get('mane_select') else 'canonical_fallback',**t})
            ids = {t['transcript_id'] for t in primary}
            printed = raw.with_suffix('.stdout.txt').read_text().splitlines()
            appendix += [f'## {c["gene"]} {c["chrom"]}:{c["pos"]}:{c["ref"]}>{c["alt"]}', '',
                '| Transcript ID | Gene | Method/Metric | Value |', '| --- | --- | --- | --- |']
            appendix += [line for line in printed if line.startswith('| ') and line.split('|')[1].strip() in ids]
            appendix += ['', 'Additional primary-transcript fields (including HGVS and any LoF flags/information):', '', '```json', json.dumps(primary, indent=2), '```', '']
        c['vep_status'] = 'completed_exact_allele_checked' if primary else 'unavailable'
        c['vep_primary'] = primary
        c['vep_json'] = str(raw) if raw.exists() else None
        p0 = primary[0] if primary else {}
        cv = c['clinvar_snapshot']
        summary = {'id':cid,'gene':c['gene'],'source_GT':'/'.join(map(str,c['source_call']['gt'])),
            'vep_status':c['vep_status'],'primary_transcript':p0.get('transcript_id'),
            'MANE_RefSeq':p0.get('mane_select'),'HGVS_c':p0.get('hgvsc'),'HGVS_p':p0.get('hgvsp'),
            'consequence':p0.get('consequence_terms'), 'loftee':p0.get('lof'),
            'clinvar_exact_status':c['clinvar_match_status'],'clinvar_variation_ids':[r['variation_id'] for r in cv],
            'clinvar_significance':[r['info'].get('CLNSIG') for r in cv],
            'clinvar_review_status':[r['info'].get('CLNREVSTAT') for r in cv],
            'clinvar_snapshot':c.get('clinvar_snapshot_date'),'clinvar_individual_assertion_evaluation_date':None,
            'exomiser_HGVS':sorted({r['HGVS'] for r in c['exomiser_rows']}),
            'selection_reason':c['selection_reason']}
        summaries.append(summary)
        exrows = c['exomiser_rows']
        if exrows:
            ex = exrows[0]
            for r in exrows:
                assert r['ALL_FREQ'] == ex['ALL_FREQ']
            for item in ex['ALL_FREQ'].split(',') if ex['ALL_FREQ'] else ['']:
                source, percent = item.split('=',1) if item else ('no_frequency_returned','')
                frequencies.append(dict(id=cid,source='Exomiser_15.1.0_data2512',subset=source,
                    raw_value=percent or None,raw_units='percent',af_fraction=percent_to_fraction(percent),
                    ac=None,an=None,homozygotes=None,qc='not_exposed_in_TSV',flags=None,
                    coverage='not_exposed_in_TSV',source_path=ex['source']))
        gp=RES/'gnomad_qc'/(cid+'.json')
        gj=json.loads(gp.read_text())
        assert not gj.get('errors'), (cid,gj.get('errors'))
        gv=gj['data']['variant']
        assert gv['variant_id'] == cid
        c['gnomad_flags']=gv.get('flags')
        for subset in ('exome','genome','joint'):
            d=gv.get(subset)
            if d is None:
                frequencies.append(dict(id=cid,source='gnomAD_r4_API',subset=subset,raw_value=None,
                    raw_units='fraction',af_fraction=None,ac=None,an=None,homozygotes=None,
                    qc='not_returned',flags=gv.get('flags'),coverage='not_returned; not_zero_AF',source_path=str(gp)))
                continue
            af=d['ac']/d['an'] if d['an'] else None
            if d.get('af') is not None:
                assert abs(d['af']-af) < 1e-12
            filters=d.get('filters')
            qc='joint_filter_field_not_requested' if subset=='joint' else ('PASS_no_filters_returned' if filters==[] else json.dumps(filters))
            frequencies.append(dict(id=cid,source='gnomAD_r4_API',subset=subset,raw_value=d.get('af',af),raw_units='fraction',
                af_fraction=af,ac=d['ac'],an=d['an'],homozygotes=d['homozygote_count'],qc=qc,flags=gv.get('flags'),
                coverage='AN is genotype callability; no depth track queried',source_path=str(gp)))
            for pop in d.get('populations',[]):
                if '_' in pop['id']:  # Preserve sex strata in raw JSON; avoid duplicated presentation.
                    continue
                frequencies.append(dict(id=cid,source='gnomAD_r4_API',subset=subset+'_'+pop['id'],raw_value=None,
                    raw_units='counts',af_fraction=pop['ac']/pop['an'] if pop['an'] else None,
                    ac=pop['ac'],an=pop['an'],homozygotes=pop['homozygote_count'],qc=qc,flags=gv.get('flags'),
                    coverage='AN is genotype callability; no depth track queried',source_path=str(gp)))
    tsv(RES/'candidate_summary.tsv',summaries)
    save(RES/'candidate_audit.json',candidates)
    tsv(RES/'frequency_audit.tsv',frequencies)
    save(RES/'primary_transcript_annotations.json',transcripts)
    allkeys=sorted(set().union(*(set(t) for t in transcripts)))
    tsv(RES/'primary_transcript_annotations.tsv',transcripts,allkeys)
    (REP/'VEP_PRIMARY_PREDICTIONS.md').write_text('\n'.join(appendix))
    # Supplied prior evidence is copied, never recalculated or treated as verified here.
    prior_phase = {}
    for label, value in MANIFEST.get('phase_evidence', {}).items():
        path = private_path(value)
        prior_phase[label] = {'source':str(path), 'sha256':sha(path),
            'evidence_status':'user_supplied_prior_evidence; not independently verified',
            'data':json.loads(path.read_text())}
    save(RES/'prior_phase_evidence.json', prior_phase)
    inherited = inherited_annotations()
    fields = sorted(set().union(*(set(row) for row in inherited))) if inherited else ['gene','annotation_source','phase_rule_interpretation']
    tsv(RES/'inheritance_audit.tsv', inherited, fields)
    save(RES/'inheritance_audit.json', inherited)
    save(OPS/'assembly_summary.json',dict(selected=len(candidates),vep_completed=sum(c['vep_status'].startswith('completed') for c in candidates),
        clinvar_exact_matches=sum(bool(c['clinvar_snapshot']) for c in candidates),frequency_rows=len(frequencies),primary_transcript_rows=len(transcripts),
        source_metadata_unchanged=json.loads((OPS/'source_before.json').read_text())==[fingerprint(p) for p in (FASTA,VCF,CLINVAR)]))


def verify():
    summary = json.loads((OPS/'assembly_summary.json').read_text())
    candidates = json.loads((RES/'candidate_audit.json').read_text())
    selected = len(candidates)
    assert selected > 0 and summary['selected'] == selected
    assert summary['vep_completed'] == selected
    assert all(c.get('vep_primary') for c in candidates)
    assert summary['primary_transcript_rows'] == sum(len(c['vep_primary']) for c in candidates)
    assert summary['clinvar_exact_matches'] == sum(bool(c['clinvar_snapshot']) for c in candidates)
    assert summary['source_metadata_unchanged']
    releases = json.loads((RES/'vep_release115/ensembl_data_release.json').read_text())['releases']
    assert releases == [115]
    assert json.loads((OPS/'source_before.json').read_text()) == [fingerprint(p) for p in (FASTA,VCF,CLINVAR)]
    for c in candidates:
        assert c['reference_check'] == 'matched_exact_GRCh38_FASTA'
        response = json.loads((RES/'gnomad_qc'/(c['id']+'.json')).read_text())
        assert not response.get('errors')
        variant = response['data']['variant']
        assert variant['variant_id'] == c['id']
        for subset in ('exome', 'genome', 'joint'):
            values = variant.get(subset)
            if values and values.get('af') is not None:
                assert values['an'] > 0
                assert abs(values['af'] - values['ac']/values['an']) < 1e-12
    sources = {}
    for label, value in MANIFEST.get('supplemental_sources', {}).items():
        path = private_path(value)
        sources[label] = {'path':str(path), 'sha256':sha(path)}
    provenance = {
        'date_utc':datetime.now(timezone.utc).isoformat(), 'status':'bounded_analysis_outputs_verified',
        'species':'Homo sapiens', 'assembly':'GRCh38', 'selected_variants':selected,
        'ensembl_primary_releases':releases,
        'clinvar_snapshot_date':MANIFEST.get('clinvar_snapshot_date'),
        'gnomad_dataset':'gnomad_r4 API; exact allele-specific raw counts and QC preserved',
        'large_input_hash_policy':'No inherited hash claimed; metadata and selected REF checks only.',
        'private_manifest':{'path':str(MANIFEST_PATH), 'sha256':sha(MANIFEST_PATH)},
        'sources':sources, 'pysam_version':pysam.__version__, 'code_sha256':sha(Path(__file__)),
        'source_metadata_unchanged':True,
        'limitations':[
            'Clinical and inheritance annotations supplied in the private manifest are not independently verified.',
            'No new genome-wide callset, clinical classification, numeric reranking or phase assignment.',
            'ClinVar snapshot and available API annotations may be incomplete; unavailable fields are not negative evidence.',
            'Population depth coverage is not queried; AN and filter fields are reported instead.']}
    save(OPS/'provenance.json', provenance)
    excludes = {'verification.json', 'completion.json'}
    artifacts = []
    for root in (RES, REP, OPS):
        for path in [root, *root.rglob('*')]:
            assert not path.is_symlink(), str(path)
            path.chmod(0o700 if path.is_dir() else 0o600)
            assert path.stat().st_uid == os.getuid(), str(path)
            if path.is_file() and path.name not in excludes:
                if path.suffix == '.json':
                    json.loads(path.read_text())
                artifacts.append({'path':str(path), 'bytes':path.stat().st_size, 'sha256':sha(path)})
    checks = {
        'selected_variants':selected,
        'recorded_reference_matches':sum(c['reference_check'] == 'matched_exact_GRCh38_FASTA' for c in candidates),
        'vep_completed':summary['vep_completed'],
        'primary_transcript_rows':summary['primary_transcript_rows'],
        'clinvar_exact_snapshot_matches':summary['clinvar_exact_matches'],
        'gnomad_identity_and_returned_AF_checked':selected,
        'source_metadata':'unchanged',
        'permissions':'0700 directories; 0600 files; current owner',
        'json_parsing':'all generated JSON parsed',
        'inheritance_annotations':'user supplied; not independently verified'}
    verification = {'status':'passed', 'checked_utc':datetime.now(timezone.utc).isoformat(),
        'checks':checks, 'artifact_count':len(artifacts), 'artifacts':artifacts}
    save(OPS/'verification.json', verification)
    report = REP/'VEP_PRIMARY_PREDICTIONS.md'
    assert report.is_file()
    save(OPS/'completion.json', dict(status='complete', report=str(report),
        verification=str(OPS/'verification.json'), verification_sha256=sha(OPS/'verification.json'),
        selected=selected, vep_completed=summary['vep_completed'],
        limitations=provenance['limitations']))
    print(json.dumps({'status':'passed', 'artifact_count':len(artifacts), 'report':str(report)}))


if __name__ == '__main__':
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['select','query','query_archive','gnomad','gnomad_qc','assemble','verify'])
    parser.add_argument('--manifest', type=Path, required=True, help='Private JSON input/configuration; see CANDIDATE_CONFIG.md')
    args = parser.parse_args()
    configure_manifest(args.manifest)
    mode = args.mode
    if mode in {'query','query_archive','gnomad','gnomad_qc'} and os.environ.get('MVA_ALLOW_NETWORK') != '1':
        raise SystemExit('set MVA_ALLOW_NETWORK=1 for opt-in public annotation queries')
    for path in (RES,REP,OPS):
        path.mkdir(parents=True, exist_ok=True, mode=0o700)
    {'select':select, 'query':query, 'query_archive':lambda:query(True), 'gnomad':gnomad, 'gnomad_qc':lambda:gnomad(True), 'assemble':assemble, 'verify':verify}[mode]()
