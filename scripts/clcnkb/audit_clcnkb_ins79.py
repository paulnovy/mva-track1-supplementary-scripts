#!/usr/bin/env python3
"""Bounded public exact-allele annotation for the reconstructed CLCNKB INS79."""
from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys

import pysam


ROOT = Path(os.environ.get('MVA_DATA_ROOT', '.')).resolve()
RUN = 'clcnkb-ins79-annotation-20260911'
OUT = {kind: ROOT / kind / RUN for kind in ('results', 'reports', 'operations')}
REF = ROOT / 'references/GRCh38_standard/GRCh38_standard.fa'
CLINVAR = ROOT / 'resources/clinvar/prepared-grch38/clinvar.normalized.vcf.gz'
API_DIR = Path(__file__).resolve().parents[1] / 'public_api'
VEP_WRAPPER = Path(os.environ.get('ENSEMBL_VEP_WRAPPER', API_DIR / 'ensembl_api.py'))
GNOMAD_WRAPPER = Path(os.environ.get('GNOMAD_WRAPPER', API_DIR / 'get_variant_frequency.py'))
DBSNP_WRAPPER = Path(os.environ.get('DBSNP_WRAPPER', API_DIR / 'dbsnp_cli.py'))
# Public gnomAD v4 allele rs1553127751 (exact identity documented in METHODS.md).
INSERTED = 'GCCATTATTTTTTCCTGCCCAGACAATGCCCATGCAGTGATCTGGGCCCCCAAGGACCCAGCTTCACCCCCACAGCACC'
ALLELES = [
    {'id': 'CLCNKB_INS79', 'chrom': 'chr1', 'pos': 16050023, 'ref': 'T', 'alt': 'T' + INSERTED,
     'role': 'reconstructed 79-bp intron-10 insertion; no calibrated genotype/copy state'},
    {'id': 'nearby_SNV_control', 'chrom': 'chr1', 'pos': 16049963, 'ref': 'T', 'alt': 'C',
     'role': 'nearby small-variant control; phase with insertion is unresolved'},
]


def save(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    path.chmod(0o600)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimal(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Remove only shared VCF padding; never shift a repeat representation."""
    while ref and alt and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref, alt


def vep_variant(a: dict) -> str:
    # VEP accepts VCF-padded insertions. Keeping the anchor makes identity audit simple.
    return f"{a['chrom'][3:]}:{a['pos']}:{a['ref']}:{a['alt']}"


def fingerprint(path: Path) -> dict:
    s = path.stat()
    return {'path': str(path), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns,
            'inode': s.st_ino, 'mode': oct(s.st_mode & 0o777)}


def check_inputs() -> list[dict]:
    fa = pysam.FastaFile(str(REF))
    out = []
    for a in ALLELES:
        observed = fa.fetch(a['chrom'], a['pos'] - 1, a['pos'] - 1 + len(a['ref']))
        assert observed == a['ref'], (a['id'], observed, a['ref'])
        out.append({**a, 'reference_check': 'matched_exact_GRCh38_available_no-alt_reference',
                    'minimal_event': minimal(a['pos'], a['ref'], a['alt'])})
    # This insertion has no one-base equivalent left shift at the submitted anchor.
    preceding = fa.fetch('chr1', 16050021, 16050022)
    assert INSERTED[-1] != preceding
    assert len(INSERTED) == 79
    fa.close()
    return out


def local_clinvar(a: dict) -> dict:
    vf = pysam.VariantFile(str(CLINVAR))
    matches = []
    for rec in vf.fetch(a['chrom'], a['pos'] - 1, a['pos']):
        if rec.pos == a['pos'] and rec.ref == a['ref'] and rec.alts == (a['alt'],):
            matches.append({'variation_id': rec.id, 'chrom': rec.chrom, 'pos': rec.pos,
                            'ref': rec.ref, 'alt': rec.alts[0], 'info': dict(rec.info)})
    vf.close()
    return {'snapshot': 'ClinVar GRCh38 2026-09-05 locally prepared exact-normalized VCF',
            'exact_match_count': len(matches), 'matches': matches,
            'interpretation': 'exact record present' if matches else 'no exact record in this snapshot; not evidence of absence from ClinVar or benignity'}


def call_vep(a: dict, archive: bool) -> dict:
    spec = importlib.util.spec_from_file_location('ensembl_skill_wrapper', VEP_WRAPPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    endpoint = 'https://sep2025.rest.ensembl.org' if archive else mod.BASE_URL
    # Current Ensembl had known failures in the preceding bounded audit; do not
    # let one probe block the pinned archive fallback.
    mod._CLIENT_REGULAR = mod.http_client.HttpClient(endpoint, qps=2, timeout=12 if not archive else 25, max_retries=0,
        user_agent='CLCNKB-INS79-local-annotation/1.0', default_headers={'Content-Type': 'application/json'})
    suffix = 'release115' if archive else 'current'
    raw = OUT['results'] / 'vep' / suffix / f"{a['id']}.json"
    log = raw.with_suffix('.stdout.txt')
    raw.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    args = argparse.Namespace(variant_str=vep_variant(a), species='human', assembly='GRCh38', output=str(raw))
    buffer = io.StringIO()
    result = {'endpoint': endpoint, 'wrapper': str(VEP_WRAPPER), 'wrapper_sha256': digest(VEP_WRAPPER),
              'variant_input': args.variant_str, 'raw': str(raw)}
    try:
        with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
            mod.cmd_vep(args)
        result['status'] = 'completed'
    except Exception as exc:
        result.update(status='failed', error=repr(exc))
    log.write_text(buffer.getvalue())
    log.chmod(0o600)
    if raw.exists(): raw.chmod(0o600)
    return result


def run_cmd(command: list[str], output: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    run = subprocess.run(command, text=True, capture_output=True)
    log = output.with_suffix(output.suffix + '.stdout.txt')
    log.write_text(run.stdout + ('\nSTDERR:\n' + run.stderr if run.stderr else ''))
    log.chmod(0o600)
    if output.exists(): output.chmod(0o600)
    return {'command': command, 'status': 'completed' if run.returncode == 0 else 'failed',
            'returncode': run.returncode, 'raw': str(output), 'log': str(log)}


def run() -> None:
    os.umask(0o077)
    for p in OUT.values(): p.mkdir(parents=True, exist_ok=True, mode=0o700)
    source_before = [fingerprint(p) for p in (REF, CLINVAR)]
    alleles = check_inputs()
    for a in alleles: a['clinvar_local'] = local_clinvar(a)
    save(OUT['results'] / 'allele_identity.json', alleles)

    commands = {'vep_current': [], 'vep_release115': [], 'gnomad_r4': [], 'dbsnp': []}
    for a in alleles:
        commands['vep_current'].append({'id': a['id'], **call_vep(a, archive=False)})
    # Archive fallback is run for both exact alleles, providing a single release provenance even if current succeeds.
    for a in alleles:
        commands['vep_release115'].append({'id': a['id'], **call_vep(a, archive=True)})
    for a in alleles:
        out = OUT['results'] / 'gnomad_r4' / f"{a['id']}.json"
        commands['gnomad_r4'].append({'id': a['id'], **run_cmd(
            ['uv', 'run', str(GNOMAD_WRAPPER), '--variant_id', f"1-{a['pos']}-{a['ref']}-{a['alt']}", '--dataset', 'gnomad_r4', '--output', str(out)], out)})
    for a in alleles:
        out = OUT['results'] / 'dbsnp' / f"{a['id']}.json"
        commands['dbsnp'].append({'id': a['id'], **run_cmd(
            ['uv', 'run', str(DBSNP_WRAPPER), 'resolve-variant', '1', str(a['pos']), a['ref'], a['alt'], '--output', str(out)], out)})
    save(OUT['operations'] / 'commands.json', commands)
    save(OUT['operations'] / 'source_before.json', source_before)
    save(OUT['operations'] / 'source_after.json', [fingerprint(p) for p in (REF, CLINVAR)])
    assemble(alleles, commands)
    verify(alleles, commands, source_before)


def vep_summary(path: Path, a: dict) -> dict:
    if not path.exists(): return {'status': 'no_response_file'}
    raw = json.loads(path.read_text())
    if not isinstance(raw, list) or not raw: return {'status': 'empty_or_nonlist_response'}
    top = raw[0]
    expected_pos, expected_ref, expected_alt = minimal(a['pos'], a['ref'], a['alt'])
    observed_ref, observed_alt = top.get('allele_string', '/').split('/', 1)
    identity = (top.get('seq_region_name') == '1' and top.get('assembly_name') == 'GRCh38' and
                top.get('start') == expected_pos and observed_ref.replace('-', '') == expected_ref and
                observed_alt.replace('-', '') == expected_alt)
    tc = [x for x in top.get('transcript_consequences', []) if x.get('gene_symbol') == 'CLCNKB']
    # The archived VEP response does not carry a canonical flag. The reconstruction
    # independently fixed ENST00000375679.9 as the Ensembl canonical transcript.
    primary = [x for x in tc if x.get('transcript_id') == 'ENST00000375679']
    return {'status': 'returned_identity_checked' if identity else 'returned_identity_mismatch_or_incomplete',
            'most_severe_consequence': top.get('most_severe_consequence'), 'input': top.get('input'),
            'primary_CLCNKB_transcripts': primary}


def gnomad_summary(path: Path, a: dict) -> dict:
    if not path.exists(): return {'status': 'query_failed_or_no_file; not_zero_AF'}
    raw = json.loads(path.read_text())
    v = raw.get('data', {}).get('variant')
    if v is None: return {'status': 'not_returned_by_query; not_zero_AF', 'errors': raw.get('errors')}
    exact = v.get('variant_id') == f"1-{a['pos']}-{a['ref']}-{a['alt']}"
    subsets = {}
    for name in ('exome', 'genome', 'joint'):
        d = v.get(name)
        if not d: subsets[name] = {'status': 'not_returned; not_zero_AF'}; continue
        an, ac = d.get('an'), d.get('ac')
        subsets[name] = {'status': 'returned', 'ac': ac, 'an': an,
                         'af': (ac / an if an else None), 'homozygote_count': d.get('homozygote_count'),
                         'faf95': d.get('faf95'), 'faf99': d.get('faf99')}
    return {'status': 'returned_exact_identity' if exact else 'returned_identity_mismatch', 'rsids': v.get('rsids'), 'subsets': subsets}


def assemble(alleles: list[dict], commands: dict) -> None:
    rows = []
    for a in alleles:
        archive = vep_summary(OUT['results'] / 'vep' / 'release115' / f"{a['id']}.json", a)
        current = vep_summary(OUT['results'] / 'vep' / 'current' / f"{a['id']}.json", a)
        pop = gnomad_summary(OUT['results'] / 'gnomad_r4' / f"{a['id']}.json", a)
        db = OUT['results'] / 'dbsnp' / f"{a['id']}.json"
        db_value = json.loads(db.read_text()) if db.exists() else {'status': 'query_failed_or_no_file'}
        rows.append({'allele': a, 'vep_current': current, 'vep_release115': archive, 'gnomad_r4': pop, 'dbsnp': db_value})
    save(OUT['results'] / 'annotation_summary.json', rows)
    ins = rows[0]
    primary = ins['vep_release115'].get('primary_CLCNKB_transcripts', [])
    transcript_text = json.dumps(primary, indent=2) if primary else 'No archive primary-transcript response available.'
    joint = ins['gnomad_r4'].get('subsets', {}).get('joint', {})
    frequency_text = (f"joint AC {joint.get('ac')}, AN {joint.get('an')}, AF {joint.get('af')}, "
                      f"homozygote count {joint.get('homozygote_count')}" if joint.get('status') == 'returned'
                      else 'not returned; not frequency 0')
    dbsnp_rsids = [f"rs{str(rsid).removeprefix('rs')}" for rsid in ins['dbsnp'].get('rsids', [])]
    dbsnp_text = (f"Exact coordinate/REF/ALT resolution returned {', '.join(dbsnp_rsids)}."
                  if dbsnp_rsids else 'No rsID was returned; this does not establish novelness.')
    report = f'''# CLCNKB INS79 exact-allele annotation and population audit

**Scope:** local reconstruction plus bounded public exact-allele queries; only explicitly configured alleles are sent by online modes; raw reads and identifiers are not sent. This script does not submit competition entries or infer genotype/copy state.

## Exact event identity

- GRCh38 available no-alt reference: `chr1:16050023 T>{'T' + INSERTED}` (79 inserted bases), VCF-padded insertion representation.
- Local FASTA REF matched. Shared-padding minimisation and the immediate left-base check found no one-base equivalent left shift at this anchor.
- The reconstruction places the boundary 107 bp after CLCNKB exon 10 and 493 bp before exon 11 (intron 10). This is an annotation location, **not** evidence of no splice effect, benignity, or a calibrated genotype.
- `chr1:16049963 T>C` was queried separately as a nearby control. No phase is assumed.

## Sources and results

- **ClinVar:** exact local normalized GRCh38 snapshot dated 2026-09-05. INS79 exact matches: {ins['allele']['clinvar_local']['exact_match_count']}. A zero means no exact record in that snapshot, not that the allele is clinically benign or absent from all ClinVar material.
- **Ensembl VEP:** current official endpoint was attempted first; pinned fallback is the official September-2025 archive (Ensembl release 115). The raw response, endpoint, inputs, and identity checks are preserved in results. Archive primary CLCNKB transcript rows:

```json
{transcript_text}
```

- **gnomAD r4:** exact submitted event only; identity status {ins['gnomad_r4'].get('status')}, returned rsIDs {ins['gnomad_r4'].get('rsids')}; {frequency_text}. The supplied frequency wrapper did not request per-record filter fields, so QC-filter status is **not returned**, not PASS. A non-return is **not frequency 0**, particularly for this 79-bp insertion at a complex paralogous locus. AN is callability, not a coverage assay.
- **dbSNP:** {dbsnp_text} The rsID can include other alternate alleles; evidence here is restricted to the exact 79-bp insertion above.

## Interpretation limits

Interpret frequency only after checking the returned exact-allele identity and population fields above. No frequency or clinical classification is assumed when a query fails or returns a different allele. The INS79 evidence is sequence-level support for an intronic insertion. It does not establish splicing impact, clinical significance, chromosome-scale phase, zygosity, or CLCNKB copy context. Population catalogs may underrepresent long insertions and complex paralogous loci; this audit does not treat nearby/homologous alleles as equivalent to the exact event.

## Files

- `results/.../allele_identity.json`: exact sequence, local REF and ClinVar snapshot checks.
- `results/.../vep`, `gnomad_r4`, and `dbsnp`: raw bounded public responses and logs.
- `operations/.../commands.json`: endpoint, wrapper and command provenance.
- `operations/.../verification.json`: artifact hashes, excluding the manifest itself to avoid a circular hash.
'''
    p = OUT['reports'] / 'CLCNKB_INS79_EXACT_ALLELE_ANNOTATION.md'
    p.write_text(report)
    p.chmod(0o600)


def verify(alleles: list[dict], commands: dict, before: list[dict]) -> None:
    assert len(INSERTED) == 79
    assert alleles[0]['alt'] == alleles[0]['ref'] + INSERTED
    assert tuple(alleles[0]['minimal_event']) == (16050024, '', INSERTED)
    assert tuple(alleles[1]['minimal_event']) == (16049963, 'T', 'C')
    assert before == json.loads((OUT['operations'] / 'source_after.json').read_text())
    artifacts = []
    manifest = OUT['operations'] / 'verification.json'
    for root in OUT.values():
        for p in root.rglob('*'):
            assert not p.is_symlink()
            if p.is_dir(): p.chmod(0o700)
            else:
                p.chmod(0o600)
                if p != manifest:
                    artifacts.append({'path': str(p), 'sha256': digest(p), 'bytes': p.stat().st_size})
    save(manifest, {
        'status': 'passed', 'checked_utc': datetime.now(timezone.utc).isoformat(),
        'checks': ['79-bp inserted sequence length', 'VCF anchor/ALT identity', 'minimal event representation',
                   'local FASTA REF', 'distinct nearby-SNV control', 'source metadata unchanged',
                   'no symlinks; output permissions 0700 directories/0600 files'],
        'manifest_exclusion': 'verification.json is not hashed inside itself',
        'artifacts': artifacts, 'commands': {k: [x['status'] for x in v] for k, v in commands.items()}})


def repair_vep() -> None:
    """Save the successful archive response after a pre-save directory error."""
    alleles = json.loads((OUT['results'] / 'allele_identity.json').read_text())
    commands = json.loads((OUT['operations'] / 'commands.json').read_text())
    for a in alleles:
        if a['id'] != 'CLCNKB_INS79':
            continue
        repaired = call_vep(a, archive=True)
        assert repaired['status'] == 'completed', repaired
        for index, item in enumerate(commands['vep_release115']):
            if item['id'] == a['id']:
                previous = {k: v for k, v in item.items() if k != 'prior_attempts'}
                commands['vep_release115'][index] = {
                    'id': a['id'], **repaired,
                    'prior_attempts': item.get('prior_attempts', []) + [previous]}
    save(OUT['operations'] / 'commands.json', commands)
    assemble(alleles, commands)
    before = json.loads((OUT['operations'] / 'source_before.json').read_text())
    verify(alleles, commands, before)


def reassemble() -> None:
    alleles = json.loads((OUT['results'] / 'allele_identity.json').read_text())
    commands = json.loads((OUT['operations'] / 'commands.json').read_text())
    assemble(alleles, commands)
    verify(alleles, commands, json.loads((OUT['operations'] / 'source_before.json').read_text()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['run','repair-vep','reassemble'], default='run')
    mode = parser.parse_args().mode
    if mode != 'reassemble' and os.environ.get('MVA_ALLOW_NETWORK') != '1':
        raise SystemExit('set MVA_ALLOW_NETWORK=1 for opt-in public annotation queries')
    if mode == 'repair-vep':
        repair_vep()
    elif mode == 'reassemble':
        reassemble()
    else:
        run()
