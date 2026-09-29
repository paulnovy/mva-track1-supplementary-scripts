#!/usr/bin/env python3
"""Verify reference/control sequences for a proposed protein-function experiment.

This does not predict stability, run an assay, or create expression-ready constructs.
"""
import os
import argparse
import hashlib
import json
from pathlib import Path


def fasta(path):
    records = {}
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            records[name] = ''
        elif line.strip():
            records[name] += line.strip()
    return records


def translate(seq):
    bases = 'TCAG'
    aa = ('FFLLSSSSYY**CC*W' 'LLLLPPPPHHQQRRRR'
          'IIIMTTTTNNKKSSRR' 'VVVVAAAADDEEGGGG')
    codons = [a+b+c for a in bases for b in bases for c in bases]
    table = dict(zip(codons, aa))
    if len(seq) % 3:
        raise ValueError('CDS length is not divisible by 3')
    return ''.join(table[seq[i:i+3]] for i in range(0, len(seq), 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--protein-fasta', required=True, type=Path)
    ap.add_argument('--cds-fasta', required=True, type=Path)
    ap.add_argument('--uniprot-json', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()
    os.umask(0o077)
    name = 'NM_001211.6_reference_GRCh38_local'
    protein = fasta(args.protein_fasta)[name]
    cds = fasta(args.cds_fasta)[name]
    uni = json.loads(args.uniprot_json.read_text())
    if 'results' in uni:
        uni = uni['results'][0]
    assert uni['primaryAccession'] == 'O60566'
    assert protein == uni['sequence']['value']
    assert len(protein) == 1050 and translate(cds) == protein + '*'
    # Genomic case allele is c.3006T>G; control codons are explicitly selected
    # protein-matched synthetic CDS changes, not claims about patient genotypes.
    controls = [('WT', None, None, None, 'reference'),
                ('N1002K', 1002, 'AAT', 'AAG', 'case protein hypothesis; untested'),
                ('L1012P', 1012, 'CTT', 'CCT', 'published low-abundance control'),
                ('Q921H', 921, 'CAG', 'CAT', 'published assay-normal control')]
    args.out.mkdir(parents=True, exist_ok=False, mode=0o700)
    manifest, protein_records, cds_records = [], [], []
    for label, pos, old, new, role in controls:
        mutant = cds
        if pos is not None:
            start = (pos-1)*3
            assert cds[start:start+3] == old
            mutant = cds[:start] + new + cds[start+3:]
        translated = translate(mutant)
        assert translated.endswith('*') and '*' not in translated[:-1]
        differences = [(i+1, a, b) for i,(a,b) in enumerate(zip(protein, translated[:-1])) if a != b]
        if pos is not None:
            assert differences == [(pos, label[0], label[-1])]
        else:
            assert not differences
        edits = [(i+1,a,b) for i,(a,b) in enumerate(zip(cds,mutant)) if a != b]
        manifest.append(dict(label=label,role=role,protein_changes=differences,cds_changes=edits,
                             protein_length=len(translated)-1,cds_length=len(mutant)))
        desc = 'in_silico_reference_derivative_NOT_experiment_or_expression_construct'
        protein_records.append(f'>{label} {desc}\n{translated[:-1]}\n')
        cds_records.append(f'>{label} {desc}\n{mutant}\n')
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
              (args.protein_fasta,args.cds_fasta,args.uniprot_json)}
    result = dict(status='sequence_checks_passed',assay_performed=False,stability_predicted=False,
                  uniprot_accession=uni['primaryAccession'],uniprot_version=uni['entryAudit'],
                  reference_identity=True,controls=manifest,
                  N1002_to_L1012_residue_difference=10,
                  N1002_in_CENPE_terminal_1031_1050_region=False,input_sha256=hashes,
                  caution='Protein-function controls do not establish endogenous splicing, phase, or patient causality.')
    (args.out/'control_map.json').write_text(json.dumps(result,indent=2)+'\n')
    (args.out/'proposed_control_proteins.fa').write_text(''.join(protein_records))
    (args.out/'proposed_control_cds.fa').write_text(''.join(cds_records))
    for p in args.out.iterdir():
        if p.is_file():
            p.chmod(0o600)
    print(json.dumps({'status':result['status'],'protein_identity':True,'controls':len(manifest)}))


if __name__ == '__main__':
    main()
