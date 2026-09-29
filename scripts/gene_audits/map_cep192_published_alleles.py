#!/usr/bin/env python3
"""Map the two Guo et al. CEP192 alleles through the exact MANE transcript.

Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10716027/ (NM_032142.4).
This maps literature alleles, not genotypes in the current patient.
"""
import os
import argparse
import json
from pathlib import Path

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    os.umask(0o077)
    gene = json.loads(Path(args.metadata).read_text())
    transcripts = [t for t in gene['Transcript'] if any(
        m.get('refseq_match') == 'NM_032142.4' for m in t.get('MANE', []))]
    assert gene['id'] == 'ENSG00000101639' and gene['assembly_name'] == 'GRCh38'
    assert len(transcripts) == 1
    tx = transcripts[0]
    assert tx['strand'] == 1
    tr = tx['Translation']
    cds = sorted((max(e['start'],tr['start']),min(e['end'],tr['end']))
                 for e in tx['Exon'] if
                 max(e['start'],tr['start']) <= min(e['end'],tr['end']))
    # Ensembl genomic translation bounds include the terminal stop codon.
    assert sum(e-s+1 for s,e in cds) == 3*(tr['length']+1)
    rows = []
    with pysam.FastaFile(args.fasta) as ref:
        sequence = ''.join(ref.fetch('chr18',s-1,e) for s,e in cds)
        assert sequence[:3] == 'ATG' and sequence[-3:] in {'TAA','TAG','TGA'}
        for c, base, alt in ((1912,'C','T'), (5750,'A','G')):
            offset = 0
            positions = []
            for s,e in cds:
                if offset < c <= offset+e-s+1:
                    positions.append(s+c-offset-1)
                offset += e-s+1
            assert len(positions) == 1
            pos = positions[0]
            assert ref.fetch('chr18',pos-1,pos) == base
            rows.append(dict(transcript='NM_032142.4',
                ensembl=tx['id']+'.'+str(tx['version']), cdna=c,
                chrom='chr18',pos=pos,ref=base,alt=alt,reference_base=base,
                source='https://pmc.ncbi.nlm.nih.gov/articles/PMC10716027/',
                method='Sequential CDS positions from MANE expanded exons on positive strand; retained GRCh38 reference checked'))
    with Path(args.output).open('x') as target:
        target.write(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows))


if __name__ == '__main__':
    main()
