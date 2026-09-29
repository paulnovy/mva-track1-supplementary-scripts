#!/usr/bin/env python3
"""Build canonical-panel inputs from saved Ensembl expanded gene JSON, offline."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene', action='append', required=True, help='SYMBOL=expanded-JSON')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    prepared = []
    for item in args.gene:
        symbol, filename = item.split('=', 1)
        path = Path(filename)
        gene = json.loads(path.read_text())
        if gene.get('assembly_name') != 'GRCh38' or gene.get('display_name') != symbol:
            raise ValueError(f'{symbol}: expected matching GRCh38 gene metadata')
        canonical_id = gene.get('canonical_transcript', '').split('.')[0]
        selected = [tx for tx in gene['Transcript'] if tx.get('is_canonical') or tx['id'] == canonical_id]
        if len(selected) != 1 or not selected[0].get('Translation'):
            raise ValueError(f'{symbol}: expected one translated canonical transcript')
        prepared.append((symbol, path, gene, selected[0]))
    args.out.mkdir(parents=True, exist_ok=False)
    provenance = []
    with (args.out / 'panel_gene_coordinates.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        for symbol, path, gene, tx in prepared:
            writer.writerow([symbol, gene['id'], gene['seq_region_name'], gene['start'], gene['end'], gene['strand'], gene.get('description', '')])
            (args.out / f'ensembl_{symbol}_canonical_structure.json').write_text(json.dumps(tx, indent=2) + '\n')
            provenance.append({'gene': symbol, 'source': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'canonical_transcript': tx['id'], 'transcript_version': tx.get('version')})
    (args.out / 'metadata_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
