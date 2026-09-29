#!/usr/bin/env python3
"""Plot already-computed synthetic mixtures, never patient mosaic fractions."""
import os
import argparse
import csv
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--table', required=True, type=Path)
p.add_argument('--out', required=True, type=Path)
a = p.parse_args()
os.umask(0o077)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
with a.table.open() as f:
    rows = {r['scenario']: r for r in csv.DictReader(f, delimiter='\t')}
keys = ['coherent gain, f=0', 'coherent gain, f=0.2',
        'coherent loss, f=0.2', 'equal AAB/ABB/A/B mixture, f=0.40',
        'near-cancellation, g=0.11 l=0.09 balanced haplotypes']
labels = ['All diploid', '20% single gain clone', '20% single loss clone',
          '40% balanced gains/losses', '11% gains + 9% losses\n(haplotype-balanced)']
colors = ['#596579', '#be6433', '#be6433', '#147d87', '#147d87']
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
for ax, field, title, baseline, limits in zip(
        axes, ['expected_relative_RD', 'expected_BAF_B'],
        ['Relative chromosome depth', 'B-allele fraction'],
        [1., .5], [(.84, 1.18), (.40, .53)]):
    values = [float(rows[k][field]) for k in keys]
    ax.axvline(baseline, color='#9ca3af', linestyle='--', linewidth=1)
    ax.scatter(values, range(len(keys)), c=colors, s=80, zorder=3)
    for y, val in enumerate(values):
        ax.annotate(f'{val:.3f}', (val, y), xytext=(8, 5), textcoords='offset points', fontsize=9)
    ax.set_xlim(*limits)
    ax.set_title(title, fontsize=11)
    ax.grid(axis='y', alpha=.16)
    ax.spines[['top','right','left']].set_visible(False)
axes[0].set_yticks(range(len(keys)), labels)
axes[0].invert_yaxis()
fig.suptitle('Different cell mixtures can produce the same bulk-DNA signal', fontsize=13)
fig.text(.5, .025, 'Synthetic copy-count model — not patient data, a burden estimate, or a detection limit.',
         ha='center', fontsize=9)
fig.tight_layout(rect=(0, .065, 1, .93))
a.out.mkdir(parents=True, exist_ok=False, mode=0o700)
for ext in ('png', 'pdf'):
    target = a.out / f'bulk_mixture_examples.{ext}'
    fig.savefig(target, dpi=170, bbox_inches='tight')
    target.chmod(0o600)
