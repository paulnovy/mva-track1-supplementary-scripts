#!/usr/bin/env python3
"""Local, read-only source review; only new interpretation artifacts are written."""
import csv
import collections
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

os.umask(0o077)
ROOT = Path(os.environ.get('MVA_DATA_ROOT', '.')).resolve()
OPS = ROOT / 'operations/mocha-interpretation-20260910'
OUT = ROOT / 'results/mocha-interpretation-20260910'
SRC = ROOT / 'results/mocha-wgs-20260909/result'
if OUT.exists():
    raise SystemExit(f"refusing existing output directory: {OUT}")

def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def num(value):
    return float(value)

def merge(intervals):
    result = []
    for start, end in sorted(intervals):
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return result

def overlap(start, end, intervals):
    return sum(max(0, min(end, b)-max(start, a)) for a, b in intervals)

calls = table(SRC / 'calls.tsv')
stats = table(SRC / 'stats.tsv')[0]
rdroot = ROOT / 'results/rd-baf-followup-20260909/revision2'
bins = table(rdroot / 'bins.masked.tsv')
segments = table(rdroot / 'segments.tsv')
masks = collections.defaultdict(list)
for line in (ROOT / 'results/mocha-wgs-20260909/work/excluded-regions.sorted.bed').read_text().splitlines():
    chrom, start, end = line.split('\t')[:3]
    masks[chrom].append((int(start), int(end)))
masks = {chrom: merge(items) for chrom, items in masks.items()}

# These are the special GRCh38 regions explicitly suggested by pinned README.
special = {'MHC': ('chr6', 27518931, 33480487), 'KIR': ('chr19', 54071492, 54992731)}
sample_ok = num(stats['call_rate']) >= .97 and num(stats['baf_auto']) <= .03
annotated = []
for i, raw in enumerate(calls, 1):
    x = dict(raw)
    chrom = x['chrom']
    # MoChA emits these same boundaries to its BED writer, with length=end-beg.
    # Preserve the native half-open-style span rather than silently subtract 1.
    start, end = int(x['beg_GRCh38']), int(x['end_GRCh38'])
    assert end-start == int(x['length'])
    bdev, rc, length = num(x['bdev']), num(x['rel_cov']), int(x['length'])
    phase, conc = num(x['lod_baf_phase']), num(x['lod_baf_conc'])
    gates = {
        'sample_qc': sample_ok,
        'not_CNP': not x['type'].startswith('CNP'),
        'marker_spacing_or_phase_concordance': (
            ('X' in chrom and x['computed_gender']=='M') or bdev < .1
            or num(x['n50_hets']) < 2e5 or (math.isfinite(conc) and conc > 10)),
        'bdev_precision_or_strong_phase': (
            math.isfinite(num(x['bdev_se'])) or (math.isfinite(phase) and phase > 10)),
        'not_germline_like_gain': (rc < 2.1 or bdev < .05
            or (length > 5e5 and bdev < .1 and rc < 2.5)
            or (length > 5e6 and bdev < .15)),
    }
    matching_bins = [b for b in bins if b['chrom']==chrom and overlap(start,end,[(int(b['start']),int(b['end']))])]
    matching_segments = [s for s in segments if s['chrom']==chrom and overlap(start,end,[(int(s['start']),int(s['end']))])]
    mask_bp = overlap(start,end,masks.get(chrom,[]))
    special_names = [k for k,(c,a,b) in special.items() if c==chrom and overlap(start,end,[(a,b)])]
    x.update({
        'review_id': f'MC{i:02d}',
        'author_example_filter_pass': all(gates.values()),
        'filter_failed_gates': ';'.join(k for k,v in gates.items() if not v),
        'span_masked_bp': mask_bp,
        'span_masked_fraction': round(mask_bp/length,6),
        'special_region_overlap': ';'.join(special_names),
        'rd_overlapping_bins': len(matching_bins),
        'rd_overlapping_clean_bins': sum(b['clean'].lower() in ['true','1'] for b in matching_bins),
        'rd_overlapping_segments': len(matching_segments),
        'rd_joint_supported_overlaps': sum(s['joint_empirical_support'].lower() in ['true','1'] for s in matching_segments),
        'review_category': ('known_CNP_context_not_mCA_evidence' if x['type'].startswith('CNP')
                            else 'loss_without_heterozygous_support' if x['type']=='Loss' and int(x['n_hets'])==0
                            else 'gain_germline_like_or_mapping_context_unresolved'),
    })
    annotated.append(x)

# Reproduce the actual author's AWK example, not the looser README prose.
readme = (OPS / 'upstream-README.md').read_text()
block = readme.split('Filter callset\n==============',1)[1].split('```',2)[1]
program = block.split("awk -F \"\\t\" '",1)[1].split("' \\\n",1)[0]
(OPS / 'author-filter.awk').write_text(program+'\n')
run = subprocess.run(['awk','-F','\t','-f',str(OPS/'author-filter.awk'),str(SRC/'stats.tsv'),str(SRC/'calls.tsv')],capture_output=True,text=True,check=True)
filtered = list(csv.DictReader(run.stdout.splitlines(),delimiter='\t'))
assert len(filtered)==sum(x['author_example_filter_pass'] for x in annotated)
(OUT / 'calls.author-example-filter.tsv').write_text(run.stdout)
with (OUT / 'calls.review.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(annotated[0]),delimiter='\t');writer.writeheader();writer.writerows(annotated)

summary = {
    'status':'scientific_review_complete',
    'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'raw_calls':len(calls),
    'type_counts':dict(collections.Counter(x['type'] for x in calls)),
    'author_example_filter_pass_count':len(filtered),
    'author_example_is_not_clinical_validation':True,
    'filter_verification':'pinned upstream AWK example equals independent numeric Python gates',
    'non_CNP_failure_reasons':dict(collections.Counter(x['filter_failed_gates'] for x in annotated if not x['type'].startswith('CNP'))),
    'phase_lod_gt5':sum(num(x['lod_baf_phase'])>5 for x in calls),
    'phase_lod_gt10':sum(num(x['lod_baf_phase'])>10 for x in calls),
    'zero_heterozygote_calls':sum(int(x['n_hets'])==0 for x in calls),
    'calls_below_10_heterozygotes':sum(int(x['n_hets'])<10 for x in calls),
    'special_region_counts':dict(collections.Counter(x['special_region_overlap'] for x in annotated if x['special_region_overlap'])),
    'calls_with_majority_span_masked':sum(x['span_masked_fraction']>.5 for x in annotated),
    'rd_overlapping_call_count':sum(x['rd_overlapping_segments']>0 for x in annotated),
    'rd_joint_supported_overlapping_call_count':sum(x['rd_joint_supported_overlaps']>0 for x in annotated),
    'stats':{k:v for k,v in stats.items() if k not in ['sample_id','computed_gender']},
    'no_confirmed_mosaic_events_from_this_review':True,
    'MVA_excluded':False,
    'cell_fraction_clinically_calibrated':False,
    'source_calls_modified':False,
    'core_or_other_analysis_rerun':False,
    'raw_patient_data_sent_to_external_services':False,
    'sources':[str(SRC/'calls.tsv'),str(SRC/'stats.tsv'),str(rdroot/'summary.json'),str(rdroot/'bins.masked.tsv'),str(rdroot/'segments.tsv')],
    'method_reference':'https://github.com/freeseek/mocha/blob/95686b7/README.md#filter-callset',
}
(OUT / 'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['sources','stats']},ensure_ascii=False))
print('Non-CNP review aggregates (no genotype/reads):')
for x in annotated:
    if not x['type'].startswith('CNP'):
        print(json.dumps({k:x[k] for k in ['review_id','chrom','type','length','n_sites','n_hets','span_masked_fraction','special_region_overlap','rd_overlapping_clean_bins','rd_overlapping_segments','rd_joint_supported_overlaps','filter_failed_gates']}))
