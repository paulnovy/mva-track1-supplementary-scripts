#!/usr/bin/env bash
set -euo pipefail
umask 077

usage() {
  cat >&2 <<'EOF'
usage: run_taxonomy_control.sh --bam BAM --db KRAKEN_DB --output DIR --sample ID [--max-pairs N]

Run a bounded, non-diagnostic Kraken2 control on primary both-unmapped BAM pairs.
The output directory must not already exist. Run inside a container capped at
1 CPU/6 GiB; this script does not start Docker or systemd services.
EOF
  exit 64
}

bam=''; db=''; output=''; sample=''; max_pairs=1000000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) bam=${2:?}; shift 2 ;;
    --db) db=${2:?}; shift 2 ;;
    --output) output=${2:?}; shift 2 ;;
    --sample) sample=${2:?}; shift 2 ;;
    --max-pairs) max_pairs=${2:?}; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$bam" && -n "$db" && -n "$output" && -n "$sample" ]] || usage
[[ "$max_pairs" =~ ^[1-9][0-9]*$ ]] || { echo 'max-pairs must be positive' >&2; exit 64; }

command -v samtools >/dev/null || { echo 'samtools is required' >&2; exit 69; }
command -v kraken2 >/dev/null || { echo 'kraken2 is required' >&2; exit 69; }
[[ -s "$bam" && -s "$bam.bai" ]] || { echo 'BAM and BAI must be readable' >&2; exit 66; }
for f in hash.k2d opts.k2d taxo.k2d; do [[ -s "$db/$f" ]] || { echo "missing Kraken2 DB file: $f" >&2; exit 66; }; done
[[ ! -e "$output" ]] || { echo "refusing existing output (duplicate guard): $output" >&2; exit 73; }
parent=$(dirname "$output")
mkdir -p "$parent"
free=$(python3 -c 'import shutil,sys; print(shutil.disk_usage(sys.argv[1]).free)' "$parent")
[[ "$free" =~ ^[0-9]+$ && "$free" -ge $((60*1024**3)) ]] || { echo 'disk floor 60 GiB not available' >&2; exit 75; }
mkdir -m 700 "$output"
exec 9>"$output/.taxonomy.lock"
flock -n 9 || { echo 'taxonomy control lock is busy' >&2; exit 75; }

status="$output/status.json"; work="$output/.work"; log="$output/taxonomy.log"
mkdir -m 700 "$work"; touch "$log"; chmod 600 "$log"
database_version=$(python3 - "$db" <<'PY'
import json, pathlib, sys
p=pathlib.Path(sys.argv[1]).parent/'provenance.json'
try:
    print(json.loads(p.read_text()).get('version') or p.parent.name)
except Exception:
    print(pathlib.Path(sys.argv[1]).name)
PY
)

write_status() {
  local state=$1 detail=${2:-}
  STATUS_PATH="$status" STATUS_STATE="$state" STATUS_DETAIL="$detail" \
  STATUS_SAMPLE="$sample" STATUS_MAX_PAIRS="$max_pairs" python3 - <<'PY'
import json, os, tempfile
p=os.environ['STATUS_PATH']
d={'schema':'taxonomy_control_status_v1','status':os.environ['STATUS_STATE'],'sample':os.environ['STATUS_SAMPLE'],'max_pairs':int(os.environ['STATUS_MAX_PAIRS'])}
if os.environ.get('STATUS_DETAIL'): d['detail']=os.environ['STATUS_DETAIL']
fd,tmp=tempfile.mkstemp(prefix='.status.',dir=os.path.dirname(p),text=True); os.fchmod(fd,0o600)
with os.fdopen(fd,'w') as f: json.dump(d,f,sort_keys=True); f.write('\n')
os.replace(tmp,p)
PY
}

finish() {
  rc=$?
  if (( rc != 0 )); then
    if [[ "${NO_INPUT:-}" == 1 ]]; then write_status no_input 'no primary both-unmapped pairs available'; else write_status failed "exit_code=$rc" || true; fi
  fi
  rm -rf -- "$work"
  exit "$rc"
}
trap finish EXIT

write_status waiting_for_resources 'launch only after PRS reservation is inactive'
[[ "${MVA_TAXONOMY_RESOURCES_READY:-}" == 1 ]] || { echo 'resource gate not satisfied; set MVA_TAXONOMY_RESOURCES_READY=1 only after PRS is inactive' >&2; exit 75; }
write_status extracting_primary_both_unmapped
filtered="$work/primary-both-unmapped.bam"; collated="$work/collated.bam"
r1all="$work/all.R1.fastq"; r2all="$work/all.R2.fastq"

# -f 12 selects both-unmapped pairs; -F 2304 excludes secondary/supplementary.
samtools view -@ 1 -bh -f 12 -F 2304 "$bam" '*' > "$filtered" 2>>"$log"
records=$(samtools view -@ 1 -c "$filtered" 2>>"$log")
(( records % 2 == 0 )) || { echo 'odd primary pair record count' >&2; exit 65; }
total_pairs=$((records / 2))
samtools collate -@ 1 -u -o "$collated" "$filtered" 2>>"$log"
samtools fastq -@ 1 -n -1 "$r1all" -2 "$r2all" -0 /dev/null -s /dev/null "$collated" 2>>"$log"
pair_lines=$(wc -l < "$r1all")
(( pair_lines % 4 == 0 && pair_lines == $(wc -l < "$r2all") )) || { echo 'paired FASTQ invariant failed' >&2; exit 65; }
extracted_pairs=$((pair_lines / 4)); (( extracted_pairs == total_pairs )) || { echo 'collated extraction lost a pair' >&2; exit 65; }
if (( extracted_pairs == 0 )); then NO_INPUT=1; exit 66; fi
selected_pairs=$extracted_pairs
if (( selected_pairs > max_pairs )); then selected_pairs=$max_pairs; fi

write_status running "selected_pairs=$selected_pairs,total_primary_pairs=$total_pairs"
sample_r1="$work/sample.R1.fastq.gz"; sample_r2="$work/sample.R2.fastq.gz"
head -n $((selected_pairs * 4)) "$r1all" | gzip -n -c > "$sample_r1"
head -n $((selected_pairs * 4)) "$r2all" | gzip -n -c > "$sample_r2"

report="$output/${sample}.kraken2.report"
write_status classifying "selected_pairs=$selected_pairs,total_primary_pairs=$total_pairs"
kraken2 --db "$db" --paired --threads 1 --memory-mapping --report "$report" --output /dev/null "$sample_r1" "$sample_r2" >>"$log" 2>&1

[[ -s "$report" ]] || { echo 'Kraken2 produced no report' >&2; exit 65; }
read_length=$(python3 - "$sample_r1" <<'PY'
import gzip, sys
with gzip.open(sys.argv[1], 'rt') as f:
    f.readline()
    print(len(f.readline().rstrip('\n')))
PY
)
report_counts=$(awk -F '\t' '$5==0 {u+=$2} $5==1 {r+=$2} END {printf "%d %d %d\n",u+0,r+0,(u+r)+0}' "$report")
read -r unclassified_count root_count report_total <<<"$report_counts"
(( report_total == selected_pairs )) || { echo "Kraken report count mismatch: $report_total != $selected_pairs" >&2; exit 65; }

bracken_status=not_run
if command -v bracken >/dev/null 2>&1 && [[ "$read_length" == 150 && -s "$db/database150mers.kmer_distrib" ]]; then
  if bracken -d "$db" -i "$report" -o "$output/${sample}.bracken.tsv" -r 150 -l S -t 10 >>"$log" 2>&1; then bracken_status=applied; else bracken_status=failed_nonblocking; fi
fi

STATUS_PATH="$status" STATUS_SAMPLE="$sample" STATUS_TOTAL_PAIRS="$total_pairs" STATUS_EXTRACTED_PAIRS="$extracted_pairs" \
STATUS_SELECTED_PAIRS="$selected_pairs" STATUS_MAX_PAIRS="$max_pairs" STATUS_BRACKEN="$bracken_status" STATUS_REPORT="$report" \
STATUS_DATABASE_VERSION="$database_version" STATUS_UNCLASSIFIED="$unclassified_count" STATUS_ROOT="$root_count" STATUS_TOTAL="$report_total" STATUS_READ_LENGTH="$read_length" python3 - <<'PY'
import json, os, tempfile
p=os.environ['STATUS_PATH']
d={'schema':'taxonomy_control_status_v1','status':'complete','sample':os.environ['STATUS_SAMPLE'],'input_selection':'primary both-unmapped BAM pairs from region *; first N pairs after name-collation','database_version':os.environ.get('STATUS_DATABASE_VERSION','unknown'),'total_primary_both_unmapped_pairs':int(os.environ['STATUS_TOTAL_PAIRS']),'extracted_pairs':int(os.environ['STATUS_EXTRACTED_PAIRS']),'selected_pairs':int(os.environ['STATUS_SELECTED_PAIRS']),'max_pairs':int(os.environ['STATUS_MAX_PAIRS']),'kraken_report':os.environ['STATUS_REPORT'],'report_unclassified':int(os.environ['STATUS_UNCLASSIFIED']),'report_root':int(os.environ['STATUS_ROOT']),'report_total':int(os.environ['STATUS_TOTAL']),'read_length':int(os.environ['STATUS_READ_LENGTH']),'bracken':'none' if os.environ['STATUS_BRACKEN']=='not_run' else os.environ['STATUS_BRACKEN'],'interpretation':'non-diagnostic control; indicative microbial composition only; not specimen provenance'}
fd,tmp=tempfile.mkstemp(prefix='.status.',dir=os.path.dirname(p),text=True); os.fchmod(fd,0o600)
with os.fdopen(fd,'w') as f: json.dump(d,f,sort_keys=True); f.write('\n')
os.replace(tmp,p)
PY
