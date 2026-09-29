#!/bin/bash
set -euo pipefail
ref=${1:?usage: filter_existing_mtdna.sh REFERENCE.fa NATIVE.vcf.gz OUTPUT_DIR}
native=${2:?usage: filter_existing_mtdna.sh REFERENCE.fa NATIVE.vcf.gz OUTPUT_DIR}
output=${3:?usage: filter_existing_mtdna.sh REFERENCE.fa NATIVE.vcf.gz OUTPUT_DIR}
[[ -s "$native" && -s "$native.stats" && -s "$ref" ]]
[[ ! -e "$output" ]] || { echo "refusing existing output: $output" >&2; exit 73; }
mkdir -m 700 "$output"
gatk --java-options '-Xmx2g -XX:ActiveProcessorCount=1' FilterMutectCalls --mitochondria-mode -R "$ref" -V "$native" -O "$output/mtdna.filtered.vcf.gz"
bcftools index --csi "$output/mtdna.filtered.vcf.gz"
python3 - "$output/mtdna.filtered.vcf.gz" "$output/summary.json" "$native" <<'PY'
import json,pysam,collections
import sys
from datetime import datetime,timezone
counts=collections.Counter()
with pysam.VariantFile(sys.argv[1]) as f:
 for r in f:
  counts['records']+=1
  if tuple(r.filter)==('PASS',): counts['pass']+=1
  else: counts['filtered']+=1
  for k in r.filter: counts['filter:'+k]+=1
json.dump({'status':'complete','method':'GATK FilterMutectCalls mitochondria mode; existing Mutect2 output only','utc':datetime.now(timezone.utc).isoformat(),'counts':dict(counts),'limitations':['No shifted circular reference analysis','No contamination estimate or NUMT calibration','PASS does not establish clinical significance or validated low heteroplasmy detection'],'source':sys.argv[3]},open(sys.argv[2],'w'),indent=2)
PY
