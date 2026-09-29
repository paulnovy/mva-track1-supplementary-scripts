#!/usr/bin/env bash
set -euo pipefail

phenotype_root=${1:?prepared phenotype root required}
result_root=${2:?result root required}
operation_root=${3:?operation root required}
resource_root=${4:?Exomiser resource root required}
vcf=${5:?normalized VCF required}
adapter_root=${6:?wgs-cockpit root required}
sample=${7:-CASE}

umask 077
mkdir -p "$result_root" "$operation_root"
chmod 700 "$result_root" "$operation_root"
exec 9>"$operation_root/exomiser-sensitivity.lock"
flock -n 9 || { echo "another Exomiser sensitivity runner holds the lock" >&2; exit 75; }

manifest="$resource_root/resource-manifest.json"
jar="$resource_root/exomiser-cli-15.1.0/exomiser-cli-15.1.0.jar"
java="$resource_root/jre21/bin/java"
[[ -s "$resource_root/resources.complete" && -s "$manifest" && -s "$jar" && -x "$java" ]]
[[ -s "$vcf" ]]

available_kib=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)
(( available_kib >= 20 * 1024 * 1024 )) || { echo "less than 20 GiB RAM available" >&2; exit 76; }
jar_sha=$(sha256sum "$jar" | awk '{print $1}')
manifest_sha=$(sha256sum "$manifest" | awk '{print $1}')

for profile in baseline full axis; do
  input="$phenotype_root/$profile"
  output="$result_root/$profile"
  log="$operation_root/$profile"
  mkdir -p "$output" "$log"
  chmod 700 "$output" "$log"
  if [[ -s "$output/complete.json" ]]; then
    echo "[$profile] verified completion marker already present" >&2
    continue
  fi
  [[ -s "$input/phenopacket.json" && -s "$input/family.ped" ]]
  profile_sha=$(sha256sum "$input/phenopacket.json" | awk '{print $1}')
  work="$output/work"
  mkdir -p "$work"
  start=$(date -u +%FT%TZ)
  (
    cd "$work"
    export PATH="$(dirname "$java"):$PATH"
    export PYTHONPATH="$adapter_root/apps/api"
    export WGS_EXOMISER_JAR="$jar"
    export WGS_EXOMISER_JAR_SHA256="$jar_sha"
    export WGS_EXOMISER_VERSION=15.1.0
    export WGS_RARE_DISEASE_RESOURCE_MANIFEST="$manifest"
    export WGS_PHENOPACKET_PROFILE_SHA256="$profile_sha"
    export WGS_EXOMISER_JAVA_OPTS="-Xmx16g -XX:+UseG1GC -XX:ActiveProcessorCount=5 -Dexomiser.data-directory=$resource_root/data -Dexomiser.hg38.data-version=2512 -Dexomiser.phenotype.data-version=2512"
    /usr/bin/time -v "$adapter_root/pipelines/nextflow/scripts/run_rare_disease_stage.sh" \
      "${sample}-${profile}" "$input/phenopacket.json" "$vcf" "$input/family.ped" hg38
  ) >"$log/stdout.log" 2>"$log/stderr.log"
  end=$(date -u +%FT%TZ)
  genes="$work/rare-disease/native/${sample}-${profile}-exomiser.genes.tsv"
  variants="$work/rare-disease/native/${sample}-${profile}-exomiser.variants.tsv"
  engine="$work/${sample}-${profile}.rare_disease.engine.json"
  [[ -s "$genes" && -s "$variants" && -s "$engine" ]]
  python3 - "$profile" "$start" "$end" "$profile_sha" "$manifest_sha" "$vcf" "$genes" "$variants" "$engine" "$output/complete.json" <<'PY'
import csv, hashlib, json, os, pathlib, sys

profile, started, completed, profile_sha, manifest_sha, vcf, genes, variants, engine, target = sys.argv[1:]
def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()
def rows(path):
    with open(path, encoding='utf-8') as handle:
        return sum(1 for line in handle if line and not line.startswith('#'))
payload = {
    'schema': 'mva1_exomiser_sensitivity_completion_v1', 'status': 'complete', 'profile': profile,
    'started_utc': started, 'completed_utc': completed, 'engine': 'Exomiser', 'engine_version': '15.1.0',
    'phenotype_data_version': '2512', 'genome_data_version': '2512', 'assembly': 'GRCh38',
    'preset': 'EXOME', 'analysis_mode': 'PASS_ONLY', 'research_only': True,
    'ranking_is_not_probability_or_diagnosis': True,
    'input_sha256': {'phenopacket': profile_sha, 'resource_manifest': manifest_sha, 'vcf': digest(vcf)},
    'outputs': {
        'gene_rows': rows(genes), 'variant_rows': rows(variants),
        'genes_sha256': digest(genes), 'variants_sha256': digest(variants), 'engine_json_sha256': digest(engine),
    },
}
p = pathlib.Path(target); q = p.with_suffix('.json.tmp')
q.write_text(json.dumps(payload, indent=2) + '\n'); os.chmod(q, 0o600); q.replace(p)
PY
done

python3 - "$result_root" "$operation_root/completion.json" <<'PY'
import hashlib, json, os, pathlib, sys
root=pathlib.Path(sys.argv[1]); target=pathlib.Path(sys.argv[2])
runs={p:json.load(open(root/p/'complete.json', encoding='utf-8')) for p in ('baseline','full','axis')}
payload={'schema':'mva1_exomiser_sensitivity_bundle_completion_v1','status':'complete','profiles':runs}
q=target.with_suffix('.tmp'); q.write_text(json.dumps(payload,indent=2)+'\n'); os.chmod(q,0o600); q.replace(target)
PY
