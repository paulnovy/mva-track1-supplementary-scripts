#!/usr/bin/env bash
set -euo pipefail

result_dir=${1:?result directory required}
operation_dir=${2:?operation directory required}
vcf=${3:?VCF required}
adapter_root=${4:?wgs-cockpit root required}
umask 077
mkdir -p "$result_dir" "$operation_dir"
chmod 700 "$result_dir" "$operation_dir"

exec 9>"$operation_dir/run.lock"
flock -n 9 || exit 0
[[ ! -e "$result_dir/complete.json" ]] || exit 0
[[ -s "$result_dir/phenopacket.json" && -s "$result_dir/family.ped" && -s "$result_dir/resources/resource-manifest.json" ]] || exit 20

free_kb=$(df -Pk "$result_dir" | awk 'NR==2 {print $4}')
(( free_kb >= 60 * 1024 * 1024 )) || exit 21

manifest="$result_dir/resources/resource-manifest.json"
jar=$(python3 - "$manifest" <<'PY'
import json, pathlib, sys
m=json.load(open(sys.argv[1], encoding='utf-8'))
root=pathlib.Path(sys.argv[1]).parent
print(root / next(x['path'] for x in m['files'] if x['id']=='exomiser_jar'))
PY
)
jar_sha=$(sha256sum "$jar" | awk '{print $1}')
profile_sha=$(sha256sum "$result_dir/phenopacket.json" | awk '{print $1}')
java="$result_dir/resources/jre21/bin/java"
[[ -x "$java" ]]
work="$result_dir/work"
mkdir -p "$work"
(
  cd "$work"
  export PATH="$(dirname "$java"):$PATH"
  export PYTHONPATH="$adapter_root/apps/api"
  export WGS_EXOMISER_JAR="$jar"
  export WGS_EXOMISER_JAR_SHA256="$jar_sha"
  export WGS_EXOMISER_VERSION=15.1.0
  export WGS_RARE_DISEASE_RESOURCE_MANIFEST="$manifest"
  export WGS_PHENOPACKET_PROFILE_SHA256="$profile_sha"
  export WGS_EXOMISER_JAVA_OPTS="-Xmx10g -XX:+UseG1GC -Dexomiser.data-directory=$result_dir/resources/data -Dexomiser.hg38.data-version=2512 -Dexomiser.phenotype.data-version=2512"
  "$adapter_root/pipelines/nextflow/scripts/run_rare_disease_stage.sh" MVA1 "$result_dir/phenopacket.json" "$vcf" "$result_dir/family.ped" hg38
) >"$operation_dir/exomiser.stdout.log" 2>"$operation_dir/exomiser.stderr.log"

python3 - "$result_dir" <<'PY'
import hashlib, json, os, pathlib, sys
root=pathlib.Path(sys.argv[1])
payload={"schema":"mva1_hpo_exomiser_completion_v1","status":"complete","research_only":True,"ranking_is_not_diagnosis":True}
p=root/'complete.json.tmp'; p.write_text(json.dumps(payload, indent=2)+'\n'); os.chmod(p,0o600); p.replace(root/'complete.json')
PY
