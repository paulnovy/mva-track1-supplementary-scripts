#!/bin/bash
# Run only after the queue state says ready. Inputs stay read-only; output must be new.
set -euo pipefail
umask 077

if [[ $# -ne 5 ]]; then
  echo "usage: $0 PHASED_WGS_BCF EXCLUDED_VARIANTS_BCF CNPS_BED OUTPUT_DIR SAMPLE" >&2
  exit 64
fi
input=$1
sites=$2
cnps=$3
out=$4
sample=$5

for path in "$input" "$input.csi" "$sites" "$sites.csi" "$cnps"; do
  [[ -r "$path" ]] || { echo "missing prerequisite: $path" >&2; exit 65; }
done
input=$(realpath -e -- "$input")
sites=$(realpath -e -- "$sites")
cnps=$(realpath -e -- "$cnps")
out=$(realpath -m -- "$out")
[[ ! -e "$out" ]] || { echo "refusing existing output: $out" >&2; exit 73; }
mkdir -m 700 "$out"

image=${MOCHA_IMAGE:-mva-mocha:bcftools-1.20-95686b7}
exec docker run --rm --network none --cpus=1 --memory=3g \
  --user "$(id -u):$(id -g)" \
  -v "$(dirname "$input"):/input:ro" \
  -v "$(dirname "$sites"):/excluded:ro" \
  -v "$(dirname "$cnps"):/cnps:ro" \
  -v "$(dirname "$out"):/output:rw" \
  "$image" \
  +mocha -g GRCh38 -s "$sample" -v "^/excluded/$(basename "$sites")" -p /cnps/"$(basename "$cnps")" \
  -c /output/"$(basename "$out")"/calls.tsv -z /output/"$(basename "$out")"/stats.tsv \
  -o /output/"$(basename "$out")"/annotated.bcf -Ob --write-index /input/"$(basename "$input")"
