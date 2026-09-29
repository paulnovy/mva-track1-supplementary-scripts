#!/usr/bin/env bash
# Bounded indexed extraction and offline transcript replay. No network calls.
set -euo pipefail
if [[ ${1:-} == --help ]]; then
    echo 'Required env: AUDIT_OUTPUT (must not exist), AUDIT_SHARED (saved Ensembl JSON), VENDOR_VCF, NORMALIZED_VCF, CORE_VCF. Optional AUDIT_PYTHON, BCFTOOLS.'
    exit 0
fi
umask 077
output=${AUDIT_OUTPUT:?Set a NEW output directory}
shared=${AUDIT_SHARED:?Set saved Ensembl metadata directory}
vendor=${VENDOR_VCF:?Set indexed single-sample vendor VCF}
normalized=${NORMALIZED_VCF:?Set indexed normalized single-sample VCF}
core=${CORE_VCF:?Set indexed single-sample core VCF}
audit_python=${AUDIT_PYTHON:-python3}
bcftools=${BCFTOOLS:-bcftools}
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
[[ ! -e "$output" ]] || { echo 'Refusing existing output directory' >&2; exit 2; }
for input in "$vendor" "$normalized" "$core" "$shared/ensembl_TRIM37_transcripts.json" "$shared/ensembl_CEP192_transcripts.json"; do [[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }; done
mkdir -p "$output"
# No PASS filter: dot, LowQual and all other source FILTER states survive.
# This is the small sequence-variant path; symbolic SV/BND/* alleles require the
# whole-BCF structural path. --regions-overlap 2 uses true changed-variant overlap.
"$bcftools" view --regions-overlap 2 -v snps,indels,mnps -r 17:58982633-59106921,18:12991283-13125053 "$vendor" -Oz -o "$output/vendor-overlap.vcf.gz"
"$bcftools" view --regions-overlap 2 -v snps,indels,mnps -r chr17:58982633-59106921,chr18:12991283-13125053 "$normalized" -Oz -o "$output/normalized-overlap.vcf.gz"
"$bcftools" view --regions-overlap 2 -v snps,indels,mnps -r chr17:58982633-59106921,chr18:12991283-13125053 "$core" -Oz -o "$output/core-overlap-chr.vcf.gz"
for file in vendor-overlap normalized-overlap core-overlap-chr; do "$bcftools" index -t "$output/$file.vcf.gz"; done
"$audit_python" "$script_dir/audit_transcript_panel_rows.py" \
    --source "vendor=$output/vendor-overlap.vcf.gz" --source "normalized=$output/normalized-overlap.vcf.gz" \
    --source "core=$output/core-overlap-chr.vcf.gz" \
    --gene "TRIM37=$shared/ensembl_TRIM37_transcripts.json" --gene "CEP192=$shared/ensembl_CEP192_transcripts.json" \
    --splice-flank 8 --out "$output/all-sources-transcript-audit.tsv"
