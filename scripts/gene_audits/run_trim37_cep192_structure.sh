#!/usr/bin/env bash
# Native tools and explicit paths; reads inputs, writes only a new output directory.
set -euo pipefail
if [[ ${1:-} == --help ]]; then
    echo 'Required env: AUDIT_OUTPUT (must not exist), AUDIT_SHARED (expanded Ensembl JSON), CORE_BAM, DELLY_BCF, MASKED_BINS, REFERENCE_FASTA. Optional AUDIT_PYTHON, SAMTOOLS, BCFTOOLS, AUDIT_CANDIDATE_TSV.'
    exit 0
fi
umask 077
output=${AUDIT_OUTPUT:?Set a NEW output directory}
shared=${AUDIT_SHARED:?Set saved Ensembl metadata directory}
bam=${CORE_BAM:?Set indexed BAM}
sv=${DELLY_BCF:?Set indexed DELLY BCF}
bins=${MASKED_BINS:?Set retained masked RD/BAF bins TSV}
reference=${REFERENCE_FASTA:?Set matching indexed GRCh38 FASTA}
audit_python=${AUDIT_PYTHON:-python3}
samtools=${SAMTOOLS:-samtools}
bcftools=${BCFTOOLS:-bcftools}
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
[[ ! -e "$output" ]] || { echo 'Refusing existing output directory' >&2; exit 2; }
for input in "$bam" "$sv" "$bins" "$reference" "$shared/ensembl_TRIM37_transcripts.json" "$shared/ensembl_CEP192_transcripts.json"; do [[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 2; }; done
[[ $("$bcftools" query -l "$sv" | wc -l) -eq 1 ]] || { echo "Expected a single-sample DELLY BCF" >&2; exit 2; }
mkdir -p "$output/raw" "$output/analysis"
"$samtools" quickcheck "$bam"
"$audit_python" "$script_dir/map_cep192_published_alleles.py" --metadata "$shared/ensembl_CEP192_transcripts.json" --fasta "$reference" --output "$output/CEP192_published_alleles_mapped.json"
for item in 'TRIM37 chr17:58932633-59156921' 'CEP192 chr18:12941283-13175053'; do
    read -r gene region <<< "$item"
    for mq in 0 20 60; do
        "$samtools" depth -aa -G UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY -q 20 -Q "$mq" -r "$region" -o "$output/raw/$gene.gene.baseq20-mapq$mq.tsv" "$bam"
    done
done
# Query the entire existing BCF, not just records with POS within each gene.
# A giant deletion may contain a gene while both breakpoints lie outside it.
"$bcftools" query -u -f '%ID\t%CHROM\t%POS\t%END\t%SVTYPE\t%FILTER\t%QUAL\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/SRMAPQ\t[%GT\t%GQ\t%FT\t%DR\t%DV\t%RR\t%RV]\t%ALT\t%INFO/CONSENSUS\n' "$sv" -o "$output/raw/delly.all.tsv"
"$audit_python" - "$output/CEP192_published_alleles_mapped.json" > "$output/published-sites.tsv" <<'PY'
import json,sys
for row in json.load(open(sys.argv[1])):
    print(row['pos'], f"{row['chrom']}:{row['pos']}-{row['pos']}", sep='\t')
PY
while IFS=$'\t' read -r pos region; do
    "$samtools" view -h -q 20 -F 3844 -o "$output/raw/published-CEP192-$pos.sam" "$bam" "$region"
done < "$output/published-sites.tsv"
"$audit_python" "$script_dir/published_cep192_support.py" --raw "$output/raw" --output "$output/analysis" --mapped-alleles "$output/CEP192_published_alleles_mapped.json"
extra=()
if [[ -n ${AUDIT_CANDIDATE_TSV:-} ]]; then extra=(--candidate-tsv "$AUDIT_CANDIDATE_TSV"); fi
"$audit_python" "$script_dir/analyze_trim37_cep192_structure.py" --trim-json "$shared/ensembl_TRIM37_transcripts.json" --cep-json "$shared/ensembl_CEP192_transcripts.json" --raw "$output/raw" --bins "$bins" --output "$output/analysis" --report "$output/STRUCTURAL_AND_CALLABILITY_AUDIT.md" "${extra[@]}"
"$samtools" --version > "$output/samtools.version.txt"
"$bcftools" --version > "$output/bcftools.version.txt"
