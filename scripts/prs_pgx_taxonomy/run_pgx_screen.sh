#!/usr/bin/env bash
set -euo pipefail
umask 077

# Bounded, offline PharmCAT screen. Inputs are read-only; OUT must be new.
# No --missing-to-ref is ever used. Research-only; no drug recommendations.
bam=""; reference=""; positions=""; preprocessor=""; jar="/opt/pharmcat/pharmcat.jar"; out=""; sample=""
while (($#)); do
  case "$1" in
    --bam) bam=$2; shift 2;; --reference) reference=$2; shift 2;;
    --positions) positions=$2; shift 2;; --preprocessor) preprocessor=$2; shift 2;;
    --jar) jar=$2; shift 2;; --out) out=$2; shift 2;; --sample) sample=$2; shift 2;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done
[[ -s "$bam" && -s "$reference" && -s "$positions" && -x "$preprocessor" && -s "$jar" && -n "$out" ]] || { echo "bam/reference/positions/preprocessor/jar/out required" >&2; exit 2; }
idx="${bam}.bai"; [[ -s "$idx" ]] || idx="${bam%.bam}.bai"; [[ -s "$idx" ]] || { echo "BAM index missing" >&2; exit 3; }
[[ ! -e "$out" || -z "$(find "$out" -mindepth 1 -print -quit)" ]] || { echo "refusing non-empty output" >&2; exit 3; }
mkdir -p "$out" "$out/caller" "$out/preprocessor" "$out/pharmcat"
command -v gatk >/dev/null || { echo "GATK missing" >&2; exit 4; }
command -v bcftools >/dev/null || { echo "bcftools missing" >&2; exit 4; }
command -v java >/dev/null || { echo "Java missing" >&2; exit 4; }
export PYTHONPATH="${PHARMCAT_PYTHON_OVERLAY:-$(dirname "$preprocessor")/../python-overlay}:$(dirname "$preprocessor"):${PYTHONPATH:-}"
bgzip_tool="${PHARMCAT_BGZIP:-$(dirname "$preprocessor")/../python-overlay/bgzip}"
[[ -x "$bgzip_tool" ]] || { echo "bgzip compatibility tool missing" >&2; exit 4; }
resource_dir=$(dirname "$positions")
manifest="$resource_dir/manifest.json"
[[ -s "$manifest" ]] || { echo "PharmCAT resource manifest missing" >&2; exit 4; }
python3 - "$manifest" "$jar" "$positions" "$resource_dir/pharmcat_positions.uniallelic.vcf.bgz" <<'PY'
import hashlib,json,pathlib,sys
m=json.loads(pathlib.Path(sys.argv[1]).read_text())
checks=[(sys.argv[2],m["jar"]["sha256"]),(sys.argv[3],m["positions_vcf"]["sha256"]),(sys.argv[4],m["uniallelic_positions_vcf"]["sha256"])]
for name,want in checks:
 h=hashlib.sha256();
 with open(name,"rb") as f:
  for chunk in iter(lambda:f.read(1024*1024),b""): h.update(chunk)
 if h.hexdigest()!=want: raise SystemExit(f"resource checksum mismatch: {name}")
PY

mapfile -t bam_samples < <(samtools view -H "$bam" | awk -F '\t' '$1=="@RG" {for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i}}' | sort -u)
if [[ -z "$sample" ]]; then
  [[ ${#bam_samples[@]} -eq 1 ]] || { echo "sample label required when BAM has zero or multiple SM tags" >&2; exit 2; }
  sample=${bam_samples[0]}
else
  [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "unsafe sample label" >&2; exit 2; }
  if [[ ${#bam_samples[@]} -eq 1 && "${bam_samples[0]}" != "$sample" ]]; then echo "supplied sample does not match BAM SM" >&2; exit 2; fi
fi
alignment_for_call="$bam"
if [[ ${#bam_samples[@]} -eq 0 ]]; then
  echo "BAM has no SM read-group; refusing whole-BAM rewrite for PGx" >&2
  exit 2
fi

# Verify that every PharmCAT position is represented on the supplied indexed
# reference and record the known Cockpit reference provenance; do not silently
# liftover or download a 0.9-GB alternate FASTA.
python3 - "$positions" "$reference" "$out/reference_check.json" <<'PY'
import json,pathlib,sys
import pysam
vf=pysam.VariantFile(sys.argv[1]); fa=pysam.FastaFile(sys.argv[2]); seen=[]; bad=[]
for rec in vf:
 chrom=str(rec.contig); pos=int(rec.pos); ref=str(rec.ref); seen.append(chrom)
 try: observed=fa.fetch(chrom,pos-1,pos-1+len(ref)).upper()
 except (KeyError,ValueError): observed=""
 if observed != ref.upper(): bad.append({"contig":chrom,"position":pos,"expected_ref":ref,"reference_ref":observed})
if bad: raise SystemExit(f"positions/reference REF mismatch: {bad[:3]}")
pathlib.Path(sys.argv[3]).write_text(json.dumps({"status":"compatible_by_contig_length_and_ref","positions":len(seen),"distinct_contigs":sorted(set(seen)),"exact_vendor_reference":False,"sequence_identity_not_rechecked":True,"no_liftover":True},indent=2)+"\n")
PY

# Force-call only the PharmCAT allele resource positions; output includes only
# evidence emitted by the caller. Missing/unreadable sites remain no-call.
positions_local="$out/caller/pharmcat_positions.vcf.bgz"
cp "$positions" "$positions_local"
gatk --java-options "-Xmx1g -XX:ActiveProcessorCount=2" IndexFeatureFile -I "$positions_local" >"$out/caller/position_index.log" 2>&1 || bcftools index -f "$positions_local"
gatk --java-options "-Xmx4g -XX:ActiveProcessorCount=2" HaplotypeCaller --native-pair-hmm-threads 2 -R "$reference" -I "$alignment_for_call" -L "$positions_local" --alleles "$positions_local" --genotype-filtered-alleles true -ERC GVCF -O "$out/caller/${sample}.targeted.g.vcf.gz" >"$out/caller/haplotypecaller.log" 2>&1
gatk --java-options "-Xmx4g -XX:ActiveProcessorCount=2" GenotypeGVCFs -R "$reference" -V "$out/caller/${sample}.targeted.g.vcf.gz" -L "$positions_local" --include-non-variant-sites true -O "$out/caller/${sample}.targeted.vcf.gz" >"$out/caller/genotypegvcfs.log" 2>&1
bcftools index -f "$out/caller/${sample}.targeted.vcf.gz"

# PharmCAT does not interpret FILTER/QUAL, so remove weak or missing genotypes
# before preprocessing. A site failing QC stays absent/no-call, never 0/0.
python3 - "$out/caller/${sample}.targeted.vcf.gz" "$out/caller/${sample}.qc.vcf.gz" "$sample" <<'PY'
import pysam,sys
src=pysam.VariantFile(sys.argv[1]); out=pysam.VariantFile(sys.argv[2],"wz",header=src.header); sample=sys.argv[3]
if sample not in src.header.samples: raise SystemExit("caller sample mismatch")
counts={"input_records":0,"retained_records":0,"retained_homref":0,"retained_nonref":0,"rejected_missing":0,"rejected_depth":0,"rejected_quality":0}
for rec in src:
 counts["input_records"]+=1; call=rec.samples[sample]; gt=call.get("GT"); dp=call.get("DP"); gq=call.get("GQ") if call.get("GQ") is not None else (call.get("RGQ") if gt and all(x==0 for x in gt) else None)
 if gt is None or any(x is None for x in gt): counts["rejected_missing"]+=1; continue
 if dp is None or int(dp)<10: counts["rejected_depth"]+=1; continue
 if gq is None or float(gq)<20: counts["rejected_quality"]+=1; continue
 counts["retained_records"]+=1
 if all(x==0 for x in gt): counts["retained_homref"]+=1
 else: counts["retained_nonref"]+=1
 out.write(rec)
out.close(); src.close()
import json,pathlib
pathlib.Path(str(sys.argv[2])+".qc.json").write_text(json.dumps({**counts,"quality_fields":["GQ","RGQ for homREF only"],"thresholds":{"DP":10,"GQ_or_RGQ":20},"no_callable_positions":counts["retained_records"]==0},indent=2)+"\n")
PY
bcftools index -f "$out/caller/${sample}.qc.vcf.gz"

cp "$positions" "$out/preprocessor/pharmcat_positions.vcf.bgz"
[[ -s "$positions.csi" ]] && cp "$positions.csi" "$out/preprocessor/pharmcat_positions.vcf.bgz.csi" || true
[[ -s "$resource_dir/pharmcat_positions.uniallelic.vcf.bgz" ]] || { echo "uniallelic PharmCAT positions resource missing" >&2; exit 5; }
cp "$resource_dir/pharmcat_positions.uniallelic.vcf.bgz" "$out/preprocessor/pharmcat_positions.uniallelic.vcf.bgz"
[[ -s "$resource_dir/pharmcat_positions.uniallelic.vcf.bgz.csi" ]] && cp "$resource_dir/pharmcat_positions.uniallelic.vcf.bgz.csi" "$out/preprocessor/pharmcat_positions.uniallelic.vcf.bgz.csi" || true
python3 - "$preprocessor" "$out/caller/${sample}.qc.vcf.gz" "$out/preprocessor/pharmcat_positions.vcf.bgz" "$reference" "$sample" "$out/preprocessor" "$bgzip_tool" <<'PY'
import subprocess,sys
cmd=[sys.argv[1],"-vcf",sys.argv[2],"-s",sys.argv[5],"-o",sys.argv[6],"-bf",sys.argv[5],"-ss","-refVcf",sys.argv[3],"-refFna",sys.argv[4],"-G","-bgzip",sys.argv[7]]
subprocess.run(cmd,check=True)
PY
mapfile -t ready < <(find "$out/preprocessor" -maxdepth 1 -type f \( -name '*.preprocessed.vcf' -o -name '*.preprocessed.vcf.gz' -o -name '*.preprocessed.vcf.bgz' \) | sort)
[[ ${#ready[@]} -eq 1 ]] || { echo "expected one preprocessed VCF" >&2; exit 6; }
mapfile -t missing < <(find "$out/preprocessor" -maxdepth 1 -type f \( -name '*.missing_pgx_var.vcf' -o -name '*.missing_pgx_var.vcf.gz' -o -name '*.missing_pgx_var.vcf.bgz' \) | sort)
[[ ${#missing[@]} -eq 1 ]] || { echo "missing-position report absent" >&2; exit 6; }

compatible="$out/pharmcat/${sample}.pharmcat-compatible.vcf"
python3 "$(dirname "$(readlink -f "$0")")/sanitize_vcf_qual.py" \
  --input "${ready[0]}" --output "$compatible" \
  --provenance "$out/pharmcat/${sample}.qual-provenance.json"
compatible=$(readlink -f "$compatible")
jar=$(readlink -f "$jar")
# PharmCAT's logger writes pharmcat.log in its working directory. Keep that
# writable and separate from captured stdout/stderr; preserve upstream inputs.
(
  cd "$out/pharmcat"
  java -Xmx4g -XX:ActiveProcessorCount=2 -jar "$jar" -vcf "$compatible" -o . -bf "$sample" -reporterJson -reporterCallsOnlyTsv
) >"$out/pharmcat/runner.log" 2>&1
calls="$out/pharmcat/${sample}.report.tsv"; phenotype="$out/pharmcat/${sample}.phenotype.json"
[[ -s "$calls" && -s "$phenotype" ]] || { echo "PharmCAT calls/phenotype output missing" >&2; exit 7; }

python3 - "$calls" "${missing[0]}" "$phenotype" "$out/caller/${sample}.qc.vcf.gz.qc.json" "$out/summary.json" <<'PY'
import csv,json,pathlib,sys,gzip
calls=pathlib.Path(sys.argv[1]); missing=pathlib.Path(sys.argv[2]); phenotype=pathlib.Path(sys.argv[3]); qc=json.loads(pathlib.Path(sys.argv[4]).read_text()); genes={}; missing_genes={}
with calls.open(newline="") as h:
 lines=h.readlines(); header_idx=next((i for i,x in enumerate(lines) if x.startswith("Gene\t")),None)
 if header_idx is None: raise SystemExit("PharmCAT calls TSV header missing")
 for row in csv.DictReader(lines[header_idx:],delimiter="\t"):
  g=(row.get("Gene") or row.get("gene") or "UNKNOWN").strip(); genes.setdefault(g,{"rows":0,"diplotypes":set(),"phenotypes":set()}); genes[g]["rows"]+=1
  for k,d in (("diplotypes","Source Diplotype"),("phenotypes","Phenotype")):
   if row.get(d): genes[g][k].add(row[d])
mh=gzip.open(missing,"rt") if missing.suffix in {".gz",".bgz"} else missing.open()
for line in mh:
 if "PX=" in line:
  g=line.split("PX=",1)[1].split(";",1)[0].split("\t",1)[0]; missing_genes[g]=missing_genes.get(g,0)+1
mh.close()
summary={"status":"complete","engine":"PharmCAT","engine_version":"3.3.0","aggregate_only":True,"qc":qc,"genes":{},"missing_positions_by_gene":missing_genes,"raw_calls":"pharmcat/"+calls.name,"raw_phenotype":"pharmcat/"+phenotype.name,"limitations":["research-only; no recommendations or doses","unphased genotype ambiguity retained; review competing diplotypes and per-gene missingness","missing positions remain no-call, never homREF","CYP2D6 not called without Cyrius/StellarPGx CNV/hybrid evidence","VCF QUAL/FILTER are not interpreted by PharmCAT; caller QC is upstream"]}
for g,v in genes.items():
 if g=="CYP2D6": continue
 summary["genes"][g]={"rows":v["rows"],"diplotype_count":len(v["diplotypes"]),"phenotype_count":len(v["phenotypes"])}
summary["cyp2d6"]={"status":"not_called","reason":"specialized CNV/hybrid caller required"}
pathlib.Path(sys.argv[5]).write_text(json.dumps(summary,indent=2)+"\n")
PY
printf '%s\n' '{"status":"complete","module":"pgx","research_only":true}' > "$out/complete.json"
echo "PGx screen complete: $out" >&2
