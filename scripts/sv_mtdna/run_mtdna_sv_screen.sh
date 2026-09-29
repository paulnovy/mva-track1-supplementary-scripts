#!/usr/bin/env bash
set -euo pipefail

# Independent, immutable-output secondary screen. Research evidence only.
mode=""; bam=""; reference=""; scope=""; out=""
while (($#)); do
  case "$1" in
    --mode) mode=$2; shift 2;; --bam) bam=$2; shift 2;;
    --reference) reference=$2; shift 2;; --scope) scope=$2; shift 2;;
    --out) out=$2; shift 2;; *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done
[[ "$mode" == mtdna || "$mode" == sv ]] || { echo "--mode mtdna|sv required" >&2; exit 2; }
[[ -s "$bam" && -s "$reference" && -s "$scope" && -n "$out" ]] || { echo "bam/reference/scope/out required" >&2; exit 2; }
if [[ -e "$out/complete.json" ]]; then echo "verified complete; skipping: $out" >&2; exit 0; fi
if [[ -e "$out" ]] && find "$out" -mindepth 1 -print -quit | grep -q .; then echo "refusing to overwrite non-empty output: $out" >&2; exit 3; fi
mkdir -p "$out"
idx="${bam}.bai"; [[ -s "$idx" ]] || idx="${bam%.bam}.bai"; [[ -s "$idx" ]] || { echo "BAM index missing" >&2; exit 3; }
command -v samtools >/dev/null || { echo "samtools missing" >&2; exit 4; }
command -v bcftools >/dev/null || { echo "bcftools missing" >&2; exit 4; }

mt=""
while IFS=$'\t' read -r c _; do case "$c" in chrM|MT|M|chrMT) mt=$c; break;; esac; done < <(samtools idxstats "$bam")
[[ -n "$mt" ]] || { echo "mitochondrial contig not present in BAM index" >&2; exit 5; }
region=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["region"])' "$scope")

python3 - "$scope" "$out" "$mt" "$mode" <<'PY'
import json, pathlib, sys
s=json.loads(pathlib.Path(sys.argv[1]).read_text()); out=pathlib.Path(sys.argv[2])
(out/"run.json").write_text(json.dumps({"schema":"mva_secondary_screen_v2","mode":sys.argv[4],"research_only":True,"mitochondrial_contig":sys.argv[3],"sv_scope":s.get("region"),"scope_source":s.get("scope_source")},indent=2)+"\n")
PY

if [[ "$mode" == mtdna ]]; then
  # Keep the default overlap-removal algorithm. samtools -x DISABLES it.
  samtools mpileup -aa -q 30 -Q 25 --ff 3844 -d 100000 -r "$mt" -f "$reference" "$bam" > "$out/mtdna.pileup.tsv"
  python3 - "$out/mtdna.pileup.tsv" "$out" "$mt" <<'PY'
import json,pathlib,sys
p=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]); mt=sys.argv[3]
depths=[]; calls=[]
for line in p.open():
 f=line.rstrip().split("\t")
 if len(f)<6: continue
 pos,ref,reported,bases=int(f[1]),f[2].upper(),int(f[3]),f[4]; counts={b:[0,0] for b in "ACGT"}; i=0; total=0
 while i<len(bases):
  ch=bases[i]
  if ch=="^": i+=2; continue
  if ch=="$": i+=1; continue
  if ch in "+-":
   j=i+1
   while j<len(bases) and bases[j].isdigit(): j+=1
   i=j+int(bases[i+1:j] or 0); continue
  if ch in ".,": total+=1
  elif ch.upper() in counts: counts[ch.upper()][0 if ch.isupper() else 1]+=1; total+=1
  i+=1
 depths.append(total)
 for alt in "ACGT":
  if alt==ref: continue
  fw,rv=counts[alt]; ad=fw+rv
  if total>=20 and ad>=3 and ad/total>=.01:
   calls.append({"contig":mt,"position":pos,"ref":ref,"alt":alt,"depth":total,"alt_depth":ad,"af":round(ad/total,6),"forward":fw,"reverse":rv,"both_strands_2plus":fw>=2 and rv>=2,"quality":"screen_pass" if total>=100 and ad>=5 and fw>=2 and rv>=2 else "low_depth_or_strand"})
(out/"mtdna.depth_summary.json").write_text(json.dumps({"method":"samtools_mpileup","positions_reported":len(depths),"positions_with_depth":sum(x>0 for x in depths),"low_depth_lt20":sum(x<20 for x in depths),"depth_cap":100000,"mapq_min":30,"baseq_min":25,"excluded_flags":3844,"overlap_handling":"samtools_ignore_overlaps","candidate_thresholds":{"depth":20,"alt_reads":3,"af":.01},"candidate_variant_table":"mtdna.pileup.calls.json","no_negative_inference":True,"no_validated_heteroplasmy_sensitivity":True},indent=2)+"\n")
(out/"mtdna.pileup.calls.json").write_text(json.dumps({"calls":calls,"circular_numt_corrected":False,"no_diploid_zygosity":True},indent=2)+"\n")
p.unlink()
PY
  gatk_status="not_run"; filtered_status="not_run"
  if command -v gatk >/dev/null; then
    gatk_status="failed"
    if gatk --java-options "-Xmx4g -XX:ActiveProcessorCount=2" Mutect2 --mitochondria-mode -R "$reference" -I "$bam" -L "$mt" -O "$out/mtdna.native.vcf.gz" >"$out/mtdna.gatk.log" 2>&1; then
      gatk_status="succeeded"; bcftools index -f "$out/mtdna.native.vcf.gz"
      if gatk --java-options "-Xmx4g -XX:ActiveProcessorCount=2" FilterMutectCalls --mitochondria-mode -R "$reference" -V "$out/mtdna.native.vcf.gz" -O "$out/mtdna.filtered.vcf.gz" >"$out/mtdna.filter.log" 2>&1; then
        filtered_status="succeeded"; bcftools index -f "$out/mtdna.filtered.vcf.gz"
      else filtered_status="failed"; fi
    fi
  fi
  python3 - "$out" "$gatk_status" "$filtered_status" <<'PY'
import json,pathlib,sys
o=pathlib.Path(sys.argv[1]); (o/"summary.json").write_text(json.dumps({"status":"complete","method":"gatk_mutect2_plus_pileup" if sys.argv[2]=="succeeded" else "pileup_fallback","gatk_status":sys.argv[2],"filter_mutect_calls_status":sys.argv[3],"native_vcf":"mtdna.native.vcf.gz" if (o/"mtdna.native.vcf.gz").exists() else None,"filtered_vcf":"mtdna.filtered.vcf.gz" if (o/"mtdna.filtered.vcf.gz").exists() else None,"limitations":["not circular-reference corrected","NUMT resources not applied by this screen","AF exploratory; heteroplasmy sensitivity unvalidated","unfiltered GATK records are not screen_pass"],"aggregate_only":True},indent=2)+"\n")
PY
  printf '%s\n' '{"status":"complete","mode":"mtdna","research_only":true}' > "$out/complete.json"
else
  command -v delly >/dev/null || { printf '%s\n' '{"status":"failed","reason":"delly_unavailable"}' > "$out/summary.json"; exit 6; }
  if ! delly call -g "$reference" -o "$out/sv.delly.bcf" "$bam" >"$out/sv.delly.log" 2>&1; then
    printf '%s\n' '{"status":"failed","caller":"delly"}' > "$out/summary.json"; exit 7
  fi
  bcftools index -f "$out/sv.delly.bcf"
  bcftools view -r "$region" -Oz -o "$out/sv.bub1b.vcf.gz" "$out/sv.delly.bcf"
  bcftools index -f "$out/sv.bub1b.vcf.gz"
  python3 - "$out" "$region" <<'PY'
import json,pathlib,sys
o=pathlib.Path(sys.argv[1]); (o/"summary.json").write_text(json.dumps({"status":"complete","caller":"delly","scope":"BUB1B","region":sys.argv[2],"limitations":["short-read SV screen","not exclusionary for complex/repeat-mediated/long insertion events","single-caller research evidence"],"variant_table":"sv.bub1b.vcf.gz","aggregate_only":True},indent=2)+"\n")
PY
  printf '%s\n' '{"status":"complete","mode":"sv","research_only":true}' > "$out/complete.json"
fi
echo "secondary $mode screen complete: $out" >&2
