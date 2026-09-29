# Population-panel, RGC, and long-range feasibility replay

This note captures the **public query definitions and preparation parameters**, not the historical observations or any patient-level data. It distinguishes executable public-resource queries from manual portal inspection. A missing record is only a release/callset-specific non-return; it must never be reported as absence from a cohort or population.

## Exact alleles and coordinate systems

| Label | GRCh38, 1-based | GRCh37/hs37d5, 1-based | dbSNP |
|---|---|---|---|
| BUB1B Leu737Ter | `chr15:40209701:T>G` | `chr15:40501902:T>G` | `rs759242053` |
| BUB1B Asn1002Lys | `chr15:40220612:T>G` | `chr15:40512813:T>G` | `rs2542593804` |

Always verify chromosome naming, assembly, position, REF, and exact ALT. Do not substitute another nucleotide allele producing the same amino-acid consequence. For zero-based APIs the GRCh38 positions are `40209700` and `40220611`.

## Executable public panel queries

The historical HGDP+1KG Shapeit5 object was:

```text
https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr15.filtered.SNV_INDEL.phased.shapeit5.bcf
```

With `bcftools`/htslib built for remote HTTP range requests, the bounded replay is:

```bash
panel='https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chr15.filtered.SNV_INDEL.phased.shapeit5.bcf'
bcftools view -H -r chr15:40209701-40209701 "$panel"
bcftools view -H -r chr15:40220612-40220612 "$panel"
bcftools view -h "$panel" | bcftools query -l | wc -l
```

The HPRC Release 2 GRCh38 graph-derived VCF and index were:

```text
https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.wave.vcf.gz
https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.wave.vcf.gz.tbi
```

Replay the exact sites and a bounded containing interval:

```bash
hprc='https://human-pangenomics.s3.amazonaws.com/pangenomes/freeze/release2/minigraph-cactus/hprc-v2.0-mc-grch38.wave.vcf.gz'
bcftools view -H -r chr15:40209701-40209701 "$hprc"
bcftools view -H -r chr15:40220612-40220612 "$hprc"
bcftools view -H -r chr15:40209701-40220612 "$hprc" > hprc.bounded.public.tsv
```

Inspect returned REF/ALT, FILTER, INFO/AC/AN/NS, FORMAT/GT, and phasing provenance. A phased delimiter describes the panel record; it does not phase another sample.

## SGDP 2021 discovery boundary

The public resource documentation and linked directory were:

- `https://reichdata.hms.harvard.edu/pub/datasets/sgdp/`
- `https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/`

The preserved bounded discovery command was equivalent to:

```bash
curl -fL --silent --show-error --connect-timeout 10 --max-time 25 \
  --max-filesize 10485760 \
  'https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021/' \
  -o sgdp-phased-directory.html
```

Only after a current filename and index are discovered may `bcftools view -r` be used with the **GRCh37/hs37d5** coordinates above. Failure to list the directory or open an index means *unqueried*, not allele-negative.

## dbSNP/ALFA aggregate queries

The bundled Apache-2.0 dbSNP client supports exact public RefSNP resolution. Network execution is deliberately opt-in:

```bash
MVA_ALLOW_NETWORK=1 python scripts/public_api/dbsnp_cli.py get-variant rs759242053 --output leu737ter.refsnp.json
MVA_ALLOW_NETWORK=1 python scripts/public_api/dbsnp_cli.py get-variant rs2542593804 --output n1002k.refsnp.json
MVA_ALLOW_NETWORK=1 python scripts/public_api/dbsnp_cli.py resolve-variant 15 40209701 T G --output leu737ter.resolve.json
MVA_ALLOW_NETWORK=1 python scripts/public_api/dbsnp_cli.py resolve-variant 15 40220612 T G --output n1002k.resolve.json
```

RefSNP frequency annotations are aggregate and release-specific. Empty ALFA/TOPMed/PAGE annotations are not proof that the allele is absent from all contributing participants. Do not add overlapping database denominators or allele counts as independent samples.

## RGC-ME: manual portal inspection only

The preserved source pages were:

- `https://rgc-research.regeneron.com/me/variant/15:40209701:T:G`
- `https://rgc-research.regeneron.com/me/variant/15:40220612:T:G`
- `https://rgc-research.regeneron.com/me/gene/BUB1B`

The historical evidence was obtained through the public portal's own page context. No stable, documented public command-line API contract or reusable request-header specification was preserved, so this repository does **not** claim an algorithmic RGC replay. Reinspect the pages manually and record portal/release metadata, exact allele identity, FILTER, AC/AN/AF, transcript, and metric version. Aggregate RGC values cannot establish pair co-occurrence or phase, and RGC/UKB/gnomAD contributions may overlap.

## Other manual or unavailable portals

UK Biobank AFB, BRAVO, FinnGen, All of Us, Genomics England, Geno2MP, jMorp, SweGen, regional catalogues, and other large-cohort routes were portal/access checks rather than preserved executable algorithms. Authorization failures, timeouts, empty portal responses, or controlled-access boundaries are **unavailable/unqueried**, not negative findings. This bundle makes no claim to reproduce them from undocumented internal endpoints.

## Long-range assay feasibility preparation

The computational preflight used the inclusive GRCh38 target span `chr15:40209701-40220612` (10,912 bases; target separation 10,911 bp) and a reference-only design window with 1,000 bp flanks: `chr15:40208701-40221612` (12,912 bases). Prepare the reference interval without patient reads:

```bash
printf 'chr15\t40208700\t40221612\n' > BUB1B_design_window.bed
samtools faidx GRCh38.fa chr15:40208701-40221612 > BUB1B_GRCh38_design_window.reference.fa
python3 - <<'PY'
from pathlib import Path
seq=''.join(x.strip() for x in Path('BUB1B_GRCh38_design_window.reference.fa').read_text().splitlines() if not x.startswith('>'))
assert len(seq) == 12912
assert seq[40209701-40208701] == 'T'
assert seq[40220612-40208701] == 'T'
print({'window_bases': len(seq), 'target_separation_bp': 40220612-40209701})
PY
```

This produces a reference design window, not primers, an amplicon, patient haplotypes, or an executable laboratory result. The cited ONT amplicon workflow (`https://github.com/j-jamshidi/ONT_amp_phase`) requires new barcoded ONT amplicon reads and a sample sheet. Existing short reads must not be relabelled as ONT input. Any future assay requires retained DNA, validated primer specificity, allele-dropout and chimera controls, and reads spanning both exact bases; no success probability follows from feasibility alone.
