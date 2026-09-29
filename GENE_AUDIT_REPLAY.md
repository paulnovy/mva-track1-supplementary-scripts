# September 12–14 exploratory and gene-audit replay

These additions supplement [Carbon-DeCoder](https://github.com/paulnovy/Carbon-DeCoder),
not a replacement alignment/calling pipeline. They contain code and public reference
identifiers, **not case genotypes, reads, annotation responses, or recorded results**.
No new API, whole-WGS, patient RNA, or cellular experiment was run to package this
revision. Full historical reproducibility requires authorized access to the pinned
inputs; this repository alone cannot recreate unavailable patient data or API responses.

## Inputs and snapshot boundary

Use Python 3.11+ (`pysam` 0.24.1 observed historically), `samtools`/`bcftools` 1.20,
and optionally matplotlib for the synthetic figure. Programs are selected by
`AUDIT_PYTHON`, `SAMTOOLS`, and `BCFTOOLS`; no root access, named Docker volume or
host-specific path is assumed. A container may wrap these tools externally.

Inputs are read-only, indexed, **single-sample** BAM/VCFs, the matching GRCh38
FASTA, an existing whole-genome DELLY BCF, and saved RD/BAF bins. The latter use
the `rd_baf_followup.py` TSV contract: `chrom,start,end,clean,corrected_rd,
median_abs_baf_minus_half` plus retained mask fields. Vendor VCF uses bare numeric
contigs; normalized/core BAM/VCF/BCF and FASTA use `chr` prefixes. If yours differ,
make separate renamed copies; do not rewrite historical inputs.

Annotation inputs:

- Expanded Ensembl `lookup/id?expand=1;mane=1` JSON including `Transcript`, `Exon`,
  `Translation` and `MANE`, saved as `ensembl_SYMBOL_transcripts.json`.
- September 14 audit: TRIM37 ENSG00000108395.16, MANE ENST00000262294.12 /
  NM_015294.6; CEP192 ENSG00000101639.21, MANE ENST00000506447.5 / NM_032142.4.
  The live `/info/data` request failed; **no release number was established**.
  Do not label this snapshot Ensembl 115 based on an earlier, separate audit.
- Normalized/split ClinVar GRCh38 snapshot dated 2026-09-05, prepared against the
  same reference. A missing exact match is not a benign classification.
- Optional saved exact-allele gnomAD r4 and VEP JSON/stdout/stderr. Save timestamps,
  provider/version if returned, exact request, exit status and SHA-256 alongside
  each response. A current provider response is a new annotation, not a byte-identical
  historical replay. Missing release, FILTER or plugin fields remain missing.

Historical expanded-transcript SHA-256 (public metadata, not patient data):

| gene | SHA-256 |
|---|---|
| TRIM37 | `f0f1475f093b3f2e448129041efc68a8db951e7a1a92b658d2a1bc527107c691` |
| CEP192 | `3b574cdeae88abb9faaba2a4f5df0135e252cd40785b01c1fdb977f9a19f61f9` |

Use retained JSON for offline replay. If intentionally collecting new public
metadata, the bundled clients require explicit network opt-in. For example, in
a **new** directory, not the historical snapshot:

```bash
test ! -e work/new-ensembl-snapshot
mkdir -p work/new-ensembl-snapshot
MVA_ALLOW_NETWORK=1 python scripts/public_api/ensembl_api.py transcripts \
  ENSG00000108395 --assembly GRCh38 \
  --output work/new-ensembl-snapshot/ensembl_TRIM37_transcripts.json
MVA_ALLOW_NETWORK=1 python scripts/public_api/ensembl_api.py transcripts \
  ENSG00000101639 --assembly GRCh38 \
  --output work/new-ensembl-snapshot/ensembl_CEP192_transcripts.json
sha256sum work/new-ensembl-snapshot/*.json > work/new-ensembl-snapshot/SHA256SUMS
```

The wrapper uses human GRCh38 here. Review current Ensembl/gnomAD service terms
before a new query. Do not send sample labels, phenotype or genotype files.
No credentials are needed for these metadata requests.

## TRIM37 / CEP192 small variants, exact ClinVar and published alleles

Set paths using `config/example.env`. `AUDIT_OUTPUT` **must not already exist**;
an interrupted run is preserved and a retry needs a fresh directory.

```bash
AUDIT_OUTPUT="$MVA_DATA_ROOT/replay-new/variants" \
  bash scripts/gene_audits/run_trim37_cep192_variants.sh

python scripts/gene_audits/audit_panel_clinvar_snapshot.py \
  --metadata "$AUDIT_SHARED/ensembl_TRIM37_transcripts.json" \
    "$AUDIT_SHARED/ensembl_CEP192_transcripts.json" \
  --vcf "$NORMALIZED_VCF" --clinvar "$CLINVAR_VCF" \
  --output "$MVA_DATA_ROOT/replay-new/clinvar-exact.json"
```

The runner extracts by **true variant overlap** (`--regions-overlap 2`), not
just POS, retains source FILTER `.` separately from PASS and does not apply a
PASS filter. Sequence SNP/MNP/indel calls feed the small-variant classifier.
Symbolic/BND/spanning-deletion `*` alleles are not meaningful input to the
sequence-span helper; review them in the whole structural-call path instead.
Malformed or mixed unsupported alleles fail explicitly, not as negative findings.
The derived VCFs retain reference/no-call rows. Only called ALT alleles enter the
transcript classification TSV; an ALT merely listed in VCF is not a patient allele.

`audit_transcript_panel_rows.py` reports MANE/canonical and secondary translated
transcripts separately, using each called ALT's minimal changed reference span;
insertions retain both flanking bases. Coding overlap takes precedence. The ±8-bp
splice windows are **intronic flanks of internal exon junctions**, not terminal
pseudo-junctions. This is not deep-intronic/regulatory or genome-wide discovery.

The structural runner below invokes `map_cep192_published_alleles.py` with the
saved MANE transcript, confirms CDS length and genomic REF, and maps the two
alleles in [Guo et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10716027/).
This is reference/literature mapping, not a patient genotype assertion.

## Coverage and all-span structural review

```bash
AUDIT_OUTPUT="$MVA_DATA_ROOT/replay-new/structure" \
  AUDIT_CANDIDATE_TSV="$MVA_DATA_ROOT/replay-new/variants/all-sources-transcript-audit.tsv" \
  bash scripts/gene_audits/run_trim37_cep192_structure.sh
```

This extracts both genes ±50 kb at MAPQ 0/20/60 with BQ≥20, excluding duplicate,
QC-fail, secondary and supplementary alignments. Depth counts read bases, **not
independent templates** (`samtools depth -s` is deliberately not used). Primary
MANE CDS and the union of translated CDS bases are summarized separately.

The whole existing DELLY BCF is queried before filtering by gene-span overlap or
either BND endpoint. This includes large events containing a gene with both
breakpoints outside. All FILTER/GT/evidence fields are preserved; nonreference
genotype spans are compared with usable normalized RD/BAF bins inside and up to
2 Mb outside. The producer derives event IDs, genotypes, statistics and labels
from input files; **historical numerical conclusions are not hard-coded**.
Missing BAF values remain null. Inspect masks and supporting reads before any CN
conclusion. A near-diploid dosage profile is not an inversion test: inversions
remain unvalidated, and direct breakpoint overlap is distinct from enclosure.

The mapped CEP192 sites are separately extracted and counted by independent
template, requiring agreeing eligible mate bases at BQ≥20/MAPQ≥20. Discordant
templates are counted separately. These counts do not establish absence of all
other alleles, low mosaicism, RNA effects or gene dysfunction.

For a **saved-table-only** replay, use a new output directory and invoke
`published_cep192_support.py --raw SAVED_RAW --mapped-alleles MAPPING.json
--output NEW_ANALYSIS`, then `analyze_trim37_cep192_structure.py --help` with
the same saved raw depth/whole-DELLY tables and bins. The analyzer accepts the
optional candidate TSV; it never needs the large BAM for this mode.

## Exact external annotations

For an intentionally selected allele, specify `CHROM-POS-REF-ALT`, not an rsID
that could resolve to a different alternate allele:

```bash
# Set EXACT_VARIANT_ID and VEP_VARIANT to the reviewed allele, not a sample name.
# VEP_VARIANT uses CHROM:POS:REF:ALT, GRCh38 forward-reference alleles.
MVA_ALLOW_NETWORK=1 python scripts/public_api/get_variant_frequency.py \
  --variant_id "$EXACT_VARIANT_ID" --dataset gnomad_r4 --output "$NEW_GNOMAD_JSON"
MVA_ALLOW_NETWORK=1 python scripts/public_api/ensembl_api.py vep \
  "$VEP_VARIANT" --assembly GRCh38 --output "$NEW_VEP_JSON"
```

These are opt-in collection examples, not part of offline replay. Use nonexisting
output filenames and save stderr/exit status. Parse stored JSON without recalling
the APIs; manually confirm its returned exact allele matches the request before
using AC/AN/AF. The gnomAD client does **not request QC FILTER fields**: absent
filters mean unknown, not PASS. Never infer frequency units or joint genotypes
from rsIDs. VEP timeout/no response is not a negative prediction. The September
14 audit had missing full VEP/plugin responses for chr18:13095610 T>C and
chr18:13116433 G>T after bounded attempts; preserve these as unavailable. The
historical rejected rsID-derived alternative must not be reused for a different ALT.

## September 12 nine-gene pilot

The canonical pilot included BUB1, BUB1B, CENATAC, CEP57, MAD1L1, MAD2L1BP,
SLF2, SMC5 and TRIP13. It did not include TRIM37 or CEP192. Expanded saved
Ensembl metadata for each gene can be converted without network calls:

```bash
python scripts/gene_audits/prepare_panel_metadata.py \
  --gene "BUB1=$PILOT_SNAPSHOTS/ensembl_BUB1_transcripts.json" \
  --gene "BUB1B=$PILOT_SNAPSHOTS/ensembl_BUB1B_transcripts.json" \
  --gene "CENATAC=$PILOT_SNAPSHOTS/ensembl_CENATAC_transcripts.json" \
  --gene "CEP57=$PILOT_SNAPSHOTS/ensembl_CEP57_transcripts.json" \
  --gene "MAD1L1=$PILOT_SNAPSHOTS/ensembl_MAD1L1_transcripts.json" \
  --gene "MAD2L1BP=$PILOT_SNAPSHOTS/ensembl_MAD2L1BP_transcripts.json" \
  --gene "SLF2=$PILOT_SNAPSHOTS/ensembl_SLF2_transcripts.json" \
  --gene "SMC5=$PILOT_SNAPSHOTS/ensembl_SMC5_transcripts.json" \
  --gene "TRIP13=$PILOT_SNAPSHOTS/ensembl_TRIP13_transcripts.json" \
  --out "$MVA_DATA_ROOT/replay-new/pilot-metadata"
python scripts/gene_audits/extract_panel_vcf_rows.py --vcf "$NORMALIZED_VCF" \
  --coordinates "$MVA_DATA_ROOT/replay-new/pilot-metadata/panel_gene_coordinates.tsv" \
  --out "$MVA_DATA_ROOT/replay-new/pilot-rows.tsv"
python scripts/gene_audits/classify_canonical_panel_rows.py \
  --records "$MVA_DATA_ROOT/replay-new/pilot-rows.tsv" \
  --coordinates "$MVA_DATA_ROOT/replay-new/pilot-metadata/panel_gene_coordinates.tsv" \
  --structures "$MVA_DATA_ROOT/replay-new/pilot-metadata" \
  --out "$MVA_DATA_ROOT/replay-new/pilot-classified.tsv"
```

If only the original `panel_gene_coordinates.tsv` and individual
`ensembl_SYMBOL_canonical_structure.json` snapshots are available, use those
directly in the last two commands; metadata conversion is unnecessary. Fresh
canonical selection may differ from September 12. Historical extraction
selected anchors within gene intervals; this revision deliberately extends the
initial selection to true changed-span overlap so a deletion starting outside
a gene is not missed. Record this method change, not a claim of byte-identical
TSVs. Interpretation remains canonical coding/intronic ±2 only, not full gene
exclusion. Exomiser “without BUB1B” is a filtered display, not reranking evidence.

## N1002K sequence-control construction and synthetic mixtures

Generate local reference mRNA/CDS/protein FASTA first using
`scripts/alphagenome/reconstruct_bub1b_transcripts.py` as documented in
`RECIPES.md`. The control helper expects the FASTA record
`NM_001211.6_reference_GRCh38_local`, plus a saved UniProt O60566 full entry JSON
(raw entry or a `results` array wrapper). The public entry can be obtained from
`https://rest.uniprot.org/uniprotkb/O60566.json` in a separately authorized request;
record the entry/sequence versions, date and hash. It is not bundled.

```bash
python scripts/exploratory/map_bubr1_experimental_controls.py \
  --protein-fasta "$BUB1B_PROTEIN_FASTA" --cds-fasta "$BUB1B_CDS_FASTA" \
  --uniprot-json "$UNIPROT_O60566_JSON" --out "$MVA_DATA_ROOT/replay-new/controls"
python scripts/exploratory/bulk_mixture_identifiability.py \
  --out "$MVA_DATA_ROOT/replay-new/mixtures"
python scripts/exploratory/plot_bulk_mixture_examples.py \
  --table "$MVA_DATA_ROOT/replay-new/mixtures/bulk_mixture_expectations.tsv" \
  --out "$MVA_DATA_ROOT/replay-new/mixture-figures"
```

Sequence controls are WT, N1002K, L1012P and Q921H reference derivatives, **not
expression-ready constructs, patient RNA, or functional experiments**. Exact
c.3006T>G sequence and terminal-stop checks are retained. No stability prediction,
phase or biological rescue is inferred. Protein and RNA/NMD experiments are
distinct; L1012P proximity does not establish a mechanism for N1002K.

The mixture script is entirely synthetic and requires no case file. Its optional
`--bins` only records file-structure counts/hash, not patient RD/BAF values or a
detection threshold. Balanced cancellation is a mathematical possibility, not
the patient's burden, likelihood, sensitivity or causal timing.

## Verification and honest limits

`python -m unittest discover -s tests -v` runs focused arithmetic, strand/allele
geometry, called-ALT/FILTER, spanning-SV/bin and overwrite checks plus the existing
small bundle suite. Source SHA-256 relationships are in `SOURCE_PROVENANCE.tsv`.
Pure-Python stages can be compiled/help-checked without running genomic analyses;
plotting needs matplotlib. Hashes establish source identity, not independent
biological validation. No live annotation or full biological replay is claimed
for this packaging revision. Historical research reports remain separate from
the future public code repository.
