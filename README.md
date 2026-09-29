# MVA Track 1 supplementary scripts

**[Carbon-DeCoder](https://github.com/paulnovy/Carbon-DeCoder) provided the automated foundation for the Track 1 analyses**, including primary whole-genome processing and workflow orchestration. This repository contains the supplementary, case-focused analysis and audit scripts that build on those automated outputs. The bundle contains methods, not patient data or biological results. It is research software, not a diagnostic pipeline.

## Start here

1. Read [METHODS.md](METHODS.md) for analysis scope, recorded resources, data boundaries and limitations.
2. Use [REGISTRY.tsv](REGISTRY.tsv) to select an analysis family and its input/dependency contract.
3. Follow [RECIPES.md](RECIPES.md) and the linked ordered replay documents. Supply authorized private inputs and immutable resource snapshots separately.
4. Use [config/example.env](config/example.env) for filesystem paths only. Phenotype curation and optional candidate interpretations are private inputs, not bundled observations.

This supplies the documented scripts and configuration contracts requested by the
[Track 1 specification](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/c9b4a7e6247574124cd5bbfd249af9af6f1bbac3/tabs/submit_track1.py).
A full methods/results report is still submitted separately. Public submission
CSVs and individual-level results are optional and deliberately **not** included.
No FASTQ/BAM/VCF, clinical document, credential, or private report archive is
published. Source contains bounded public/reference target coordinates; they are
not a comprehensive individual callset. Do not copy private inputs into this repo.

## Published supplement — 29 September 2026

The code bundle now includes the September 12 nine-gene canonical pilot,
reference-derived BUBR1 experimental controls and synthetic bulk-mixture model,
plus the September 14 TRIM37/CEP192 transcript, exact-ClinVar, published-allele,
callability and all-span structural audits. Follow [GENE_AUDIT_REPLAY.md](GENE_AUDIT_REPLAY.md)
for ordered commands and snapshot requirements. This public release adapts the reviewed local sources for independent
reproduction. It does not submit an entry to the competition.

New audit runners require a fresh output directory and explicit input paths.
Structural summaries calculate from supplied event/bin tables; they do not embed
the historical sample's variants or conclusions. Missing VEP/QC fields remain
unavailable rather than negative/PASS. The original September 11 bundle and
historical data are preserved.

## Screening-integrity revision

See [PARTIAL_GT_REPLAY.md](PARTIAL_GT_REPLAY.md) and
[MOCHA_WRAPPER_REPLAY.md](MOCHA_WRAPPER_REPLAY.md) for retained-data audits.
ClinVar aggregation now keeps historical complete-GT counts and adds explicit
called-ALT/partial-called-ALT classification strata. P/LP triage retains a partial
call only when an ALT is explicitly present and marks it for original genotype
review; it does not manufacture a complete genotype, phase, or likelihood.
Original multiallelic allele identities and GT/AD/PL remain authoritative.
This is a prospective selection change, not a claim that historical outputs
already used it. The MoChA wrapper's excluded-site argument now includes `^`,
matching the historical command. No caller rerun is required for these audits.

## Primary workflow dependency

Alignment, coverage, small-variant calling, standard phasing, SV calling, rare-disease analysis, mtDNA, PRS, PGx, and taxonomy orchestration are provided by [Carbon-DeCoder](https://github.com/paulnovy/Carbon-DeCoder), pinned here to commit `cd03b3bdb6a29847916429d072f27de94c88065c`. They are not duplicated. This repository adds the audits, targeted follow-ups, offline interpretation, and case-specific reconstruction performed for Track 1.

## Install

The scripts target Python 3.11+; the pinned publication-verification environment uses Python 3.14.
Create a project-local environment. `requirements.txt` records the tested Python
dependencies, not the historical analysis environment:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
# Optional model-prediction module (not needed for offline tests):
# python -m pip install alphagenome==0.9.0
# API wrapper scripts invoked through uv also require the uv CLI.
```

External command-line tools are module-specific: `bcftools`, `samtools`, `bgzip/tabix`, Java 21 + Exomiser, PharmCAT/GATK, Kraken2, or AlphaGenome 0.9.0. See `REGISTRY.tsv` for exact contracts. Reference files and tool databases are deliberately not bundled.

Build the pinned MoChA image when that module is selected:

```bash
docker build -t mva-mocha:bcftools-1.20-95686b7 -f docker/Dockerfile.mocha .
```

The build downloads public bcftools/MoChA sources; the subsequent caller runs with networking disabled.

## Configuration boundary

Copy `config/example.env` outside version control, set local paths, and source it. No credential belongs in this file. Most modern scripts expose `--help`; preserved report-stage scripts use `MVA_DATA_ROOT` and the documented directory contract. Use a fresh output directory and retain input checksums.

```bash
cp config/example.env config/local.env
${EDITOR:-vi} config/local.env
. config/local.env
python scripts/qc/audit_vcf_reference.py --help
python scripts/copy_number/rd_baf_followup.py --help
python scripts/phase/prepare_regional_phasing.py --help
```

## Networked stages and offline replay

All public API work is opt-in. Preparation, saved-response parsing, scoring, plotting, and verification can run offline. AlphaGenome prediction modes and exact-allele annotation wrappers require `MVA_ALLOW_NETWORK=1` after reviewing the bounded request. Standalone AlphaGenome additionally requires `MVA_ALPHAGENOME_TRANSPORT=native` and an externally provisioned environment credential; gateway sentinels are rejected. See [METHODS.md](METHODS.md#offline-and-online-execution-boundary) for the SDK retry/timeout boundary. Credentials never belong in command arguments or committed configuration. `query_alphagenome_zipstored.py` and the analysis modules replay existing local tables without network access. Atlas ZIP replay requires local `tabix` and temporary disk equal to the selected compressed member.

## Reproducing families

`REGISTRY.tsv` maps every analysis family to its implementation or to the pinned Carbon-DeCoder stage. Start with the row for the desired family, satisfy its input contract, run its listed entry point, and compare output checksums and summary counts. Scripts under `candidates/`, `clcnkb/`, and some report-stage modules preserve the original bounded calculations and expected directory layouts; replace identifiers/alleles only when intentionally adapting the method to another case.

`run_bounded_prs.py --helper-dir` must reference Carbon-DeCoder's `pipelines/nextflow/scripts` directory. In particular, it invokes the upstream `parse_plink_pgs_scores.py` contract with metadata, optional PLINK2 `.sscore`, output path, reference build/path, contig style, sample build, and manifest version; this helper is intentionally reused rather than copied.

## Limits

Read-backed phase is informative only for molecules spanning targets; absence of a phase block is unresolved, not evidence of cis or trans. RD/BAF, MoChA, SV, mtDNA, Exomiser, PRS, PGx, taxonomy, and model predictions require independent interpretation and do not establish or exclude a diagnosis. External databases are versioned inputs and may no longer reproduce historical annotations exactly.

## Verification

Run the existing synthetic/offline tests without private data or API calls:

```bash
python -B -m pytest -q -p no:cacheprovider tests
sha256sum -c SHA256SUMS
```

Release verification: **20 offline tests passed**; syntax checks passed for
**70 Python files, nine shell scripts and one C++ scanner**.
Tests verify software contracts, not historical biological findings. The release
was not validated by rerunning the genome workflow or querying a clinical answer
key. Online predictions and full private-data replay remain unexecuted here.

## Licensing and attribution

See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Apache-2.0 notices for the
bundled third-party clients are retained; no unsupported blanket license is
imposed on the remaining sources. The organizer acknowledgement and dataset
source link are in [METHODS.md](METHODS.md#acknowledgement).
