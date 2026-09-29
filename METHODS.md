# Track 1 supplementary computational methods

This is the public code/methods supplement for the
[MVA Hackathon 2026, Track 1](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026).
It is not the individual-level clinical report, a diagnosis, or a competition
submission. Predictions and the full report are supplied separately to the
organizers. No competition answer key or private sequencing/phenotype files are
distributed here.

## Scope and provenance

The source baseline is the screening-reviewed supplementary release of
14 September 2026, internal source revision
`be1b9fae8766a91e33d661e99a9955a5467b5cf5`. Public-release adaptations remove
embedded historical clinical observations, require private phenotype/candidate
configuration, and replace a host-specific AlphaGenome transport with an
explicitly selected native SDK client. These are portability and publication
changes, not new biological analyses or retrospective confirmation of results.
`SOURCE_PROVENANCE.tsv` records source lineage and current bundled file hashes.

Primary sequencing orchestration remains in
[Carbon-DeCoder](https://github.com/paulnovy/Carbon-DeCoder/tree/cd03b3bdb6a29847916429d072f27de94c88065c),
pinned to `cd03b3bdb6a29847916429d072f27de94c88065c` for public reproduction.
This is not a claim that the public commit was the historical execution checkout:
retained historical checkout identification was
`2d2ac14b24b44486a9a6141bd916c6dc72d772f6`, and some execution-version fields
were not populated. The upstream implementation is reused rather than copied.

## Analysis families

| Family | Method and replay contract |
| --- | --- |
| Reference/VCF integrity | Exact REF audit, left normalization, original multiallelic lineage and GT/AD/PL preservation. `scripts/qc/`, `PARTIAL_GT_REPLAY.md`. |
| Clinical annotation | Allele-exact matching against a pinned, normalized ClinVar snapshot; complete and partial-called-ALT strata remain separate. `scripts/clinvar/`. |
| Phenotype prioritization | Exomiser with privately curated HPO profiles, ontology/source auditing, and controlled phenotype sensitivity. `scripts/exomiser/`, `config/phenotype-sensitivity.md`. |
| Phase and candidate follow-up | Regional read/fragment review, public population-response replay, candidate annotation and bounded SV/MEI measurements. `scripts/phase/`, `scripts/candidates/`. |
| Dosage and mosaicism | Masked depth/BAF analysis and pinned MoChA calling with explicit excluded-site selection. `scripts/copy_number/`, `MOCHA_WRAPPER_REPLAY.md`. |
| Structural and mitochondrial variation | Bounded SV/mtDNA calling, filtering and gene annotation; no complete repeat-expansion or MEI sensitivity claim. `scripts/sv_mtdna/`. |
| Targeted reconstruction | CLCN locus read extraction/assembly, exact reference scanning, remapping and transcript context. `CLCN_REPLAY.md`. |
| Molecular prediction | Offline Atlas extraction and saved API response analysis; optional RNA/splice/chromatin predictions and reference-derived sequence contrasts. `scripts/alphagenome/`. |
| Gene-panel audits | Canonical/transcript-aware coding and splice overlap, called-ALT review, exact ClinVar and published-allele/structural audits. `GENE_AUDIT_REPLAY.md`. |
| Ancillary screens | Bounded raw polygenic sums, PharmCAT genotype completeness and unmapped-read taxonomy controls. `scripts/prs_pgx_taxonomy/`. |
| Exploratory controls | Reference-derived protein controls and synthetic bulk-mixture identifiability examples, not patient measurements. `scripts/exploratory/`. |

The detailed input/output/dependency map is [REGISTRY.tsv](REGISTRY.tsv), with
ordered commands in [RECIPES.md](RECIPES.md). Manual public-panel/population
lookups and assay-feasibility assessment are explicitly listed in
[POPULATION_AND_FEASIBILITY_REPLAY.md](POPULATION_AND_FEASIBILITY_REPLAY.md):
no executable producer was recovered for those manual steps.

## Resources and reproducibility limits

Recorded resources include GRCh38 no-alt `GCA_000001405.15`, ClinVar
2026-09-05, Exomiser 15.1.0 with bundle 2512 and Java 21, VEP 115 / REST 15.11,
saved gnomAD v4 responses, Eagle 2.4.1, GATK 4.6.2.0, PharmCAT 3.3.0
(matcher 2.0.0; data 2026-07-05-00-25), PLINK v2.0.0-a.6.9LM,
Kraken2 database 20260226, and AlphaGenome SDK 0.9.0.
The MoChA image separately pins bcftools 1.20 and MoChA `95686b7`.
Module-specific versions must not be conflated into one historical environment.
Historical executable versions not established by retained records are not
invented. The public Python dependency file pins the publication-verification
environment, not the historical research environment.

Reference assemblies, access-controlled inputs, database snapshots, masks,
runtime tools and saved API responses must be supplied by the authorized user.
Several preserved analysis scripts require the documented dated directory
layout; they are replay tools, not a turnkey generic diagnostic pipeline.
Exact re-execution of mutable web services and server-side model versions is not
guaranteed. Source inspection and synthetic tests do not reproduce an individual's
biological results or establish screening completeness.

### Public CLCNKB query target

The 79-base insertion literal in the CLCNKB annotation/prediction scripts is the
exact public gnomAD v4 allele associated with **rs1553127751**, anchored at
GRCh38 `chr1:16050023 T>T+insertion`. Its identity was cross-checked against
the retained public gnomAD response. It is a public alternate allele query
target, **not reference-derived sequence and not a patient read fixture**.
See the [dbSNP public record](https://www.ncbi.nlm.nih.gov/snp/rs1553127751).
Public allele identity does not establish that any individual carries it,
its phase, copy context, or clinical significance.

## Offline and online execution boundary

Default use is local preparation and saved-response replay. Online modes require
`MVA_ALLOW_NETWORK=1`. The bundled native AlphaGenome factory additionally
requires `MVA_ALPHAGENOME_TRANSPORT=native` and an `ALPHAGENOME_API_KEY`
provided through the user's external credential mechanism. No key is included,
printed, passed on the command line, or requested for the offline tests.
Gateway sentinel values are rejected; there is no automatic fallback from a
protected transport to a native connection.

The initial API-stage runner can select the bundled factory with
`--client-factory alphagenome_native_client:create_client`. Other molecular
prediction scripts use the same factory. Install `alphagenome==0.9.0` only
when that module is needed. This public transport differs from the historical
host-specific adapter: the SDK can retry transient RPC failures, and its
`timeout` parameter controls channel readiness, not per-request deadlines.
Task-level checkpoints still prevent automatic restart of an incomplete task.
No live prediction or credential use was performed to prepare this release.

Online requests expose the selected allele/reference-window inputs to the
chosen service; local storage alone is not an exclusively local-processing
guarantee. Review each bounded request and applicable data-access terms first.
Keep generated manifests, source wording, source paths, variant evidence and
reports private. The gated dataset is accessed through the
[organizer dataset entry](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data);
it is not redistributed by this repository.

## Interpretation and reporting

Prediction scores are not calibrated causal probabilities. Phase, exact-allele
function, disease-specific evidence and analytical coverage require independent
interpretation. Absence from a query or a bounded screen is not a benign
classification. Raw PRS sums are not validated risk estimates; unphased PGx
matches are not confirmed diplotypes; taxonomy output is not an infection assay.

This code was assembled and reviewed with AI assistance under human direction.
The historical provider/model/plan/data-handling declaration belongs in the
separate submission report and cannot be inferred from this repository's current
environment. No undocumented historical privacy setting is asserted here.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
