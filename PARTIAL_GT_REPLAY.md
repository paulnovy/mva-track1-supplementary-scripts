# Offline partial-GT screening-integrity audit

This audit traces a split call back to its original allele identities. It does
not rewrite missing GT as reference, infer phase, call variants, align reads,
rerun Exomiser/PRS, or make annotation/prediction API requests.

Run from this checkout using an existing environment with `pysam`. Both runners
require fresh private output directories; all patient-level outputs belong
outside Git. `PYTHONDONTWRITEBYTECODE=1` avoids modifying imported historical
environments. A run uses one sequential stream, normally under 1 GiB RAM.

```bash
python scripts/qc/audit_partial_gt.py \
  --original /private/supported-original.bcf \
  --normalized /private/normalized.vcf.gz \
  --source /private/clinvar-source/annotated.vcf.gz \
  --core /private/clinvar-core/annotated.vcf.gz \
  --fasta /reference/exact.fa --aliases /reference/aliases.tsv \
  --transcripts /reference/2512_hg38_transcripts_ensembl.ser \
  --exomiser /private/baseline.variants.tsv /private/full.variants.tsv /private/axis.variants.tsv \
  --output /private/new-audit/verified

python scripts/qc/audit_partial_gt_followup.py \
  --audit /private/new-audit/verified \
  --source /private/clinvar-source/annotated.vcf.gz \
  --core-original /private/core.raw.vcf.gz \
  --clinvar /reference/clinvar.normalized.vcf.gz \
  --reports /private/retained-reports --exomiser-log /private/exomiser.stdout.log \
  --panel nine_gene=/private/panel_alt_present_canonical_classified.tsv \
          expanded_transcripts=/private/all-sources-transcript-audit.tsv \
          noncoding=/private/noncoding_eligibility.tsv \
  --output /private/new-audit/followup
```

The first runner reuses `clinvar_exact.terms`, `clinvar_triage.genotype_state`,
`audit_panel_clinvar_snapshot.called_alts`, and the canonical audit's
`changed_span`. A small representation-only bcftools replay adds an original
allele index, then verifies GT/AD/PL/DP/GQ/phasing against both retained files.
The original multiallelic BCF remains intact alongside the JSONL lineage.
Exact allele comparisons deliberately preserve REF/ALT identity and reject
ambiguous duplicate partial alleles. Core recovery is exact representation
recovery, not proof of independent sample or biological validation.

The transcript reader supports the retained Exomiser 15.1.0 `JTPB` gzip/protobuf
snapshot, using the published `jannovar.proto` field numbers and strand-oriented
zero-based intervals. It streams transcript messages rather than loading their
sequences en masse. Coding overlap uses only transcripts present in that snapshot;
it is not a new consequence annotation or comprehensive transcript survey.

The second runner resolves core partial candidates to original raw records using
`NORMALIZATION_ORIGINAL`, reads exact saved ClinVar evidence, and joins normalized
or provenance-linked original alleles to the supplied follow-up TSVs. Its report
search detects coordinate mentions only in retained MD/JSON/TSV <=5 MB. Neither
table presence nor a coordinate mention establishes clinical review. No mention
does not prove absence from every possible document.

## Interpretation boundary

The screening-reviewed release changes `clinvar_exact.py` prospectively: original
complete-GT counters remain, with additional called-ALT and partial-called-ALT
classification counters. `clinvar_triage.py` includes explicitly carried partial
P/LP matches, retaining GT and a required-original-genotype-review flag. It never
fills a dot or treats split PL/AD as reconstructed diploid evidence. Historical
VCFs and summary files are not rewritten; use the lineage artifacts above for
original-context recovery, and keep uncertain/conflicting candidates in the
follow-up review ledger. No new Exomiser scores are implied by that ledger.

Completeness and called-ALT presence are distinct predicates. A complete original
`1/2` can become `1/.` and `./1`; the other ALT is not reference. Original AD has
one slot per original allele, and original diploid PL has triangular genotype
ordering. Split AD/PL must not be used to invent a complete `0/1` call.

Exomiser `PASS_ONLY` retains variants passing the configured analysis filters;
it does not mean literal VCF `FILTER=PASS`. `failedVariantFilter` admits both PASS
and dot. Input allele presence, parsed unknown alleles, configured quality and
inheritance filters, and final TSV presence require separate checks. Retained
outputs are evidence of recovery; absent output is not by itself evidence that
partial GT caused exclusion. A counterfactual full reanalysis is not performed.

Primary format/implementation references (15.1.0):

- https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/proto/jannovar.proto
- https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/java/org/monarchinitiative/exomiser/core/genome/jannovar/JannovarProtoConverter.java
- https://github.com/exomiser/Exomiser/blob/15.1.0/exomiser-core/src/main/java/org/monarchinitiative/exomiser/core/genome/VariantContextSampleGenotypeConverter.java
- https://exomiser.readthedocs.io/en/stable/advanced_analysis.html

Focused verification:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p test_partial_gt_audit.py -v
git diff --check
```
