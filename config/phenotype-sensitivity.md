# Private phenotype sensitivity configuration

`prepare_mva1_phenotype_sensitivity.py` requires `--config` pointing to a private
JSON file; no patient phenotype or family context is embedded in the script.
Keep that file, the source document, VCF, and generated outputs outside this
public repository. Outputs include source wording, sample IDs, and local paths.

```bash
python scripts/exomiser/prepare_mva1_phenotype_sensitivity.py \
  --docx /private/phenotype.docx \
  --config /private/phenotype-sensitivity.json \
  --hpo /resources/hp.obo --vcf /private/input.vcf.gz \
  --output /private/work/phenotypes
```

The following is a **synthetic schema example**, not clinical input. Replace the
invented HPO identifier and label with the exact pinned ontology entry and curate
the private source document yourself. Paragraph numbers are 1-based Word XML
paragraph positions, including empty paragraphs.

```json
{
  "curation": [
    {
      "key": "synthetic_feature",
      "hpo_id": "HP:9999999",
      "label": "Synthetic example only",
      "feature_paragraphs": [1],
      "onset": "not stated",
      "subject": "proband",
      "justification": "Replace with a source-supported private rationale."
    }
  ],
  "profiles": {
    "baseline": ["HP:9999999"],
    "full": ["HP:9999999"],
    "axis": ["HP:9999999"]
  }
}
```

`curation` and `profiles` are required and nonempty. Each curation row must contain
the illustrated keys and a unique HPO ID. Only present features with `subject`
exactly `proband` may enter scored profiles. Rows for other subjects can remain
in the curation audit but must not be included in `profiles`.

Profile names use letters, digits, hyphens, or underscores. HPO IDs within each
profile are unique and their supplied order is preserved. The existing
`run_exomiser_sensitivity.sh` runner expects profiles named `baseline`, `full`,
and `axis`; define their contents privately according to the intended analysis.

Optional `family_history` is private JSON context written unchanged to
`family-history.json`, separate from scored Phenopackets; omit it when there is
no curated family context. Use a fresh output directory for each configuration.
Each run records the configuration SHA-256 in `provenance.json` alongside the
document, VCF, ontology, and generated Phenopacket hashes.
