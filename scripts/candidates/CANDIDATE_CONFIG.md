# Private candidate-annotation configuration

`focused_candidate_annotation.py` requires `--manifest /private/manifest.json`
for every mode. Keep the manifest, phase reports, input genomes and generated
results outside the public repository. The script no longer embeds historical
case observations or fixed expected result counts.

The manifest is a JSON object. All fields are optional; `{}` retains the
historical directory conventions under `MVA_DATA_ROOT` and the published target
gene list, but adds no extra variants or inheritance interpretations.

| Field | Value |
| --- | --- |
| `input_paths` | Object with optional `reference_fasta`, `normalized_vcf`, `clinvar_vcf`, `exomiser_root` paths. Relative paths resolve beside the manifest. |
| `genes` | Array of target gene symbols used to select Exomiser rows. |
| `audit_genes` | Array of genes whose Exomiser row counts should be reported; presence is measured, not assumed absent. |
| `extra_variants` | Array of objects with `gene`, `chrom`, one-based `pos`, `ref`, `alt`. Each is checked against the supplied FASTA and exact normalized VCF allele. |
| `clinvar_snapshot_date` | Caller-supplied snapshot date string; absent means unknown, not a fabricated date. |
| `phase_evidence` | Object mapping arbitrary labels to existing private JSON reports. The reports are copied with computed hashes and labeled not independently verified. |
| `inheritance_annotations` | Array of caller-supplied annotation objects, each requiring `gene`. Optional `phase` is `cis_supported`, `trans_supported`, or `unresolved`; only boolean `pathogenic_alleles_confirmed: true` enables the trans/biallelic rule. All supplied text is labeled as caller interpretation, not calculated evidence. |
| `supplemental_sources` | Object mapping provenance labels to private file paths; hashes are computed from those files. |

Modes remain `select`, `query`, `query_archive`, `gnomad`, `gnomad_qc`,
`assemble`, and `verify`. Run `select`, `query_archive`, `gnomad_qc`, `assemble`,
then `verify` to produce the assembled archive-based report. Network modes still
require the explicit `MVA_ALLOW_NETWORK=1` opt-in and send selected variant
coordinates to the relevant public annotation API. Outputs use the existing
`results`, `reports`, and `operations` subdirectories beneath `MVA_DATA_ROOT`;
use a fresh private root for a new analysis to avoid reusing historical outputs.

`analyze_targeted_sv_mei.py` keeps its existing CLI and calculation methods.
Target intervals now have generic labels (`clcn_window_a`, `clcn_window_b`,
`clcn_target`) rather than historical caller IDs. Reports contain calculated
counts and neutral interpretation-needed statuses, not static positive/negative
findings. `target_to_mean_control_depth_ratio_mapq20` replaces the historically
misnamed depth-ratio key; it is null when the control mean is zero. The
`--delly-records` input is fingerprinted only, not parsed or clinically assessed.
Target-specific coordinates are not a genome-wide assay.
