# Reproduction recipes

All paths below are placeholders. Use versioned, checksum-verified inputs and a new output directory. Analysis scripts require `MVA_ALLOW_NETWORK=1` for public API requests. Installation and image-building commands download public software independently; provision these separately before offline analysis.

## Baseline VCF and ClinVar

```bash
python scripts/qc/audit_vcf_reference.py --vcf case.vcf.gz --fasta GRCh38.fa --aliases contig-aliases.tsv --output work/ref-audit
python scripts/qc/normalize_vcf.py --audit work/ref-audit/summary.json --aliases contig-aliases.tsv --output work/normalized
python scripts/qc/vcf_qc.py --normalization work/normalized/summary.json --output work/vcf-qc
python scripts/clinvar/prepare_clinvar.py --help  # select snapshot/build arguments supported by the historical implementation
python scripts/clinvar/clinvar_exact.py --help
```

## Phase, RD/BAF, and MoChA

```bash
python scripts/phase/prepare_regional_phasing.py prepare --bam case.bam --vcf case.vcf.gz --fasta GRCh38.fa --region chr15:START-END --target chr15:POSITION --out work/phase
python scripts/phase/audit_regional_fragment_phase.py --help
python scripts/copy_number/chromosome_rd_baf.py --regions mosdepth.regions.bed.gz --fasta GRCh38.fa --vcf case.vcf.gz --out work/chromosome-rd-baf
python scripts/copy_number/rd_baf_followup.py --prior-windows work/windows.tsv --prior-summary work/summary.json --regions mosdepth.regions.bed.gz --umap-bed umap.bed --superdups genomicSuperDups.bed --gaps gaps.bed --centromeres centromeres.bed --out work/rd-baf-followup
bash scripts/copy_number/run_mocha_wgs.sh phased.bcf excluded-variants.bcf cnps.bed work/mocha SAMPLE_ID
MVA_DATA_ROOT="$PWD/workspace-layout" python scripts/copy_number/interpret_mocha.py
```

MoChA preparation expects one dense, phased BCF containing common SNPs, a CSI index, MoChA's excluded-variant BCF/CSI, and a CNP BED for the same genome build. The caller script pins the historical container tag `mva-mocha:bcftools-1.20-95686b7` and passes the excluded sites as `^/excluded/<file>`; build or substitute an equivalent image containing bcftools 1.20 and MoChA commit `95686b7`. See [MOCHA_WRAPPER_REPLAY.md](MOCHA_WRAPPER_REPLAY.md) for the correction provenance and offline argv test.

## SV, mtDNA, phenotype sensitivity, PRS, PGx, taxonomy

```bash
bash scripts/sv_mtdna/run_mtdna_sv_screen.sh --mode sv --bam case.bam --reference GRCh38.fa --scope scope.json --out work/sv
bash scripts/sv_mtdna/filter_existing_mtdna.sh GRCh38.fa mtdna.native.vcf.gz work/mtdna-filtered
python scripts/exomiser/prepare_mva1_phenotype_sensitivity.py --docx phenotype.docx --config /private/phenotype-sensitivity.json --hpo hp.obo --vcf case.vcf.gz --output work/phenotypes
bash scripts/exomiser/run_exomiser_sensitivity.sh work/phenotypes work/exomiser-sensitivity work/exomiser-ops /path/to/exomiser-resources case.vcf.gz /path/to/Carbon-DeCoder SAMPLE_ID
bash scripts/exomiser/run_exomiser_followup.sh work/exomiser work/ops case.vcf.gz /path/to/Carbon-DeCoder
python scripts/prs_pgx_taxonomy/run_bounded_prs.py --help
bash scripts/prs_pgx_taxonomy/run_pgx_screen.sh --bam case.bam --reference GRCh38.fa --positions pharmcat.vcf.bgz --preprocessor pharmcat-preprocessor --jar pharmcat.jar --out work/pgx --sample SAMPLE_ID
bash scripts/prs_pgx_taxonomy/run_taxonomy_control.sh --bam case.bam --db kraken-db --output work/taxonomy --sample SAMPLE_ID
```

Historical target versions were Exomiser 15.1.0 with its 2512 phenotype bundle, Java 21, bcftools/samtools 1.20-era interfaces, PharmCAT as recorded by the local JAR manifest, and Kraken2 with a separately versioned database. Exact external resource versions are part of the input contract, not downloaded automatically.

## Offline API replay and targeted mechanisms

The cross-family CLCNKB scripts import the AlphaGenome analysis helpers; add that
source directory to the process search path:

```bash
export PYTHONPATH="$PWD/scripts/alphagenome${PYTHONPATH:+:$PYTHONPATH}"
```

```bash
python scripts/phase/query_gnomad_schema.py --variants 15-EXAMPLE-A-G 15-EXAMPLE-C-T --saved-response saved-gnomad.json --output work/cooccurrence.json
python scripts/alphagenome/query_alphagenome_zipstored.py --help
MVA_DATA_ROOT="$PWD/workspace-layout" python scripts/alphagenome/alphagenome_mechanism.py --mode prepare
MVA_DATA_ROOT="$PWD/workspace-layout" python scripts/clcnkb/clcnkb_alphagenome.py --mode analyze
```

Online AlphaGenome 0.9.0 prediction, VEP, gnomAD, and dbSNP stages additionally require `MVA_ALLOW_NETWORK=1`; standalone AlphaGenome additionally requires `MVA_ALPHAGENOME_TRANSPORT=native` and an externally provisioned API key. The initial API-stage runner accepts `--client-factory alphagenome_native_client:create_client`. See `METHODS.md` for retry and timeout semantics. Never put credentials in these commands or configuration files.

## BUB1B hidden-allele and CLCN reconstruction

```bash
python scripts/candidates/analyze_bub1b_hidden_alleles.py \
  --sam bounded-region.sam --depth-q0 depth.q0.tsv --depth-q20 depth.q20.tsv \
  --depth-q60 depth.q60.tsv --reference region.fa --mane mane.json \
  --alu alu-consensus.json --delly regional-delly.tsv --chrom chr15 \
  --start START --end END --output work/bub1b-hidden

# The ordered CLCN recipe includes required reference scans and remapping.
# Follow CLCN_REPLAY.md; do not run the stage names as an uninterrupted loop.
```

The BUB1B SAM must contain the bounded region and mates; depth tables correspond to MAPQ 0/20/60. [CLCN_REPLAY.md](CLCN_REPLAY.md) gives the container paths and complete interstage producer order. Its fixed unitig identifiers intentionally reproduce this case; they are not a generic allele-discovery interface.

The abbreviated CLCN loop above is not the full interstage chain. Use `CLCN_REPLAY.md`, which includes the exact-reference scanner, mandatory `count` stage, and minimap2 producers.

Produce the BUB1B transcript metadata and FASTA consumed by the AlphaGenome profile scripts:

```bash
python scripts/alphagenome/reconstruct_bub1b_transcripts.py --reference GRCh38.fa \
  --annotation ncbiRefSeqCurated.txt.gz --outdir work/transcripts
cp work/transcripts/transcript_reconstruction.json \
  "$MVA_DATA_ROOT/reports/alphagenome-mechanism-20260910/transcript_reconstruction.json"
```

For the fixed candidate-set annotation workflow:

```bash
export MVA_DATA_ROOT="$PWD/private-case-layout"
python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json select
MVA_ALLOW_NETWORK=1 python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json query
MVA_ALLOW_NETWORK=1 python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json query_archive
MVA_ALLOW_NETWORK=1 python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json gnomad_qc
python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json assemble
python scripts/candidates/focused_candidate_annotation.py --manifest /private/candidate-manifest.json verify
```

Bundled Ensembl, gnomAD, and dbSNP clients retain their Apache-2.0 notices and use the `polite-http` package. They are overridable through `config/example.env` but no private helper tree is required.

## Preserved directory-layout analyses

The historical MoChA interpretation, SV gene review, focused annotation, and AlphaGenome mechanism scripts intentionally preserve their calculations. They resolve inputs below `MVA_DATA_ROOT` using `results/`, `reports/`, `operations/`, `resources/`, and `references/` subdirectories named in the script constants. Create symlinks or copies inside a private case workspace to satisfy this schema; do not edit the scripts to point at protected data when producing a public archive. Layout-bound report stages refuse an existing output directory. The registry marks them as layout-bound rather than universally portable.
# September 12–14 additions

See [GENE_AUDIT_REPLAY.md](GENE_AUDIT_REPLAY.md) for the nine-gene canonical
pilot, BUBR1 sequence-control construction, synthetic mixture demonstration,
and TRIM37/CEP192 complete offline annotation/callability/SV replay. Snapshot
collection is optional and network-gated; the replay itself uses saved inputs.

## Private configuration for public-release replay

Phenotype sensitivity requires the private schema documented in
[config/phenotype-sensitivity.md](config/phenotype-sensitivity.md).
Candidate annotation requires [CANDIDATE_CONFIG.md](scripts/candidates/CANDIDATE_CONFIG.md);
its private phase/inheritance annotations are not independently verified by the script.
These explicit inputs replace historical observations formerly embedded in local source.

SV gene annotation also requires a private JSON object mapping HPO identifiers
to labels (synthetic format example: `{"HP:9999999": "Synthetic term"}`):

```bash
python scripts/sv_mtdna/annotate_sv_genes.py --hpo-terms /private/hpo-terms.json
```

Its decoded SV, RefSeq, Exomiser and HPO database inputs follow the existing
`MVA_DATA_ROOT` directory contract. No patient HPO combination or fixed sample
callset count is supplied by the public code.
