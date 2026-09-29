# CLCNKA/CLCNKB reconstruction replay

This documents the complete interstage contract omitted from the shorter overview. `clcn_reconstruct.py` preserves the historical `/ref/GRCh38_standard.fa` and `/out` container paths. Work in a new private output directory containing `raw/locus-with-mates.bam` and its index.

Create `work/` in this repository and compile the scanner before entering the runtime container. The runtime needs Python 3.11+, `pysam`, `minimap2` and a compatible C++ runtime. No special image is implied.

Build the exact reference scanner:

```bash
mkdir -p work
g++ -O3 -std=c++17 scripts/clcnkb/clcn_exact_reference_scan.cpp -o work/clcn_exact_reference_scan
```

For container execution, prepare a fresh `work/clcn/raw` and `work/clcn/analysis`; place the independently extracted indexed locus-with-mates BAM in `raw/`. Mount the entire reference directory read-only (FASTA plus `.fai`), and this repository read-only. Replace the illustrative image with your provisioned environment:

```bash
docker run --rm -it --network none \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/app:ro" -w /app \
  -v "$PWD/reference:/ref:ro" -v "$PWD/work/clcn:/out:rw" \
  your-provisioned-clcn-image bash
```

Inside that environment, run the producer chain in order (the Python environment needs `pysam`; minimap2 used `-x sr -c --secondary=yes -N 50 -p 0.1`):

```bash
python scripts/clcnkb/clcn_reconstruct.py prepare
work/clcn_exact_reference_scan /out/analysis/scan-probes.tsv /ref/GRCh38_standard.fa /out/analysis/whole-reference-kmers.tsv
python scripts/clcnkb/clcn_reconstruct.py assemble
python scripts/clcnkb/clcn_reconstruct.py count
minimap2 -x sr -c --secondary=yes -N 50 -p 0.1 /out/analysis/candidate-haplotypes.fa /out/raw/prior-junction-templates.fa > /out/analysis/junction-read-models.paf
minimap2 -x sr -c --secondary=yes -N 50 -p 0.1 /out/analysis/candidate-haplotypes.fa /out/analysis/unitigs-k51.fa > /out/analysis/unitigs-models.paf
python scripts/clcnkb/clcn_reconstruct.py refine
work/clcn_exact_reference_scan /out/analysis/refined-scan-probes.tsv /ref/GRCh38_standard.fa /out/analysis/refined-whole-reference-kmers.tsv
minimap2 -x sr -c --secondary=yes -N 50 -p 0.1 /out/analysis/refined-haplotypes.fa /out/raw/all-primary-reads.fa > /out/analysis/all-read-refined-models.paf
python scripts/clcnkb/clcn_reconstruct.py analyze
python scripts/clcnkb/clcn_reconstruct.py depth
```

Run every command in the same container namespace; do not create global host `/out` or `/ref` paths merely to reproduce this analysis. Retain the original bounded region/mate selection and full reference; arbitrary re-extraction or different ordering may change the fixed unitig identifiers in `refine` (a case replay constraint). The `count` stage is required after the first whole-reference scan. Both scanner outputs are required by `analyze`. The reconstruction remains a bounded sequence hypothesis and does not infer genotype, copy number, or phase.
