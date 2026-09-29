# MoChA WGS wrapper replay note

The WGS wrapper passes the excluded-variant BCF to MoChA with its exclusion
prefix: `-v ^/excluded/<file>`. This matches the retained historical runner
(`operations/mocha-priority-20260910/run_mocha_priority.sh`) and the retained
annotated BCF header documented in the September 14 filter audit. The leading
`^` is part of MoChA's filter syntax; omitting it changes the site-selection
semantics.

The wrapper change is source-only. `tests/test_run_mocha_wgs.py` uses synthetic
placeholder inputs and a fake `docker` executable to capture argv, asserting
the exact prefixed path and that no container, WGS input, or network/API work
is run. Reproduce it with:

```bash
python3 -m unittest tests/test_run_mocha_wgs.py
```

Historical execution outputs and operation manifests remain read-only. The
original source provenance hash in `SOURCE_PROVENANCE.tsv` is preserved; its
bundled hash is updated for this correction.

The corrected wrapper SHA-256 is
`abde78c078942411774849f8a71915b500658b804e2b069e4de03e6fc45867c4`.

## Audit the retained post-filter without calling MoChA

`scripts/copy_number/audit_mocha_filter.py` replays the saved filter and all
32 combinations of its gates. It requires the original protected directory
layout under `DATA_ROOT`, including the historical calls, stats, AWK, pinned
source copies, review and execution records. It does not run a caller or
change genotypes. This is a case-specific historical replay, not a generic
constitutional-CNV acceptance filter.

Create a **fresh** private output directory under `DATA_ROOT/results/` and
copy the already-verified public-source snapshots into it before running:

```bash
export DATA_ROOT=/path/to/protected/case
export AUDIT_OUT="$DATA_ROOT/results/mocha-filter-replay-new"
mkdir -m 700 "$AUDIT_OUT"
cp "$DATA_ROOT/operations/mocha-interpretation-20260910/upstream-README.md" "$AUDIT_OUT/public-pinned-README.md"
cp "$DATA_ROOT/operations/mocha-interpretation-20260910/upstream-mocha.c" "$AUDIT_OUT/public-pinned-mocha.c"
python3 -B scripts/copy_number/audit_mocha_filter.py --data-root "$DATA_ROOT" --output "$AUDIT_OUT"
```

The source names preserve the audit's comparison contract. Copying them for
replay is not a fresh network verification; pin and verify them against
MoChA commit `95686b7b65f53a490513be76bb120e5fc20a8bcf` separately if needed.
Some protected historical files may need the existing read-only sudo access.
Outputs are exclusively created; do not reuse a populated output directory.
All rejected rows and predicates remain visible. No zero-survivor result is
interpreted as excluding constitutional CNVs or cellular MVA.
