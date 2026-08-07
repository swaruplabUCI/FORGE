# T2 tutorial — measured resource baseline

**Authoritative run: `results_tutorial_remeasure`, SLURM job 55119808,
2026-08-07, on the FIXED tier.** All numbers below are measured from
`trace.tsv`, not estimated. Raw artifacts are preserved beside this file in
`remeasure_trace/` because `results*/` is gitignored.

Run 1 (the broken-tier run that produced the fixes) is kept in `run1_trace/`
for comparison only. Do not cite its resource numbers.

## How it was run

```bash
sbatch -A <account> -p <partition> launch_tutorial.sh \
    --no-resume --outdir results_tutorial_remeasure
```

Cold run (no `-resume`), 8 CPUs / 48 GB allocation, local executor inside one
SLURM job, commit `c2415ab` (tier fix) + uncommitted launcher.

## Headline

| Metric | Value |
|---|---|
| **Wall-clock** | **1 h 42 m 45 s** (8 CPUs) |
| **CPU-hours** | **6.6** |
| Tasks | **94 / 94 COMPLETED**, 0 failed, 0 cached |
| Summed task realtime | 2 h 49 m |
| **Peak single-task RSS** | **8.80 GB** (`ATAC_FINAL_PIPELINE`) |
| `results_tutorial_remeasure/` | **3.88 GB** |
| `work/` for this run alone | **4.81 GB** |
| Recommended free disk | **~15 GB** |

More cores help less than expected: run 1 on 10 CPUs took 1 h 33 m, only ten
minutes faster. Wall-clock is dominated by a few long serial tasks
(`HDWGCNA_ENRICHMENT` 53 m, `CELLBENDER` 23 m, `RUN_CELLCHAT` 22 m,
`MULTIVI_INTEGRATE` 19 m). Extra cores buy the 45-way hdWGCNA and 25-way Cicero
fan-outs, nothing else.

## Per-process

| process | n | total | longest | peak RSS | allocated | headroom |
|---|---:|---:|---:|---:|---|---:|
| `RNA:HDWGCNA_ENRICHMENT` | 2 | 58m11s | 53m30s | 7.80 GB | 2c / 12 GB | 35% |
| `RNA:CELLBENDER` | 1 | 23m15s | 23m15s | 2.90 GB | 2c / 4 GB | 28% |
| `RNA:RUN_CELLCHAT` | 1 | 22m30s | 22m30s | 1.10 GB | 2c / 6 GB | 82% |
| `MULTIOME_INTEGRATION:MULTIVI_INTEGRATE` | 1 | 18m49s | 18m49s | 7.80 GB | 4c / 12 GB | 35% |
| `RNA:HDWGCNA_PER_CELLTYPE` | 45 | 24m47s | 5m45s | 6.60 GB | 2c / 10 GB | 34% |
| `REGULATORY_ANALYSIS:CICERO_ESTIMATE_DP` | 1 | 2m36s | 2m36s | 1.60 GB | 2c / 4 GB | 60% |
| `ATAC_INITIAL:ATAC_INITIAL_QC` | 1 | 2m05s | 2m05s | 7.10 GB | 4c / 8 GB | 11% |
| `ATAC_FINAL:ATAC_FINAL_PIPELINE` | 1 | 1m47s | 1m47s | 8.80 GB | 4c / 12 GB | 27% |
| `RNA:RUN_CELLTYPIST` | 1 | 1m34s | 1m34s | 0.57 GB | 4c / 8 GB | 93% |
| `MULTIOME_INTEGRATION:MOFA_INTEGRATE` | 1 | 1m17s | 1m17s | 1.00 GB | 4c / 8 GB | 88% |
| `RNA:CONCAT_BATCHES` | 1 | 47s | 47s | 1.10 GB | 4c / 12 GB | 91% |
| `MULTIOME_INTEGRATION:MOFA_VISUALIZE` | 1 | 44s | 44s | 0.92 GB | 4c / 12 GB | 92% |
| `RNA:PLOT_POST_SCANVI` | 1 | 36s | 36s | 1.20 GB | 4c / 8 GB | 85% |
| `MULTIOME_INTEGRATION:MULTIVI_VISUALIZE` | 1 | 36s | 36s | 4.90 GB | 4c / 8 GB | 39% |
| `REGULATORY_ANALYSIS:CICERO_JOIN` | 1 | 28s | 28s | 1.30 GB | 2c / 4 GB | 68% |
| `RNA:RNA_QC` | 1 | 28s | 28s | 0.58 GB | 4c / 8 GB | 93% |
| `REGULATORY_ANALYSIS:CICERO_FULL_CHROM` | 25 | 7m43s | 27s | 1.20 GB | 2c / 4 GB | 70% |
| `RNA:CONVERT_H5AD_TO_SEURAT` | 1 | 26s | 26s | 1.30 GB | 2c / 6 GB | 78% |
| `MULTIOME_INTEGRATION:EXPORT_MUDATA_RNA` | 1 | 18s | 18s | 0.59 GB | 2c / 6 GB | 90% |
| `REGULATORY_ANALYSIS:CICERO_TRIPLETS` | 1 | 15s | 15s | 0.64 GB | 2c / 4 GB | 84% |
| `ATAC_FINAL:MERGE_ANNOTATIONS` | 1 | 11s | 11s | 0.74 GB | 2c / 4 GB | 82% |
| `ATAC_DESCRIPTIVE_REPORT` | 1 | 11s | 11s | 0.76 GB | 2c / 4 GB | 81% |
| `MULTIOME_INTEGRATION:BUILD_MUDATA` | 1 | 10s | 10s | 0.87 GB | 4c / 12 GB | 93% |
| `ATAC_FINAL:ATAC_CELLTYPE_ANNOTATION` | 1 | 2s | 2s | 0.06 GB | 2c / 4 GB | 98% |
| `ATAC_INITIAL:ATAC_MAKE_THRESHOLDS` | 1 | 2s | 2s | 0.06 GB | 2c / 4 GB | 98% |

**No process exceeded its allocation.** Two are tighter than intended and should
be raised before publication, because peak RSS is NOT stable between runs:

| Process | run 1 peak | re-measure peak | allocated | note |
|---|---|---|---|---|
| `ATAC_INITIAL_QC` | 5.80 GB | **7.10 GB** | 8 GB | only 11% headroom — raise to 10 GB |
| `CELLBENDER` | 1.50 GB | **2.90 GB** | 4 GB | nearly doubled run to run — raise to 6 GB |

That run-to-run variance is the single most important caveat in this document.
Sizing memory to a single observation is not safe; use ~50% headroom minimum.

## Structural results — REPRODUCIBLE

Every structural count is **identical** between run 1 and the re-measure,
including a float to 15 significant figures. This is direct evidence that the
scvi-tools seeding (`params.random_seed = 42`) works end to end.

| Quantity | Value | Reproducible? |
|---|---|---|
| ATAC cells past initial QC | 944 of 1,000 | identical |
| ATAC median fragments | 1,419.5 | identical |
| ATAC median TSS enrichment | 16.300101023624922 | identical to 15 s.f. |
| ATAC final (peak matrix) | 817 cells x 12,085 peaks | identical |
| Cicero input | 817 cells / 12,085 peaks | identical |
| MuData (RNA ∩ ATAC) | 767 cells, 21,014 genes | identical |
| MOFA+ | 3 factors, seed 42 | identical |
| CellTypist labels | 45 | identical |
| hdWGCNA cell types with modules | 45 | identical |

## GPU stripping — verified by ground truth

```
tasks with '--nv' in .command.run  : 0
tasks with 'gres=gpu' in .command.run : 0
```

This is the only trustworthy check. See the verification section below.

## The tier fix landed

Allocation distribution, broken tier vs fixed:

| | broken (run 1) | fixed (re-measure) |
|---|---|---|
| 2c / 8 GB | **88 tasks** | 0 |
| 4c / 16 GB | 6 tasks | 0 |
| 2c / 10 GB | — | 45 |
| 2c / 4 GB | — | 33 |
| 4c / 8 GB | — | 6 |
| 4c / 12 GB | — | 5 |
| 2c / 6 GB | — | 3 |
| 2c / 12 GB | — | 2 |

Six distinct tiers instead of one flat wildcard value. See commit `c2415ab` for
the Nextflow selector-precedence rule that caused it.

## KNOWN REPORTING BUG — do not trust `atac_pipeline_summary.json`

`results_*/atac/final/atac_pipeline_summary.json` reports

```json
"filtering_thresholds": {"min_counts": 5000, "min_tsse": 6, "max_counts": 100000}
```

Those are `bin/atac_consolidated_pipeline.py`'s **argparse defaults**, not what
was applied. `modules/atac/consolidated_pipeline.nf` never passes
`--min_counts` / `--min_tsse` / `--max_counts`; the script filters using
`thr.get(key, args.key)` where `thr` comes from `--thresholds_file`. The
thresholds actually applied in this run were the computed per-sample values:

```json
{"TUTORIAL_PBMC_tutorial": {"min_counts": 660.15,
                            "max_counts": 14886.7,
                            "min_tsse": 6.30572243680171}}
```

which is what takes 944 cells down to 817. Filtering is correct; the summary
field is misleading. Read `atac/initial_qc/sample_thresholds.json` instead.

## How to verify the tier — and how NOT to

**Counting strings in `nextflow config` output is an INVALID check.** The Step 2
verification used `gres=gpu 0 / accelerator=1 0 / --nv 0` on the resolved config
and read that as "GPU fully stripped". Those counts were 0 only because the
tutorial tier's GPU-strip selector was *textually identical* to
nextflow.config's, so the two blocks merged and the later value overwrote the
earlier one in the dump. Rewriting the selector to the qualified-safe form made
it a distinct block, nextflow.config's `accelerator = 1` / `'--nv'` reappeared in
the dump, and the counts went 0 -> 1 — with no change whatsoever in what a task
actually gets. The dump lists blocks; it does not resolve them per process.

Valid checks, in increasing order of strength:

1. Structural, on the resolved config: no `cpus`/`memory`/`time` inside the
   `withName: '.*'` block; all 8 `withLabel` blocks present.
2. Empirical precedence, via a throwaway two-process `.nf` plus a `trace` scope
   emitting `name,cpus,memory`. This is how the
   `'ATAC_.*'`-shadows-`'ATAC_FINAL_PIPELINE'` hazard was caught.
3. Ground truth, after a real run: read `trace.tsv` for cpus/memory per process,
   and grep `.command.run` for `--nv`. Only this exercises the real merge for
   every process in the DAG. **This is what was done above.**

## Regenerating

1. Runtime/memory — `<outdir>/pipeline_info/trace.tsv`. NOTE: the trace only
   lands in `<outdir>` when `--outdir` is passed on the CLI. With the profile
   default it goes to `results/pipeline_info/` instead, because
   `${params.outdir}` is interpolated when `nextflow.config` is parsed, before
   the `tutorial` profile sets it. A stale `results/pipeline_info/trace.tsv`
   from an earlier run is a real trap — check the mtime before trusting it.
2. Structural counts — the per-stage JSON summaries under the output tree.
3. Disk — sum `stat -c %s`; `du -sh` under-reports fresh writes on BeeGFS.
4. Allocations — always confirm from `trace.tsv`. The failure mode is silent.
