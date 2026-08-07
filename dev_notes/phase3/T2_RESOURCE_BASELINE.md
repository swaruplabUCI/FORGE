# T2 tutorial — measured resource baseline (run 1)

**Measured 2026-08-07.** First successful end-to-end execution of the tutorial
dataset. Numbers here are measured from `pipeline_info/trace.tsv`, not estimated.
Re-measure after any change to `configs/resource_tiers/tutorial.config`.

## Run identity

| | |
|---|---|
| Repo commit | `01a7485` (+ uncommitted `configs/datasets/tutorial_pbmc.config`) |
| Command | `nextflow run main.nf -profile tutorial,singularity -c configs/datasets/tutorial_pbmc.config -resume` |
| Host | HPC3 `hpc3-14-17`, 11 cores / 188 GB node, inside a 66 GB SLURM cgroup |
| Executor | `local`, capped to 10 CPUs / 58 GB for this run (session-only `-c` file, not in the repo) |
| Dataset | 20,000 barcodes → 1,000 RNA cells / 944 ATAC cells, chr21+chr22 |

## Headline

| Metric | Value |
|---|---|
| **Wall-clock** | **1h 33m 05s** |
| **CPU-hours** | **6.4** |
| Tasks succeeded | **94 / 94** (0 failed, 0 cached — effectively a cold run) |
| Summed task realtime | 3h 07m |
| Peak single-task RSS | **8.20 GB** (`ATAC_FINAL_PIPELINE`) |
| `results_tutorial/` on disk | **3.85 GB** |

Wall-clock is well under summed realtime because the local executor ran up to 10
tasks concurrently. On a 4-core laptop expect substantially longer — the
45-way `HDWGCNA_PER_CELLTYPE` fan-out and the 25-way `CICERO_FULL_CHROM`
fan-out are what parallelism buys you here.

## Per-process

| process | n | total | longest | peak RSS | alloc |
|---|---:|---:|---:|---:|---:|
| `RNA:HDWGCNA_ENRICHMENT` | 2 | 1h02m | 56m54s | 7.70 GB | 2c / 8 GB |
| `RNA:RUN_CELLCHAT` | 1 | 27m17s | 27m17s | 0.99 GB | 2c / 8 GB |
| `MULTIOME_INTEGRATION:MULTIVI_INTEGRATE` | 1 | 21m38s | 21m38s | 7.80 GB | 2c / 8 GB |
| `RNA:CELLBENDER` | 1 | 19m28s | 19m28s | 1.50 GB | 2c / 8 GB |
| `RNA:HDWGCNA_PER_CELLTYPE` | 45 | 32m07s | 6m45s | 6.20 GB | 2c / 8 GB |
| `MULTIOME_INTEGRATION:MOFA_INTEGRATE` | 1 | 3m27s | 3m27s | 0.98 GB | 2c / 8 GB |
| `REGULATORY_ANALYSIS:CICERO_ESTIMATE_DP` | 1 | 2m51s | 2m51s | 1.50 GB | 2c / 8 GB |
| `ATAC_FINAL:ATAC_FINAL_PIPELINE` | 1 | 1m34s | 1m34s | 8.20 GB | 4c / 16 GB |
| `RNA:RUN_CELLTYPIST` | 1 | 1m33s | 1m33s | 0.57 GB | 2c / 8 GB |
| `ATAC_INITIAL:ATAC_INITIAL_QC` | 1 | 1m29s | 1m29s | 5.80 GB | 4c / 16 GB |
| `MULTIOME_INTEGRATION:MOFA_VISUALIZE` | 1 | 1m02s | 1m02s | 0.91 GB | 2c / 8 GB |
| `RNA:CONCAT_BATCHES` | 1 | 48s | 48s | 1.10 GB | 2c / 8 GB |
| `MULTIOME_INTEGRATION:MULTIVI_VISUALIZE` | 1 | 40s | 40s | 4.90 GB | 2c / 8 GB |
| `RNA:PLOT_POST_SCANVI` | 1 | 37s | 37s | 1.10 GB | 2c / 8 GB |
| `RNA:RNA_QC` | 1 | 31s | 31s | 0.78 GB | 2c / 8 GB |
| `RNA:CONVERT_H5AD_TO_SEURAT` | 1 | 29s | 29s | 1.30 GB | 2c / 8 GB |
| `REGULATORY_ANALYSIS:CICERO_JOIN` | 1 | 28s | 28s | 1.30 GB | 2c / 8 GB |
| `REGULATORY_ANALYSIS:CICERO_FULL_CHROM` | 25 | 7m37s | 27s | 1.20 GB | 2c / 8 GB |
| `MULTIOME_INTEGRATION:EXPORT_MUDATA_RNA` | 1 | 24s | 24s | 0.67 GB | 2c / 8 GB |
| `REGULATORY_ANALYSIS:CICERO_TRIPLETS` | 1 | 14s | 14s | 0.62 GB | 2c / 8 GB |
| `MULTIOME_INTEGRATION:BUILD_MUDATA` | 1 | 11s | 11s | 0.87 GB | 2c / 8 GB |
| `ATAC_DESCRIPTIVE_REPORT` | 1 | 11s | 11s | 0.74 GB | 4c / 16 GB |
| `ATAC_FINAL:MERGE_ANNOTATIONS` | 1 | 10s | 10s | 0.66 GB | 4c / 16 GB |
| `ATAC_FINAL:ATAC_CELLTYPE_ANNOTATION` | 1 | 2s | 2s | 0.06 GB | 4c / 16 GB |
| `ATAC_INITIAL:ATAC_MAKE_THRESHOLDS` | 1 | 2s | 2s | 0.06 GB | 4c / 16 GB |

TOTAL tasks 94 | summed realtime 3h07m | peak single-task RSS 8.20 GB

## Output footprint

| Stage | Size |
|---|---|
| `multiome/` | 3,124.6 MB |
| `hdwgcna/` | 247.3 MB |
| `rna/` | 169.5 MB |
| `atac/` | 144.8 MB |
| `mofa_visualization/` | 142.8 MB |
| `cell_annotation/` | 82.5 MB |
| everything else | < 30 MB |
| **TOTAL** | **3.85 GB** |

3.08 GB of that is three files:

- `multiome/multivi/multivi_integrated.h5mu` — 1.97 GB
- `multiome/multivi/multivi_model/model.pt` — 967 MB
- `multiome/mudata/integrated.h5mu` — 138 MB

None of this ships. The release asset is the 79 MB input dataset plus a small
metrics/figures bundle (Step 3b decision: `expected_results.json` + reference
PNGs, not the objects). This is a **local disk requirement** to document:
budget ~10 GB free including `work/`.

## TODO — RE-MEASURE (blocking for publication)

**The tier was fixed after this run; these numbers are from the BROKEN tier.**
Wall-clock, CPU-hours and scheduling behaviour all change with the corrected
allocations. Peak RSS per process is still valid — it measures actual demand,
not what was requested — and is what the new tier is sized from.

Before publishing any tutorial timing:

1. Re-run end to end on the fixed tier (fresh outdir, not `-resume`).
2. Confirm every process gets its intended allocation — the failure mode is
   silent, so check `trace.tsv`, do not assume.
3. Regenerate the tables in this file from the new `trace.tsv`.
4. Then do the run-2 reproducibility diff (Step 3 of PHASE3_HANDOFF.md).

Everything below documents the run that produced the fix.

## FIXED DEFECT — the tier's per-process resource blocks were mostly inert

**88 of 94 tasks ran at 2 CPU / 8 GB — the `withName: '.*'` wildcard's values —
not the 4 CPU / 16 GB the specific blocks declare.** Only 6 tasks got 16 GB.

Cause: Nextflow 25.10.0 selector precedence. A selector matching a process's
**fully-qualified** name outranks one matching only its **simple** name, and
declaration order is irrelevant. Verified in isolation:

| Config | Result |
|---|---|
| `withName: '.*'` then `withName: 'FOO'` | 2c / 8 GB — wildcard wins |
| `withName: 'FOO'` then `withName: '.*'` | 2c / 8 GB — wildcard wins *regardless of order* |
| `withName: '.*'` then `withName: 'SUB:FOO'` | 4c / 16 GB — specific wins |

`.*` matches the qualified name `SUB:FOO`; a bare `FOO` does not. So the
wildcard silently outranks every simple-name selector for any process nested in
a named workflow. That is why `ATAC_FINAL:MERGE_ANNOTATIONS` got 16 GB (its
qualified name matches `ATAC_.*`) while `RNA:CELLBENDER` got 8 GB.

Containers were unaffected: the `.*` block sets no `container`, so the
simple-name container selectors were uncontested.

**Production tiers are NOT affected** — `small`/`medium`/`large` contain no
`withName: '.*'` wildcard, so their simple-name selectors are the only match.
This is specific to `tutorial.config`.

### Why it mattered less than it could have

The local executor treats `memory` as advisory for scheduling, so nothing was
OOM-killed. Under SLURM, `MULTIVI_INTEGRATE` at **7.80 GB peak against an 8 GB
request (97.5%)** would be a coin-flip, and `ATAC_FINAL_PIPELINE` at 8.20 GB
would have been killed outright had it landed on the 8 GB default.

### Fix applied

Naively moving `cpus`/`memory` to bare `process` scope would have been WORSE:
`nextflow.config` assigns 8-500 GB via `withLabel` (`process_medium` 64 GB,
`process_high_memory` 256 GB, `hugemem` 500 GB), and `withLabel` outranks bare
process scope. A workstation run would have requested 500 GB, never scheduled,
and hung rather than failed.

What was actually done in `configs/resource_tiers/tutorial.config`:

1. Bare `process` scope holds a 2 CPU / 4 GB floor for unlabelled processes.
2. **All eight `withLabel` blocks are overridden** to tutorial-scale values —
   required, not optional, per the trap above.
3. `withName: '.*'` keeps ONLY `clusterOptions = ''` and `accelerator = null`.
   It must stay a wildcard to outrank nextflow.config's per-process GPU
   settings; it must not declare resources, or it shadows everything again.
4. Every specific selector is written `'(.*:)?NAME'` so it matches the
   qualified name too. Verified empirically that a broad `'ATAC_.*'` (which
   matches `ATAC_FINAL:ATAC_FINAL_PIPELINE`) otherwise beats a bare
   `'ATAC_FINAL_PIPELINE'` override and hands a 12 GB process 4 GB.
5. General buckets are declared BEFORE specific overrides, since with all
   selectors qualified-safe plain declaration order applies.

Right-sizing targets from measured peak RSS (~50% headroom):

| Process | Peak RSS | Suggested |
|---|---|---|
| `ATAC_FINAL_PIPELINE` | 8.20 GB | 12 GB |
| `MULTIVI_INTEGRATE` | 7.80 GB | 12 GB |
| `HDWGCNA_ENRICHMENT` | 7.70 GB | 12 GB |
| `HDWGCNA_PER_CELLTYPE` | 6.20 GB | 10 GB |
| `ATAC_INITIAL_QC` | 5.80 GB | 8 GB |
| `MULTIVI_VISUALIZE` | 4.90 GB | 8 GB |
| everything else | < 1.6 GB | 4 GB |

A corrected tier changes allocations, so **these wall-clock numbers must be
re-measured before publication**.

## Other findings from this run

- **`trace.tsv` / `report.html` / `timeline.html` land in `results/`, not
  `results_tutorial/`.** `nextflow.config` interpolates
  `${params.outdir}/pipeline_info/...` at parse time, before the `tutorial`
  profile sets `params.outdir`. Same class as the `scprinter.gtf_human` trap.
  Pre-flight prints the `logs/nextflow/trace.txt` path, which is a third
  location and also wrong.
- **`RUN_CELLTYPIST` downloads all 61 CellTypist models at runtime** from
  `celltypist.cog.sanger.ac.uk` to use one. A hard network dependency mid-run,
  and it contradicts the "nothing to download" framing. ~1m33s of the wall time.
- **The 45 CellTypist labels are not biology.** `Immune_All_Low.pkl` on 1,000
  cells yields `ELP`, `Double-negative_thymocytes`, `CD8a_a`,
  `Age-associated_B_cells` — thymocyte/progenitor classes impossible in PBMC.
  Retained deliberately: T2 is a wiring-and-plumbing demo, not a biological one.
  State this plainly in the tutorial.

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
2. Empirical precedence, via the throwaway harness pattern used for this fix —
   a two-process `.nf` plus a `trace` scope emitting `name,cpus,memory`. This is
   how the `'ATAC_.*'`-shadows-`'ATAC_FINAL_PIPELINE'` hazard was caught, and how
   the GPU-strip block was confirmed to outrank both competitors.
3. Ground truth, after a real run: read `trace.tsv` for cpus/memory per process,
   and grep the task's `.command.run` for `--nv` to confirm no GPU was attached.
   Only this exercises the actual merge for every process in the DAG.

Do (3) as part of the re-measure TODO above.
