# Refactor5 Pipeline: Prioritized Breakage Risk List (Round 1)

**Date:** 2026-03-25
**Reviewer:** Senior Bioinformatics Engineering Review
**Branch:** dev (commit 9877678)

## Ranking Criteria

Each issue is scored on four dimensions and then sorted by composite severity:

1. **Likelihood** -- How easily can a real user (especially a new collaborator) trigger this?
2. **Blast radius** -- Silent data corruption (worst) > wrong results with no warning > loud crash (least bad)
3. **Fix complexity** -- One-line config change < small refactor < architectural change
4. **Who it affects** -- Only the original dev < any user on this HPC < any user on any HPC < published results

## Priority Tier Definitions

| Tier | Meaning |
|------|---------|
| P0   | Fix before any production run -- pipeline will fail or silently corrupt results |
| P1   | Fix before sharing with any collaborator -- blocks portability or reproducibility |
| P2   | Fix before publication or release -- correctness edge cases, hardening |
| P3   | Nice-to-have cleanup -- tech debt, defensive coding, minor UX |

## Ranked Issues

| Rank | Orig # | Tier | Issue | Rationale |
|------|--------|------|-------|-----------|
| 1 | 2 | P0 | Global `errorStrategy = 'ignore'` silently swallows critical failures | The cluster profile sets the default to `ignore`, meaning any core process (CONCAT_BATCHES, RNA_QC, ATAC_FINAL) that fails will be silently skipped; downstream steps receive empty channels and produce incomplete results with no error -- the worst possible failure mode for a scientific pipeline. |
| 2 | 18 | P0 | Reading from publishDir instead of process output channels (race condition) | If any process reads files from publishDir rather than Nextflow output channels, results depend on filesystem timing; a slow NFS can produce silent data loss or stale inputs that are impossible to debug. |
| 3 | 12 | P0 | `ad.concat(index_unique=None)` can silently merge barcodes across samples | Found in 3 scripts including the critical concat_batches.py; duplicate barcodes across samples will be silently collapsed, corrupting all downstream per-cell analyses with no warning. |
| 4 | 1 | P0 | Undefined legacy params (`june_metadata`, `july_souporcell_dir`) still referenced in main.nf | RNA_QC is called with these on every run (line 221-223); for any non-BD-brain dataset they resolve to `NO_FILE` sentinels, which only works if the QC script correctly ignores them -- a fragile implicit contract that will break on script changes. |
| 5 | 10 | P0 | Fragile `NO_FILE` sentinel pattern across main.nf | Used in 15+ locations with inconsistent sentinel names (`NO_FILE`, `NO_FILE_JULY`, `NO_FILE_dorc`, etc.); any process that does not explicitly check for these strings will attempt to open a nonexistent file, and `errorStrategy = 'ignore'` will hide the crash. |
| 6 | 5 | P1 | User-specific conda path in launch.sh (`/pub/lesolano/miniconda3`) | Any collaborator who clones the repo will get an immediate launch failure; trivial to fix but guaranteed to hit every new user on first run. |
| 7 | 24 | P1 | SLURM account hardcoded to `vswarup_lab` | Every clusterOptions string uses `params.slurm_account` which defaults to `vswarup_lab`; any external collaborator's jobs will be rejected by the scheduler immediately. |
| 8 | 4 | P1 | Hardcoded `/dfs7` absolute paths throughout all configs | Binds the entire pipeline to UCI HPC3 filesystem layout; affects launch.sh, singularity runOptions, and tools PATH -- no collaborator at a different institution can run this without editing multiple files. |
| 9 | 14 | P1 | Hardcoded SLURM partition names (`highmem`, `hugemem`, `maxmem`) | Found in 20+ clusterOptions lines in both resource tier configs; these are UCI-specific partition names that will cause job submission failures at any other site. |
| 10 | 6 | P1 | Missing `medium` resource tier (auto-detected but no config exists) | launch.sh assigns `medium` for 6-50 samples but nextflow.config only loads `small` or `large`; 6-50 sample datasets silently fall back to `small` tier, causing OOM kills on moderately sized data. |
| 11 | 23 | P1 | `pycistopic.gtf = null` default but required when enabled | pycisTopic is enabled by default (`run = true`) but its required GTF is null; the process will fail at runtime with an unhelpful null-pointer error rather than a clear config validation message. |
| 12 | 7 | P1 | Container name mismatch between launch.sh check and actual config | launch.sh checks for `cicero.sif` but config defines `r_cicero` pointing to `cicero.sif`; if container naming drifts, the pre-flight check passes but the process fails. |
| 13 | 16 | P1 | GPU profile applies to ALL processes including CPU-only ones | The `gpu` profile sets `clusterOptions` with `--gres=gpu` at the process level globally; CPU-only processes like RNA_QC and CONCAT_BATCHES will request (and waste) GPU allocations, increasing queue wait times and burning GPU-hour credits. |
| 14 | 22 | P1 | `SINGULARITY_BINDPATH` in launch.sh has no effect with `--contain` | The singularity profile uses `--contain` which ignores the environment variable; the actual bind mounts come from `--bind /dfs7 --bind /tmp` in runOptions, making the export misleading and a trap for anyone trying to add new bind paths. |
| 15 | 13 | P2 | `/tmp` as TMPDIR on shared HPC (disk exhaustion) | Both launch.sh and config set TMPDIR=/tmp; large intermediate files (CellBender, scVI checkpoints) can fill the shared /tmp partition and crash other users' jobs on the same node. |
| 16 | 15 | P2 | Hardcoded MAST DE comparisons (Group1-4 only) | The RNA_DIFFERENTIAL workflow hardcodes 6 pairwise comparisons between Group1-Group4 (lines 479-486) instead of reading from `params.differential_rna.comparisons`; any real experiment with different group names requires editing main.nf. |
| 17 | 9 | P2 | Only human MT genes detected (`MT-` hardcoded); mouse `mt-` ignored in concat_batches.py | concat_batches.py only checks `MT-` (line 99); mouse datasets will have 0% mito genes, disabling mito QC filtering and passing low-quality cells through. rna_qc.py handles both species correctly, creating an inconsistency. |
| 18 | 3 | P2 | Hardcoded batch detection (`june`/`july`/`nov` string matching) breaks non-matching datasets | The legacy params (june_metadata, etc.) are still wired into RNA_QC; while the generic batch_dirs map is the intended path forward, the old code path remains as dead-but-executed code that confuses new developers. |
| 19 | 11 | P2 | MOFA gets GPU container but CPU partition in small tier | MOFA_INTEGRATE is listed in the GPU accelerator block (line 423) but small.config may route it to a CPU partition; the job will either fail to find a GPU or waste a GPU slot depending on profile precedence. |
| 20 | 17 | P2 | `pathway_to_tfs.json` not validated before use | The CellChat-to-TF hypothesis script loads this JSON without schema checks; a malformed or empty file produces silent empty-output rather than a clear error. |
| 21 | 25 | P2 | No CSV manifest schema validation | The manifest CSV is the source of truth for the entire pipeline; a missing column (e.g., `sample_type`, `batch`, `fragment_file`) causes cryptic Groovy null-pointer errors deep in channel operations. |
| 22 | 8 | P2 | `eval "$NF_CMD"` with user-controlled EXTRA_ARGS (shell injection risk) | EXTRA_ARGS comes from positional arguments to launch.sh; on a shared HPC with sbatch this is low-risk in practice since the user is already authenticated, but it violates defense-in-depth principles. |
| 23 | 21 | P3 | hdWGCNA_DIFFERENTIAL channel shape / filename parsing fragility | The `.replaceAll(/^hdwgcna_/, '')` pattern (line 383) to extract cell type from filename is brittle; cell types with underscores or special characters will be mangled. |
| 24 | 19 | P3 | `builtins.long = int` monkey-patch fragility | Present in 3 scPRINTer scripts; works around a Python 3 incompatibility in an upstream library but will break silently if the library updates its internal checks. |
| 25 | 20 | P3 | Duplicate process definitions across config files | Container and resource assignments are duplicated between nextflow.config, small.config, and large.config (documented as FIX-51 necessity); increases maintenance burden but is a known Nextflow limitation, not a bug. |

## Summary by Tier

| Tier | Count | Action |
|------|-------|--------|
| P0 | 5 | Fix immediately -- these can silently corrupt results or mask failures |
| P1 | 9 | Fix before any collaborator runs the pipeline -- portability blockers |
| P2 | 8 | Fix before publication -- correctness edge cases and hardening |
| P3 | 3 | Cleanup when convenient -- tech debt with low blast radius |
