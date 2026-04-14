# Adversarial Audit Round 3 — Post-R4 QoL Port
**Date**: 2026-03-30
**Scope**: Full pipeline review after R4 QoL port (commit c4d08fc)
**Auditors**: Useful Idiot (UI) + Senior Dev (SD) — dual-persona parallel audit

---

## Merged & Deduplicated Issue List

Issues from both auditors merged, deduplicated, and re-ranked by the Senior Dev.

### P1 — Crash / Wrong Results (8 issues)

| ID | Title | Files | Description |
|----|-------|-------|-------------|
| **R3-1** | SCENICPLUS_GRN_VIZ missing container assignment | nextflow.config ~L406 | Process not listed in `withName` block for scenicplus container. Will run without container or fail. (UI-1, SD-1) |
| **R3-2** | Hardcoded Python path in GRN viz module | modules/multiome/scenicplus_grn_viz.nf:31 | `/opt/miniforge3/envs/gt/bin/python` assumes specific conda env path inside container. No fallback. (UI-2, SD-2) |
| **R3-3** | GTF species not conditional in EXTRACT_CCAN_ENHANCERS call | main.nf:1770 | Always passes `params.scprinter.gtf_human` regardless of `params.species`. Mouse datasets get human GTF. (UI-7, UI-8) |
| **R3-4** | params.scprinter.gtf_human may resolve to "null" string | nextflow.config ~L319 | `gtf_human = "${params.gtf_human_full}"` embeds null at parse time → string "null" passed to scripts. (UI-20, UI-9) |
| **R3-5** | Enhancer footprinting Phase 2 variables may be undefined | bin/run_enhancer_footprinting.py ~L1759-1835 | `agg_msfp`, `agg_bs_mean`, etc. initialized to None but referenced in plot code without None-checks when Phase 2 conditions aren't met. (UI-28, UI-29) |
| **R3-6** | plot_postscanvi.py Path B fallback to nonexistent key | bin/plot_postscanvi.py | If no cell_type_key auto-detected, falls back to 'scanvi_prediction' which doesn't exist in CellTypist-only (Path B) runs. (UI-15) |
| **R3-7** | mofa_vis.py cell_type_source uninitialized before loop | bin/mofa_vis.py:274 | Variable set inside for loop; if loop doesn't break, NameError at line 362 when building UMAP title. (UI-5, SD-15) |
| **R3-8** | Missing container assignments for multiple processes | nextflow.config | MOFA_VISUALIZE, CONVERT_CELLISMO, and ~10 other processes missing from `withName` container blocks. (SD-8, SD-10, SD-18) |

### P2 — Degraded Output / Wasted Compute (12 issues)

| ID | Title | Files | Description |
|----|-------|-------|-------------|
| **R3-9** | tissue_type config path mismatch | nextflow.config:187, plot_post_scanvi.nf:19 | `tissue_type` defined in `params.celltypist` block but module references `params.celltypist.tissue_type` — verify actual location matches. (UI-6, UI-26) |
| **R3-10** | Undefined resource labels (hugemem, process_gpu, process_high_memory, process_hdwgcna, process_small) | nextflow.config, various modules | 6+ processes use labels not defined in config. Default minimal resources → OOM on large data. (SD-3 through SD-7) |
| **R3-11** | Enhancer footprinting Phase 2 args not wired through .nf module | modules/scprint/enhancer_footprinting.nf | `--binding-threshold`, `--control-condition`, `--treatment-condition` not passed from module to script. Phase 2 differential mode unreachable. (SD-13) |
| **R3-12** | SCENICPLUS_GRN_VIZ outputs not captured in workflow emit | main.nf ~L1710-1735 | GRN viz pdfs/pngs created but not emitted from MULTIOME_GRN workflow. (UI-11, UI-14, SD-20) |
| **R3-13** | condition_label empty string produces odd filenames | main.nf:1771, extract_ccan_enhancers.py | `''` condition_label → output files with no condition suffix. Works but confusing naming. (UI-10) |
| **R3-14** | Cicero graph_test() may not exist in container | bin/run_cicero_full.R:361 | `graph_test()` and `plot_accessibility_in_pseudotime()` availability depends on cicero container version. Wrapped in tryCatch but silently skips. (SD-12) |
| **R3-15** | hdWGCNA cell type column inconsistency with broad mapping | bin/plot_postscanvi.py:566-578 | Adaptive broad mapping creates `cell_type_broad` column that downstream hdWGCNA may not expect. (SD-14) |
| **R3-16** | differential_rna.cell_type_key hardcoded mismatch | nextflow.config:118 | Always 'celltypist_prediction' even in Path A where 'scanvi_prediction' is primary. (UI-19) |
| **R3-17** | ENHANCER_FOOTPRINTING_RECIPES runs on empty channels if ATAC disabled | main.nf ~L1760-1800 | If `params.atac.run = false` but `enhancer_footprinting.run = true`, recipes run on empty inputs. (UI-17) |
| **R3-18** | No empty-channel guard for Cicero outputs | main.nf ~L1767-1772 | EXTRACT_CCAN_ENHANCERS assumes non-empty cicero_conns_ch / cicero_ccan_ch. (UI-24) |
| **R3-19** | HDF5 locking env var set redundantly | nextflow.config, various modules | Both global singularity.runOptions and per-process script blocks set HDF5_USE_FILE_LOCKING=FALSE. Inconsistent but harmless. (UI-21) |
| **R3-20** | Missing BD-specific param definitions | nextflow.config vs main.nf | `params.june_metadata`, `params.july_souporcell_dir`, `params.nov_souporcell_dir` used in main.nf but not defined in nextflow.config. (SD-11) |

### P3 — Cosmetic / Code Smell (2 issues)

| ID | Title | Files | Description |
|----|-------|-------|-------------|
| **R3-21** | Deprecated demux_*_dir params still in config | nextflow.config:176-181 | Dead parameters from BD-specific legacy code. Confusing but harmless. (SD-17) |
| **R3-22** | SCENICPLUS_GRN_VIZ no errorStrategy | modules/multiome/scenicplus_grn_viz.nf | Missing `errorStrategy 'retry'`. Single failure blocks pipeline vs graceful degradation. (UI-25) |

---

## Summary

| Severity | Count |
|----------|-------|
| P1 (Crash/Wrong Results) | 8 |
| P2 (Degraded Output) | 12 |
| P3 (Cosmetic) | 2 |
| **Total** | **22** |

## Recommended Fix Order

**Immediate (blocks execution):**
1. R3-1 — Add SCENICPLUS_GRN_VIZ to container withName block
2. R3-8 — Audit all processes for container assignments
3. R3-2 — Validate/fix Python path in GRN viz module

**Pre-production (data integrity):**
4. R3-3 — Species-conditional GTF in EXTRACT_CCAN_ENHANCERS
5. R3-4 — Fix scprinter.gtf_human/mouse null-to-string issue
6. R3-7 — Initialize cell_type_source before loop in mofa_vis.py
7. R3-5 — Add None-checks for Phase 2 variables in enhancer footprinting
8. R3-6 — Fix Path B fallback in plot_postscanvi.py

**Runtime stability:**
9. R3-10 — Define all missing resource labels
10. R3-11 — Wire Phase 2 enhancer args through module
11. R3-9 — Verify tissue_type config path
12. R3-16 — Make differential_rna.cell_type_key dynamic
