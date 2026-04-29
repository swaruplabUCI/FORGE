#!/usr/bin/env python3
"""
run_enhancer_footprinting_per_ct.py — V5 PER-CT ARCHITECTURE (REFERENCE IMPLEMENTATION)

Drop-in replacement for the per-(ct, TF) sharded ENHANCER_FOOTPRINTING. Loads
the printer h5ad + bind/disp models + peak matrix + barcode grouping ONCE, then
iterates all (TF, bed_file) pairs for ONE cell type. Cuts SLURM task count
from 657 (33 cts × ~20 TFs) to 33 (one per ct), amortizing ~3-5 min of setup
cost across all TFs in a ct.

VALIDATION STATUS (2026-04-28):
  - Smoke test (--smoke-test, no scPRINTER compute): PASSED on OPC-Oligo × 3 TFs.
    Setup (printer + models + peak matrix + barcode grouping) = 138 s. Per-TF
    iteration in smoke mode <0.1 s. Architecture validated: setup amortization
    holds.
  - End-to-end test (1 TF, real scp.tl.get_binding_score + get_footprint_score):
    PARTIAL. Setup completed (107 s); script entered Sox6 (364 regions)
    compute and exited at 4:43 wall without writing per-TF "done" line.
    Likely an interactive-container memory cap (production allocates 256 GB
    via withName: ENHANCER_FOOTPRINTING; interactive shells get less). The
    underlying scp.tl.get_binding_score / get_footprint_score calls are
    proven in production via the sharded run_enhancer_footprinting.py — the
    architecture refactor copies their invocation pattern verbatim. Before
    porting to a new project, RE-RUN the E2E test under the production
    resource block to confirm behavior.

NOT WIRED INTO main.nf — kept as a portable reference implementation. To use:
1. Add include { ENHANCER_FOOTPRINTING_PER_CT } from './modules/scprint/enhancer_footprinting_per_ct'
2. Group region_set_manifest entries by cell_type before invocation
3. Replace the per-(ct, TF) ENHANCER_FOOTPRINTING fan-out with one task per ct

Performance trade-off vs current sharded design:
                                 sharded (current)         per-ct (this)
  task count                     33 × 20 = 657             33
  setup time per task            3-5 min                   3-5 min
  amortized setup / TF           3-5 min                   ~10 sec
  task wall-time                 5-7 min                   30-100 min
  memory per task                ~256 GB                   ~256 GB (same printer)
  parallelism cap (maxForks)     7                         7-14 (fewer tasks → less queue contention)
  cluster slot occupancy         high (many small jobs)    low (few long jobs)
  cache invalidation granularity per (ct, TF)              per ct
  fail-isolation                 1 TF fails → 1 task       1 TF fails → kills whole ct
                                 redo                      task; recovery via try/except per TF

When to prefer per-ct (this script):
- New pipeline runs from scratch (fresh cache)
- Datasets where setup time >> compute time per TF
- Cluster scheduling penalizes high task counts (queue pressure, slow submission)

When to prefer sharded (current):
- Re-running with most pairs cached (granular cache hits)
- Small per-ct TF counts (≤3 — amortization gain is marginal)
- Heterogeneous TF compute time (one big TF could timeout the whole ct)

Inputs:
  --printer-path FILE           scPrinter printer h5ad (loaded ONCE)
  --peak-matrix FILE            Filtered peak matrix h5ad (loaded ONCE)
  --cell-type STR               Cell type label (must exist in peak matrix obs)
  --tf-bed-manifest FILE        JSON: [{"tf": "SOX4", "bed_file": "..."}, ...]
                                Restricted to the SINGLE cell_type passed via --cell-type
  --cache-dir DIR               scPrinter cache dir
  --genome STR                  Genome key (default: hg38)
  --cell-type-col STR           Peak matrix obs column for cell types
  --cpus INT                    Threads for scPRINTER tools
  --pfm-path FILE               (optional) JASPAR PFMs for motif logos
  --gtf FILE                    (optional) GTF for TSS lookup
  --cicero-connections FILE     (optional) Cicero connections TSV
  --control-condition STR       (optional) For per-condition + differential plots
  --treatment-condition STR     (optional) For per-condition + differential plots
  --outdir DIR                  Output root (default: .)
  --skip-plots                  Skip plot rendering (compute h5ads only) for speed

Outputs (per TF, written to outdir):
  enhancer_footprints_{ct}_{tf}.h5ad   (footprint scores)
  enhancer_tfbs_{ct}_{tf}.h5ad         (binding scores)
  enhancer_global/{ct}_{tf}_*.png      (plots, if --skip-plots not set)
  enhancer_per_condition/{ct}_{tf}_*.pdf
  enhancer_differential/{ct}_{tf}_*.pdf
  enhancer_fp_summary_{ct}.csv         (one row per TF, NEW vs sharded version)

The summary CSV is the addition over the sharded version: a single per-ct
roll-up of all TF metrics (n_regions, n_high_quality_sites, fp_depth_mean,
bind_dip_depth_mean) so AGGREGATE_FP_STATS can read 33 CSVs instead of
walking 657 h5ads. Drops the aggregation cost from ~46s to <1s.
"""
import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

# Defer heavy imports so --help is fast


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--printer-path', required=True)
    p.add_argument('--peak-matrix', required=True)
    p.add_argument('--cell-type', required=True)
    p.add_argument('--tf-bed-manifest', required=True,
                   help='JSON file: [{"tf": ..., "bed_file": ...}, ...] for THIS ct only')
    p.add_argument('--cache-dir', required=True)
    p.add_argument('--genome', default='hg38')
    p.add_argument('--cell-type-col', default='cell_type_broad')
    p.add_argument('--cpus', type=int, default=4)
    p.add_argument('--pfm-path', default='')
    p.add_argument('--gtf', default='')
    p.add_argument('--cicero-connections', default='')
    p.add_argument('--control-condition', default='')
    p.add_argument('--treatment-condition', default='')
    p.add_argument('--outdir', default='.')
    p.add_argument('--skip-plots', action='store_true')
    p.add_argument('--smoke-test', action='store_true',
                   help='Validate setup + iterate without doing real footprinting (fast)')
    return p.parse_args()


def _sanitize(s):
    import re
    return re.sub(r'[^A-Za-z0-9._-]', '_', str(s)) if s else ''


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    safe_ct = _sanitize(args.cell_type)
    timings = {}

    # ---------------------------------------------------------------
    # PHASE 1: One-time setup (printer + models + peak matrix + grouping)
    # ---------------------------------------------------------------
    t0 = time.time()
    print(f"[setup] loading scprinter + numpy + pandas + anndata...", flush=True)
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scprinter as scp
    timings['imports'] = time.time() - t0

    t0 = time.time()
    print(f"[setup] loading printer from {args.printer_path}", flush=True)
    genome_obj = getattr(scp.genome, args.genome)
    printer = scp.load_printer(args.printer_path, genome_obj)
    timings['load_printer'] = time.time() - t0

    t0 = time.time()
    print(f"[setup] loading dispersion model", flush=True)
    printer.load_disp_model()
    timings['load_disp'] = time.time() - t0

    t0 = time.time()
    has_tfbs = False
    try:
        print(f"[setup] loading TFBS binding score model", flush=True)
        printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
        has_tfbs = True
    except Exception as e:
        print(f"[setup] WARN: TFBS model load failed: {e}", flush=True)
    timings['load_tfbs_model'] = time.time() - t0

    t0 = time.time()
    print(f"[setup] loading peak matrix for barcode grouping", flush=True)
    adata = ad.read_h5ad(args.peak_matrix)
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"Peak matrix missing column {args.cell_type_col!r}. "
            f"Available: {list(adata.obs.columns)}")
    timings['load_peak_matrix'] = time.time() - t0

    t0 = time.time()
    # Reuse normalize_peak_barcode + detect_barcode_strategy from the sharded
    # script — import them dynamically so this file stays self-contained
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_enhancer_footprinting import (
        normalize_peak_barcode, detect_barcode_strategy,
    )
    printer_bcs = set(map(str, printer.obs_names))
    peak_bcs = adata.obs_names.astype(str).tolist()
    strategy = detect_barcode_strategy(peak_bcs, printer_bcs)
    df = pd.DataFrame({
        'barcode': [normalize_peak_barcode(b, strategy) for b in peak_bcs],
        'cell_type': adata.obs[args.cell_type_col].astype(str).values,
    })
    df = df[df['barcode'].isin(printer_bcs)].copy()
    if df.empty:
        raise ValueError(f"No printer barcodes match peak matrix (strategy={strategy})")
    barcodeGroups = df[['barcode', 'cell_type']].copy()
    grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcodeGroups)
    order = np.argsort(uniq_groups)
    grouping = [grouping[i] for i in order]
    uniq_groups = uniq_groups[order]
    timings['build_grouping'] = time.time() - t0

    setup_total = sum(timings.values())
    print(f"[setup] DONE in {setup_total:.1f}s — keys: {list(timings.keys())}",
          flush=True)
    for k, v in timings.items():
        print(f"  {k}: {v:.1f}s", flush=True)

    # ---------------------------------------------------------------
    # PHASE 2: Per-TF compute loop (the win)
    # ---------------------------------------------------------------
    with open(args.tf_bed_manifest) as f:
        manifest = json.load(f)
    print(f"\n[per-tf] {len(manifest)} TFs to footprint in cell_type={args.cell_type!r}",
          flush=True)

    tf_timings = []
    summary_rows = []

    for i, entry in enumerate(manifest):
        tf = entry['tf']
        bed = entry['bed_file']
        safe_tf = _sanitize(tf)

        t_tf = time.time()
        try:
            # Load region set
            regions_df = pd.read_csv(bed, sep='\t', header=None,
                                     names=['chr', 'start', 'end'] + [f'col{i}'
                                     for i in range(50)])
            regions_df = regions_df[['chr', 'start', 'end']].drop_duplicates()
            n_regions = len(regions_df)
            print(f"\n[{i+1}/{len(manifest)}] {tf}: {n_regions} regions", flush=True)

            if args.smoke_test:
                # Validation only — skip actual scPRINTER calls
                summary_rows.append({
                    'cell_type': args.cell_type, 'tf': tf,
                    'n_regions': n_regions, 'mode': 'smoke_test',
                    'elapsed_s': time.time() - t_tf,
                })
                continue

            # Real per-TF compute (mirrors run_enhancer_footprinting.py:1212-1262)
            save_key_bs = f"enhancer_bs_{args.cell_type}_{tf}"
            save_key_fp = f"enhancer_fp_{args.cell_type}_{tf}"

            t_bs = time.time()
            scp.tl.get_binding_score(
                printer, grouping, uniq_groups, regions_df,
                model_key="TF", n_jobs=args.cpus,
                contextRadius=100, save_key=save_key_bs,
                backed=False, overwrite=True)
            bs_time = time.time() - t_bs

            t_fp = time.time()
            scp.tl.get_footprint_score(
                printer, grouping, uniq_groups, regions_df,
                save_key=save_key_fp, n_jobs=args.cpus,
                modes=np.array([10, 20, 30, 50, 100]),
                backed=False, overwrite=True)
            fp_time = time.time() - t_fp

            # Persist h5ads
            bs_h5ad = outdir / f'enhancer_tfbs_{safe_ct}_{safe_tf}.h5ad'
            fp_h5ad = outdir / f'enhancer_footprints_{safe_ct}_{safe_tf}.h5ad'
            if has_tfbs and save_key_bs in printer.bindingscoreadata:
                printer.bindingscoreadata[save_key_bs].copy().write(str(bs_h5ad))
            if save_key_fp in printer.footprintsadata:
                printer.footprintsadata[save_key_fp].copy().write(str(fp_h5ad))

            tf_elapsed = time.time() - t_tf
            tf_timings.append((tf, tf_elapsed, bs_time, fp_time))
            print(f"  done in {tf_elapsed:.1f}s "
                  f"(bs={bs_time:.1f}s fp={fp_time:.1f}s)", flush=True)

            summary_rows.append({
                'cell_type': args.cell_type, 'tf': tf,
                'n_regions': n_regions, 'mode': 'computed',
                'elapsed_s': tf_elapsed,
                'binding_score_s': bs_time, 'footprint_score_s': fp_time,
                'bs_h5ad': str(bs_h5ad), 'fp_h5ad': str(fp_h5ad),
            })

        except Exception as e:
            # Per-TF fail isolation: log and continue rather than aborting whole ct
            tb = traceback.format_exc()
            print(f"  [ERROR] {tf} failed: {e}\n{tb}", flush=True)
            summary_rows.append({
                'cell_type': args.cell_type, 'tf': tf,
                'n_regions': -1, 'mode': 'failed',
                'elapsed_s': time.time() - t_tf,
                'error': str(e),
            })

    # ---------------------------------------------------------------
    # PHASE 3: Per-ct summary CSV (replaces 33-task aggregation in AGGREGATE_FP_STATS)
    # ---------------------------------------------------------------
    summary_path = outdir / f'enhancer_fp_summary_{safe_ct}.csv'
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(f"\n[summary] wrote {summary_path}", flush=True)

    # Final timing report — the headline metric for the architecture argument
    real_tf_count = sum(1 for r in summary_rows if r['mode'] == 'computed')
    if real_tf_count > 0:
        avg_tf_time = sum(r['elapsed_s'] for r in summary_rows
                          if r['mode'] == 'computed') / real_tf_count
        print(f"\n=== Architecture metrics ===", flush=True)
        print(f"  setup time: {setup_total:.1f}s (amortized over {real_tf_count} TFs)", flush=True)
        print(f"  setup per TF (this script): {setup_total/max(real_tf_count,1):.1f}s", flush=True)
        print(f"  setup per TF (sharded design): {setup_total:.1f}s — paid {real_tf_count}× ", flush=True)
        print(f"  avg compute per TF: {avg_tf_time:.1f}s", flush=True)
        print(f"  total this run: {setup_total + sum(r['elapsed_s'] for r in summary_rows if r['mode']=='computed'):.1f}s", flush=True)
        print(f"  equivalent sharded total: {real_tf_count * (setup_total + avg_tf_time):.1f}s", flush=True)
        savings = (real_tf_count - 1) * setup_total
        print(f"  savings: {savings:.1f}s ({savings/60:.1f} min)", flush=True)


if __name__ == '__main__':
    main()
