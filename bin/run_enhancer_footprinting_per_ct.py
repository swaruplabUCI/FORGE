#!/usr/bin/env python3
"""
run_enhancer_footprinting_per_ct.py — V6 PER-CT ARCHITECTURE

Drop-in replacement for the per-(ct, TF) sharded ENHANCER_FOOTPRINTING. Loads
the printer h5ad + bind/disp models + peak matrix + barcode grouping ONCE per
cell type, then iterates all (TF, bed_file) pairs by delegating each pair to
run_enhancer_footprinting.run_one_pair() — the same code path the sharded CLI
uses. Cuts SLURM task count from 657 (33 cts × ~20 TFs) to 33 (one per ct),
amortizing ~3-5 min of setup cost across all TFs in a ct.

V6 (2026-04-30): No longer rolls its own compute loop. Calls into
run_one_pair() so per-ct emits the SAME per-(ct, TF) artifacts the legacy
sharded path produces:
  enhancer_global/{ct}_{tf}_*.png        (composite + intermediate plots)
  enhancer_per_condition/{ct}_{tf}_*.pdf
  enhancer_differential/{ct}_{tf}_*.pdf
  enhancer_footprints_{ct}_{tf}.h5ad
  enhancer_tfbs_{ct}_{tf}.h5ad
  enhancer_fp_summary_{ct}.csv          (per-ct roll-up across TFs)

The composite PNGs are required by COMPOSITE_ENHANCER_VIZ — its
find_footprint_pngs() does a recursive glob for *_composite.png and embeds
them as insets in the per-(gene, TF, ct) browser composites.

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
  --smoke-test                  Validate setup + iterate without doing real footprinting (fast)
"""
import argparse
import json
import os
import re
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace


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
    p.add_argument('--strip-target-genes', default='',
                   help='Comma-separated gene names; regions near their TSS get full tensor in obsm.')
    p.add_argument('--strip-context-bp', type=int, default=500000,
                   help='Context window (bp) around each strip target gene TSS.')
    p.add_argument('--binding-threshold', default='otsu',
                   help='Passed through to run_one_pair (otsu|percentile_75|float)')
    p.add_argument('--smoke-test', action='store_true',
                   help='Validate setup + iterate without doing real footprinting (fast)')
    return p.parse_args()


def _sanitize(s):
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
    print(f"[setup] importing heavy deps + sharded footprinting module...", flush=True)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import pandas as pd
    import run_enhancer_footprinting as ref
    timings['imports'] = time.time() - t0

    t0 = time.time()
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir
    ref.init_scprinter(args.cache_dir)
    state = ref.setup_pipeline(args)
    timings['setup_pipeline'] = time.time() - t0

    setup_total = sum(timings.values())
    print(f"[setup] DONE in {setup_total:.1f}s — keys: {list(timings.keys())}",
          flush=True)
    for k, v in timings.items():
        print(f"  {k}: {v:.1f}s", flush=True)

    # ---------------------------------------------------------------
    # PHASE 2: Plot output dirs (per-ct, shared across all TFs)
    # ---------------------------------------------------------------
    plots_dir, per_condition_dir, differential_dir = ref._make_plot_dirs(outdir)

    # ---------------------------------------------------------------
    # PHASE 3: Per-TF compute loop — delegates to run_one_pair()
    # ---------------------------------------------------------------
    with open(args.tf_bed_manifest) as f:
        manifest = json.load(f)
    print(f"\n[per-tf] {len(manifest)} TFs to footprint in cell_type={args.cell_type!r}",
          flush=True)

    summary_rows = []

    try:
        for i, entry in enumerate(manifest):
            tf = entry['tf']
            bed = entry['bed_file']
            # Nextflow stages bed files under tf_bed_files/ (see module input);
            # manifest carries bare basenames. Try as-is first (smoke-test
            # layout), then the staged subdir.
            if not os.path.exists(bed):
                staged = os.path.join('tf_bed_files', os.path.basename(bed))
                if os.path.exists(staged):
                    bed = staged

            t_tf = time.time()
            print(f"\n[{i+1}/{len(manifest)}] {tf}  bed={bed}", flush=True)

            if args.smoke_test:
                summary_rows.append({
                    'cell_type': args.cell_type, 'tf': tf,
                    'mode': 'smoke_test',
                    'elapsed_s': time.time() - t_tf,
                })
                continue

            # Build per-pair args namespace for run_one_pair. Field set must
            # match what the sharded CLI parser produces (see parse_args() in
            # run_enhancer_footprinting.py).
            pair_args = SimpleNamespace(
                cell_type=args.cell_type,
                tf_name=tf,
                region_set=bed,
                gtf=args.gtf,
                pfm_path=args.pfm_path,
                cicero_connections=args.cicero_connections,
                control_condition=args.control_condition,
                treatment_condition=args.treatment_condition,
                binding_threshold=args.binding_threshold,
                cpus=args.cpus,
                strip_target_genes=getattr(args, 'strip_target_genes', ''),
                strip_context_bp=getattr(args, 'strip_context_bp', 500000),
            )

            try:
                row = ref.run_one_pair(
                    pair_args,
                    printer=state['printer'],
                    peak_matrix_path=state['peak_matrix_path'],
                    grouping=state['grouping'],
                    uniq_groups=state['uniq_groups'],
                    ct_col=state['ct_col'],
                    has_tfbs=state['has_tfbs'],
                    plots_dir=plots_dir,
                    per_condition_dir=per_condition_dir,
                    differential_dir=differential_dir,
                )
                tf_elapsed = time.time() - t_tf
                # If run_one_pair returned but footprint_saved is False, the
                # multiscale footprint computation hit a non-fatal error
                # (caught by inner try/except in run_enhancer_footprinting.py).
                # Reflect that in mode so the summary CSV doesn't claim success.
                final_mode = 'computed' if row.get('footprint_saved') else 'partial'
                row.update({'mode': final_mode, 'elapsed_s': tf_elapsed})
                summary_rows.append(row)
                print(f"  done in {tf_elapsed:.1f}s (mode={final_mode})", flush=True)
            except BaseException as e:
                # pyo3 PanicException inherits BaseException (not Exception);
                # `except Exception` lets it through and kills the whole ct
                # task (run 28 lost 14h of work this way at BACH1 #17).
                # Allow KI/SE/GE through so ^C / sys.exit() still work.
                if isinstance(e, (KeyboardInterrupt, SystemExit, GeneratorExit)):
                    raise
                # Per-TF fail isolation: log and continue rather than aborting whole ct
                tb = traceback.format_exc()
                etype = type(e).__name__
                print(f"  [ERROR] {tf} failed ({etype}): {e}\n{tb}", flush=True)
                summary_rows.append({
                    'cell_type': args.cell_type, 'tf': tf,
                    'mode': 'failed',
                    'elapsed_s': time.time() - t_tf,
                    'error_type': etype,
                    'error': str(e),
                })

        # ------------------------------------------------------------
        # PHASE 4: Per-ct summary CSV (drives AGGREGATE_FP_STATS)
        # ------------------------------------------------------------
        summary_path = outdir / f'enhancer_fp_summary_{safe_ct}.csv'
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
        print(f"\n[summary] wrote {summary_path}", flush=True)

        # ---- Architecture metrics ----
        real_tf_count = sum(1 for r in summary_rows if r.get('mode') == 'computed')
        if real_tf_count > 0:
            avg_tf_time = sum(r['elapsed_s'] for r in summary_rows
                              if r.get('mode') == 'computed') / real_tf_count
            total = setup_total + sum(r['elapsed_s'] for r in summary_rows
                                      if r.get('mode') == 'computed')
            print(f"\n=== Architecture metrics ===", flush=True)
            print(f"  setup time: {setup_total:.1f}s (amortized over {real_tf_count} TFs)", flush=True)
            print(f"  setup per TF (this script): {setup_total/real_tf_count:.1f}s", flush=True)
            print(f"  setup per TF (sharded design): {setup_total:.1f}s — paid {real_tf_count}×", flush=True)
            print(f"  avg compute per TF: {avg_tf_time:.1f}s", flush=True)
            print(f"  total this run: {total:.1f}s", flush=True)
            print(f"  equivalent sharded total: {real_tf_count * (setup_total + avg_tf_time):.1f}s", flush=True)
            savings = (real_tf_count - 1) * setup_total
            print(f"  savings: {savings:.1f}s ({savings/60:.1f} min)", flush=True)

    finally:
        # Always close the printer to release file handles
        if hasattr(state['printer'], 'close'):
            try:
                state['printer'].close()
            except Exception:
                pass


if __name__ == '__main__':
    main()
