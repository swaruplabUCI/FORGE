#!/usr/bin/env python3
"""
build_ct_annotation_v2.py  (FORGE generalized version)

Produce ct_annotation.csv: a barcode-level cell_type × condition mapping
that drives CICERO_TRIPLETS_PER_CT (per-CT × condition Cicero reanalysis).

Works with any FORGE h5ad that carries:
  - a cell-type obs column  (--cell-type-col, default: cell_type_broad)
  - a condition obs column  (--condition-col,  default: condition)

No tissue-specific annotation hierarchy is assumed. Users pre-select the
desired resolution by choosing --cell-type-col:
  - cell_type_broad         (always present, FORGE-canonical broad groupings)
  - celltypist_prediction   (fine CellTypist labels)
  - cell_type_prediction    (scATAnno labels)
  - any other obs column

Abundance filter (applied globally, per cell type across all conditions):
  threshold = max(--min-cells, int(--min-pct x total_cells))
  CTs below threshold -> relabeled EXCLUDED in output (skipped by Nextflow).

Explicit exclude list (--exclude-labels):
  CTs matching these labels -> relabeled EXCLUDED regardless of abundance.

Output columns:
  barcode      -- full cell barcode (matches h5ad obs index)
  condition    -- condition value
  cell_type_v2 -- cell type label, or 'EXCLUDED' if below floor / in exclude list

Usage (typical FORGE run, broad resolution):
  python build_ct_annotation_v2.py \\
      --input  results/atac/final/peak_matrix_annotated.h5ad \\
      --output results/cicero/ct_annotation.csv

Usage (fine-grained, custom filter):
  python build_ct_annotation_v2.py \\
      --input         peak_matrix_annotated.h5ad \\
      --output        ct_annotation.csv \\
      --cell-type-col celltypist_prediction \\
      --condition-col condition \\
      --min-pct       0.01 \\
      --min-cells     100 \\
      --exclude-labels "Unknown,Other,Low_Quality"
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

EXCLUDED_LABEL = "EXCLUDED"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",    required=True,
                   help="Path to annotated peak matrix h5ad")
    p.add_argument("--output",   default="ct_annotation.csv",
                   help="Output CSV path (default: ct_annotation.csv)")
    p.add_argument("--cell-type-col", default="cell_type_broad",
                   help="obs column for cell types (default: cell_type_broad)")
    p.add_argument("--condition-col", default="condition",
                   help="obs column for condition (default: condition)")
    p.add_argument("--min-pct",  type=float, default=0.0,
                   help="Min fraction of total cells for a CT to qualify (default: 0)")
    p.add_argument("--min-cells", type=int,  default=0,
                   help="Min absolute cell count for a CT to qualify (default: 0)")
    p.add_argument("--exclude-labels", default="",
                   help="Comma-separated cell-type labels to always exclude "
                        "(e.g. 'Unknown,Other,Low_Quality')")
    return p.parse_args()


def main():
    args = parse_args()

    exclude_set = {lbl.strip() for lbl in args.exclude_labels.split(",") if lbl.strip()}

    print(f"[load] {args.input}")
    adata = ad.read_h5ad(args.input, backed="r")
    total = adata.n_obs
    print(f"[info] {total:,} cells loaded")

    for col in (args.cell_type_col, args.condition_col):
        if col not in adata.obs.columns:
            sys.exit(f"ERROR: column '{col}' not found in obs "
                     f"(available: {sorted(adata.obs.columns.tolist())})")

    obs = adata.obs[[args.cell_type_col, args.condition_col]].copy()
    obs.index.name = "barcode"

    # ── abundance filter ─────────────────────────────────────────────────────
    threshold = max(args.min_cells, int(args.min_pct * total)) if (args.min_pct > 0 or args.min_cells > 0) else 0
    ct_counts = obs[args.cell_type_col].value_counts()

    excluded_floor = set()
    if threshold > 0:
        excluded_floor = set(ct_counts[ct_counts < threshold].index.tolist())

    excluded_all = exclude_set | excluded_floor

    # ── assign cell_type_v2 ──────────────────────────────────────────────────
    obs["cell_type_v2"] = obs[args.cell_type_col].astype(str)
    obs.loc[obs["cell_type_v2"].isin(excluded_all), "cell_type_v2"] = EXCLUDED_LABEL

    # ── audit ────────────────────────────────────────────────────────────────
    print(f"\n[audit] Strata (cell_type_v2 x {args.condition_col})")
    print("-" * 70)
    ct_cond = (
        obs.groupby(["cell_type_v2", args.condition_col], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    ct_cond["total"] = ct_cond.sum(axis=1)
    ct_cond["pct"]   = (ct_cond["total"] / total * 100).round(2)
    ct_cond = ct_cond.sort_values("total", ascending=False)

    conditions = [c for c in ct_cond.columns if c not in ("total", "pct")]
    header = f"  {'Cell type':<40}" + "".join(f"{c:>10}" for c in conditions)
    header += f"{'Total':>8} {'%':>6}  {'Cicero?'}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    n_runnable = 0
    for ct, row in ct_cond.iterrows():
        cond_vals = "".join(f"{int(row[c]):>10}" for c in conditions)
        if ct == EXCLUDED_LABEL:
            cicero_str = "no (excluded)"
        else:
            n_runnable += 1
            min_cond = min(int(row[c]) for c in conditions)
            warn = "  !! <250" if min_cond < 250 else ""
            cicero_str = f"yes{warn}"
        print(f"  {ct:<40}{cond_vals}{int(row['total']):>8} {row['pct']:>5.2f}%  {cicero_str}")

    n_jobs = n_runnable * len(conditions)
    excluded_total = int(ct_cond.loc[EXCLUDED_LABEL, "total"]) if EXCLUDED_LABEL in ct_cond.index else 0
    print(f"\n  Runnable strata: {n_runnable} x {len(conditions)} conditions = {n_jobs} Cicero jobs")
    if excluded_total:
        print(f"  Excluded (EXCLUDED): {excluded_total:,} cells "
              f"({100*excluded_total/total:.2f}%)")
    if excluded_floor:
        print(f"  Below floor (threshold={threshold}): {sorted(excluded_floor)}")
    if exclude_set:
        print(f"  Explicit exclude list: {sorted(exclude_set)}")

    # ── write ────────────────────────────────────────────────────────────────
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    result = pd.DataFrame({
        "barcode":      obs.index,
        "condition":    obs[args.condition_col].values,
        "cell_type_v2": obs["cell_type_v2"].values,
    })
    result.to_csv(out, index=False)
    print(f"\n[write] {out}  ({len(result):,} rows, {n_runnable} runnable cell types)")


if __name__ == "__main__":
    main()
