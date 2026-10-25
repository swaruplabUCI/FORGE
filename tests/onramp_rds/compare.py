#!/usr/bin/env python3
"""Diff the Seurat-side and h5ad-side artifact sets and report pass/fail.

Usage:
    compare.py --ref <seurat_outdir> --test <h5ad_outdir> [--tol 1e-6]

Exit code 0 if every check passes, 1 otherwise, so this can gate CI.

Tolerances are deliberately tight. The converter does not resample, rescale or
reorder anything -- the sparse slots are reused verbatim -- so the only
expected difference is float64 round-tripping through HDF5, which is exact.
A loose tolerance here would let a real regression through. The one check that
is allowed to be approximate is the per-(cell_type, gene) dot-plot block,
where R and numpy accumulate a mean in a different order.
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--ref", required=True, help="summarize_seurat.R outdir")
p.add_argument("--test", required=True, help="summarize_h5ad.py outdir")
p.add_argument("--tol", type=float, default=1e-6)
p.add_argument("--label", default="")
args = p.parse_args()

results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))


def load_json(d):
    with open(os.path.join(d, "summary.json")) as fh:
        return json.load(fh)


def rel_err(a, b):
    a, b = float(a), float(b)
    if a == b:
        return 0.0
    denom = max(abs(a), abs(b), 1e-30)
    return abs(a - b) / denom


ref = load_json(args.ref)
tst = load_json(args.test)

# ---- shape and totals ----
check("n_cells", ref["n_cells"] == tst["n_cells"], f"{ref['n_cells']} vs {tst['n_cells']}")
check("n_genes", ref["n_genes"] == tst["n_genes"], f"{ref['n_genes']} vs {tst['n_genes']}")
check("X nnz", ref["x_nnz"] == tst["x_nnz"], f"{ref['x_nnz']} vs {tst['x_nnz']}")
e = rel_err(ref["x_sum"], tst["x_sum"])
check("X sum", e <= args.tol, f"{ref['x_sum']:.6g} vs {tst['x_sum']:.6g} (rel {e:.2e})")

check("lognorm present", ref["has_lognorm"] == tst["has_lognorm"],
      f"{ref['has_lognorm']} vs {tst['has_lognorm']}")
if ref["has_lognorm"] and tst["has_lognorm"]:
    e = rel_err(ref["lognorm_sum"], tst["lognorm_sum"])
    check("lognorm sum", e <= args.tol,
          f"{ref['lognorm_sum']:.6g} vs {tst['lognorm_sum']:.6g} (rel {e:.2e})")

# ---- identity and ordering of names ----
check("first cell ids", ref["first_cells"] == tst["first_cells"],
      f"{ref['first_cells']} vs {tst['first_cells']}")
check("first gene ids", ref["first_genes"] == tst["first_genes"],
      f"{ref['first_genes']} vs {tst['first_genes']}")

# ---- cell type composition ----
rc, tc = ref["cell_type_counts"], tst["cell_type_counts"]
check("cell type labels", set(rc) == set(tc),
      f"ref-only={sorted(set(rc) - set(tc))} test-only={sorted(set(tc) - set(rc))}")
if set(rc) == set(tc):
    bad = {k: (rc[k], tc[k]) for k in rc if int(rc[k]) != int(tc[k])}
    check("cell type counts", not bad, str(bad) if bad else f"{len(rc)} types match")

# ---- per-cell QC, recomputed independently on each side ----
pc_r = pd.read_csv(os.path.join(args.ref, "per_cell.csv"))
pc_t = pd.read_csv(os.path.join(args.test, "per_cell.csv"))
check("per_cell rows", len(pc_r) == len(pc_t), f"{len(pc_r)} vs {len(pc_t)}")
if len(pc_r) == len(pc_t):
    check("per_cell order", (pc_r.cell_id.astype(str) == pc_t.cell_id.astype(str)).all(),
          "cell_id order identical")
    for col in ("total_counts", "n_genes"):
        d = np.abs(pc_r[col].to_numpy(float) - pc_t[col].to_numpy(float))
        check(f"per_cell {col}", d.max() <= args.tol, f"max abs diff {d.max():.3e}")
    check("per_cell cell_type",
          (pc_r.cell_type.astype(str) == pc_t.cell_type.astype(str)).all(),
          "labels aligned per cell")

# ---- per-gene totals ----
pg_r = pd.read_csv(os.path.join(args.ref, "per_gene.csv"))
pg_t = pd.read_csv(os.path.join(args.test, "per_gene.csv"))
check("per_gene rows", len(pg_r) == len(pg_t), f"{len(pg_r)} vs {len(pg_t)}")
if len(pg_r) == len(pg_t):
    check("per_gene order", (pg_r.gene.astype(str) == pg_t.gene.astype(str)).all(),
          "gene order identical")
    for col in ("total_counts", "n_cells"):
        d = np.abs(pg_r[col].to_numpy(float) - pg_t[col].to_numpy(float))
        check(f"per_gene {col}", d.max() <= args.tol, f"max abs diff {d.max():.3e}")

# ---- embeddings ----
em_r = pd.read_csv(os.path.join(args.ref, "embeddings.csv"))
em_t = pd.read_csv(os.path.join(args.test, "embeddings.csv"))
keys_r, keys_t = set(em_r.embedding.unique()), set(em_t.embedding.unique())
check("embedding keys", keys_r == keys_t,
      f"ref={sorted(keys_r)} test={sorted(keys_t)}")
for k in sorted(keys_r & keys_t):
    sr = em_r[em_r.embedding == k].reset_index(drop=True)
    st = em_t[em_t.embedding == k].reset_index(drop=True)
    if len(sr) != len(st):
        check(f"embedding {k}", False, f"{len(sr)} vs {len(st)} rows")
        continue
    same_order = (sr.cell_id.astype(str) == st.cell_id.astype(str)).all()
    d1 = np.abs(sr.dim1.to_numpy(float) - st.dim1.to_numpy(float))
    d2 = np.abs(np.nan_to_num(sr.dim2.to_numpy(float)) -
                np.nan_to_num(st.dim2.to_numpy(float)))
    check(f"embedding {k}", same_order and d1.max() <= args.tol and d2.max() <= args.tol,
          f"max abs diff dim1={d1.max():.3e} dim2={d2.max():.3e}")

# ---- dot plot ----
# Mean over a cell-type block: R and numpy sum in different orders, so this is
# the one block compared at float32-ish precision rather than exactly.
dot_tol = max(args.tol, 1e-9)
dp_r = pd.read_csv(os.path.join(args.ref, "dotplot.csv"))
dp_t = pd.read_csv(os.path.join(args.test, "dotplot.csv"))
m = dp_r.merge(dp_t, on=["cell_type", "gene"], suffixes=("_r", "_t"))
check("dotplot coverage", len(m) == len(dp_r) == len(dp_t),
      f"ref={len(dp_r)} test={len(dp_t)} joined={len(m)}")
if len(m):
    de = np.abs(m.mean_expr_r.to_numpy(float) - m.mean_expr_t.to_numpy(float))
    dp = np.abs(m.pct_expressed_r.to_numpy(float) - m.pct_expressed_t.to_numpy(float))
    check("dotplot mean_expr", de.max() <= dot_tol, f"max abs diff {de.max():.3e}")
    check("dotplot pct_expressed", dp.max() <= dot_tol, f"max abs diff {dp.max():.3e}")

# ---- report ----
width = max(len(n) for n, _, _ in results) + 2
title = f"ON-RAMP FIDELITY: {args.label}" if args.label else "ON-RAMP FIDELITY"
print("=" * 78)
print(title)
print("=" * 78)
for name, ok, detail in results:
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}} {detail}")
n_fail = sum(1 for _, ok, _ in results if not ok)
print("-" * 78)
print(f"  {len(results) - n_fail}/{len(results)} checks passed")
print("=" * 78)
sys.exit(1 if n_fail else 0)
