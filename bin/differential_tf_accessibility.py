#!/usr/bin/env python3
"""
differential_tf_accessibility.py

Per-celltype TF motif accessibility testing on the ChromVAR deviation matrix.
Two modes:

  --mode differential (default, Shi et al. style):
      Wilcoxon of TF deviations between two conditions WITHIN a cell type.
      Requires --treatment and --control labels.
      Output: tf_differential_<cell_type>_<trt>_vs_<ctrl>.csv

  --mode descriptive:
      Wilcoxon of TF deviations for the target cell type vs the REST of the
      cells (ignoring condition). Produces a per-celltype TF ranking — useful
      as a standalone product on single-condition datasets. Optionally
      restrict to one condition via --treatment (treated here as a subset
      filter, not a contrast).
      Output: tf_descriptive_<cell_type>[_<condition>].csv

Inputs (both modes):
  - ChromVAR deviations h5ad (cells x TF motifs; X = z-score deviations)
  - Annotated peak matrix h5ad (same barcodes) providing cell_type + condition obs
"""

import argparse
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--chromvar-dev", required=True,
                   help="ChromVAR deviations h5ad (cells x TF motifs)")
    p.add_argument("--annotated-peaks", required=True,
                   help="peak_matrix_annotated.h5ad (supplies cell_type + condition obs)")
    p.add_argument("--cell-type", required=True,
                   help="Cell type label to test (must match obs['cell_type'] or cell_type_prediction)")
    p.add_argument("--mode", choices=["differential", "descriptive"],
                   default="differential",
                   help="'differential' = trt vs ctrl (default); "
                        "'descriptive' = cell type vs rest (single-condition friendly)")
    p.add_argument("--condition-key", default="condition",
                   help="obs key holding the condition label (default: condition)")
    p.add_argument("--treatment", default=None,
                   help="Differential: treatment label. Descriptive: optional "
                        "subset filter (restrict to cells in this condition).")
    p.add_argument("--control", default=None,
                   help="Differential: control label. Ignored in descriptive mode.")
    p.add_argument("--min-cells", type=int, default=50,
                   help="Minimum cells per group required to run the test")
    p.add_argument("--fdr-cutoff", type=float, default=0.05,
                   help="FDR cutoff used for the summary count (does not filter CSV)")
    p.add_argument("--out-prefix", default=None,
                   help="Output filename prefix; default depends on mode")
    return p.parse_args()


def sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s))


def _resolve_cell_type_column(obs: pd.DataFrame) -> str:
    for col in ("cell_type", "cell_type_prediction", "celltypist_prediction", "pred_y"):
        if col in obs.columns:
            return col
    raise KeyError(
        "No cell-type column found on annotated peak matrix (looked for: "
        "cell_type, cell_type_prediction, celltypist_prediction, pred_y)"
    )


def _run_differential(dev, args):
    ct_tag   = sanitize(args.cell_type)
    trt_tag  = sanitize(args.treatment)
    ctrl_tag = sanitize(args.control)
    prefix   = args.out_prefix or "tf_differential"
    out_csv     = f"{prefix}_{ct_tag}_{trt_tag}_vs_{ctrl_tag}.csv"
    out_summary = f"{prefix}_{ct_tag}_{trt_tag}_vs_{ctrl_tag}.summary"

    mask = (dev.obs["cell_type"] == args.cell_type) & \
           (dev.obs["condition"].isin([args.treatment, args.control]))
    sub = dev[mask].copy()

    n_trt  = int((sub.obs["condition"] == args.treatment).sum())
    n_ctrl = int((sub.obs["condition"] == args.control).sum())
    print(f"[{args.cell_type}] {args.treatment}={n_trt}  {args.control}={n_ctrl}  "
          f"(min_cells={args.min_cells})", flush=True)

    if n_trt < args.min_cells or n_ctrl < args.min_cells:
        msg = (f"SKIPPED: cell_type={args.cell_type} trt({args.treatment})={n_trt} "
               f"ctrl({args.control})={n_ctrl} below min_cells={args.min_cells}")
        Path(out_summary).write_text(msg + "\n")
        print(msg, flush=True)
        return

    sub.obs["condition"] = sub.obs["condition"].astype("category")
    sc.tl.rank_genes_groups(
        sub, groupby="condition", groups=[args.treatment], reference=args.control,
        method="wilcoxon", use_raw=False, pts=True,
    )

    rg = sub.uns["rank_genes_groups"]
    df = pd.DataFrame({
        "motif":          pd.Series(rg["names"][args.treatment]).astype(str),
        "score":          pd.Series(rg["scores"][args.treatment]).astype(float),
        "logfoldchange":  pd.Series(rg["logfoldchanges"][args.treatment]).astype(float),
        "pval":           pd.Series(rg["pvals"][args.treatment]).astype(float),
        "pval_adj":       pd.Series(rg["pvals_adj"][args.treatment]).astype(float),
    })
    if "pts" in rg:
        df["pct_trt"]  = pd.Series(rg["pts"][args.treatment]).astype(float).values
        if args.control in rg["pts"].dtype.names if hasattr(rg["pts"], "dtype") else rg["pts"]:
            df["pct_ctrl"] = pd.Series(rg["pts"][args.control]).astype(float).values
    df["cell_type"]  = args.cell_type
    df["mode"]       = "differential"
    df["treatment"]  = args.treatment
    df["control"]    = args.control
    df = df.sort_values("pval_adj", kind="mergesort").reset_index(drop=True)
    df.to_csv(out_csv, index=False)

    n_sig = int((df["pval_adj"] < args.fdr_cutoff).sum())
    summary = (f"mode=differential cell_type={args.cell_type} trt={args.treatment} "
               f"ctrl={args.control} n_trt={n_trt} n_ctrl={n_ctrl} "
               f"n_motifs={len(df)} n_sig(FDR<{args.fdr_cutoff})={n_sig}")
    Path(out_summary).write_text(summary + "\n")
    print(summary, flush=True)


def _run_descriptive(dev, args):
    """Rank TFs for cell_type vs rest-of-cells. Single-condition friendly.

    Optional --treatment acts as a condition filter (subset before ranking),
    not a contrast. If None, uses all cells regardless of condition.
    """
    ct_tag = sanitize(args.cell_type)
    prefix = args.out_prefix or "tf_descriptive"

    if args.treatment:
        cond_tag = sanitize(args.treatment)
        out_csv     = f"{prefix}_{ct_tag}_{cond_tag}.csv"
        out_summary = f"{prefix}_{ct_tag}_{cond_tag}.summary"
        sub = dev[dev.obs["condition"] == args.treatment].copy()
    else:
        cond_tag = "all"
        out_csv     = f"{prefix}_{ct_tag}.csv"
        out_summary = f"{prefix}_{ct_tag}.summary"
        sub = dev.copy()

    n_ct   = int((sub.obs["cell_type"] == args.cell_type).sum())
    n_rest = int((sub.obs["cell_type"] != args.cell_type).sum())
    print(f"[{args.cell_type}] in={n_ct}  rest={n_rest}  "
          f"cond_filter={args.treatment or 'none'}  (min_cells={args.min_cells})",
          flush=True)

    if n_ct < args.min_cells or n_rest < args.min_cells:
        msg = (f"SKIPPED: cell_type={args.cell_type} n={n_ct} rest={n_rest} "
               f"below min_cells={args.min_cells}")
        Path(out_summary).write_text(msg + "\n")
        print(msg, flush=True)
        return

    sub.obs["_group"] = np.where(
        sub.obs["cell_type"] == args.cell_type, args.cell_type, "rest"
    )
    sub.obs["_group"] = sub.obs["_group"].astype("category")
    sc.tl.rank_genes_groups(
        sub, groupby="_group", groups=[args.cell_type], reference="rest",
        method="wilcoxon", use_raw=False, pts=True,
    )

    rg = sub.uns["rank_genes_groups"]
    df = pd.DataFrame({
        "motif":          pd.Series(rg["names"][args.cell_type]).astype(str),
        "score":          pd.Series(rg["scores"][args.cell_type]).astype(float),
        "logfoldchange":  pd.Series(rg["logfoldchanges"][args.cell_type]).astype(float),
        "pval":           pd.Series(rg["pvals"][args.cell_type]).astype(float),
        "pval_adj":       pd.Series(rg["pvals_adj"][args.cell_type]).astype(float),
    })
    if "pts" in rg:
        df["pct_in"]   = pd.Series(rg["pts"][args.cell_type]).astype(float).values
        df["pct_rest"] = pd.Series(rg["pts"]["rest"]).astype(float).values
    df["cell_type"] = args.cell_type
    df["mode"]      = "descriptive"
    df["condition"] = args.treatment or "all"
    df = df.sort_values("pval_adj", kind="mergesort").reset_index(drop=True)
    df.to_csv(out_csv, index=False)

    n_sig = int((df["pval_adj"] < args.fdr_cutoff).sum())
    summary = (f"mode=descriptive cell_type={args.cell_type} condition={cond_tag} "
               f"n_in={n_ct} n_rest={n_rest} n_motifs={len(df)} "
               f"n_sig(FDR<{args.fdr_cutoff})={n_sig}")
    Path(out_summary).write_text(summary + "\n")
    print(summary, flush=True)


def main():
    args = parse_args()

    if args.mode == "differential":
        if not args.treatment or not args.control:
            raise ValueError("--mode differential requires both --treatment and --control")

    # Load ChromVAR deviations and annotated peak matrix
    dev = ad.read_h5ad(args.chromvar_dev)
    ann = ad.read_h5ad(args.annotated_peaks)

    ct_col = _resolve_cell_type_column(ann.obs)
    need_condition = (args.mode == "differential") or (args.treatment is not None)
    if need_condition and args.condition_key not in ann.obs.columns:
        raise KeyError(f"condition_key='{args.condition_key}' not in annotated_peaks obs "
                       f"(available: {list(ann.obs.columns)})")

    # Intersect barcodes and join labels onto dev
    common = dev.obs_names.intersection(ann.obs_names)
    if len(common) == 0:
        raise RuntimeError("No shared barcodes between chromvar_dev and annotated_peaks")
    dev = dev[common].copy()
    cols = [ct_col]
    if need_condition:
        cols.append(args.condition_key)
    labels = ann.obs.loc[common, cols].copy()
    labels.columns = ["cell_type"] + (["condition"] if need_condition else [])
    dev.obs = dev.obs.join(labels)

    if args.mode == "differential":
        _run_differential(dev, args)
    else:
        _run_descriptive(dev, args)


if __name__ == "__main__":
    main()
