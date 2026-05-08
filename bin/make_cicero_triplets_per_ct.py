#!/usr/bin/env python3
"""
make_cicero_triplets_per_ct.py

Subset a cell x peak matrix by both cell type and condition before generating
Cicero triplets. Extends the condition-only stratification of
make_cicero_triplets_stratified.py to the full per-CT × condition design.

Cell type membership is read from ct_annotation_v2.csv (produced by
build_ct_annotation_v2.py) rather than from an obs column in the h5ad, because
the v2 label hierarchy (glia fine / neural broad / vascular) is not stored
directly as a single obs column.

Inputs:
    --peak-matrix     : cell x peak matrix .h5ad
    --annotation-csv  : ct_annotation_v2.csv  (barcode, condition, cell_type_v2)
    --cell-type       : cell_type_v2 value (e.g. "Oligo NN", "IT-ET Glut")
    --condition       : condition value      (e.g. "5xFAD", "SREBF1_OE")
    --out             : output triplets file (.tsv.gz)

Output:
    TSV with columns: cell  peak  count
"""

import argparse
import gzip
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peak-matrix",    required=True,
                   help="Path to annotated peak matrix .h5ad (peak_matrix_annotated.h5ad)")
    p.add_argument("--annotation-csv", required=True,
                   help="Path to ct_annotation_v2.csv (barcode, condition, cell_type_v2)")
    p.add_argument("--cell-type",      required=True,
                   help="cell_type_v2 value to subset (e.g. 'Oligo NN', 'IT-ET Glut')")
    p.add_argument("--condition",      required=True,
                   help="Condition value to subset (e.g. '5xFAD', 'SREBF1_OE')")
    p.add_argument("--out",            required=True,
                   help="Output triplets file (e.g. cicero_triplets_OligoNN_5xFAD.tsv.gz)")
    p.add_argument("--min-cells",      type=int, default=0,
                   help="Minimum cells per stratum for Cicero reliability. "
                        "Exit 77 if below threshold (Nextflow: errorStrategy ignore). "
                        "Default 0 = no floor.")
    return p.parse_args()


def main():
    args = parse_args()

    # ── load annotation CSV ───────────────────────────────────────────────────
    print(f"[load] annotation CSV: {args.annotation_csv}")
    ann = pd.read_csv(args.annotation_csv, dtype=str)

    required = {"barcode", "condition", "cell_type_v2"}
    missing = required - set(ann.columns)
    if missing:
        sys.exit(f"ERROR: annotation CSV missing columns: {missing}")

    mask = (ann["cell_type_v2"] == args.cell_type) & (ann["condition"] == args.condition)
    barcodes = set(ann.loc[mask, "barcode"])
    n_barcodes = len(barcodes)
    print(f"[info] {args.cell_type} × {args.condition}: {n_barcodes} barcodes in annotation")

    if n_barcodes == 0:
        sys.exit(f"ERROR: no barcodes found for "
                 f"cell_type_v2='{args.cell_type}' AND condition='{args.condition}'. "
                 f"Available cell types: {sorted(ann['cell_type_v2'].unique())}")

    # ── load h5ad and subset ──────────────────────────────────────────────────
    print(f"[load] peak matrix: {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)
    print(f"[info] full matrix: {adata.n_obs:,} cells × {adata.n_vars:,} peaks")

    obs_mask = pd.Series(adata.obs_names).isin(barcodes).values
    n_matched = int(obs_mask.sum())
    if n_matched == 0:
        sys.exit("ERROR: none of the annotation barcodes matched obs_names in the h5ad. "
                 "Check that both files use the same barcode format "
                 "(e.g. 'sample1_batch1:AAACGAAA...-1').")
    if n_matched < n_barcodes:
        print(f"[warn] {n_barcodes - n_matched} barcodes in annotation not found in h5ad "
              f"(expected 0 — investigate if large).")

    if args.min_cells > 0 and n_matched < args.min_cells:
        print(f"[skip] {args.cell_type} × {args.condition}: {n_matched} cells < "
              f"min_cells={args.min_cells}. Exiting with code 77 (Nextflow: ignore).")
        sys.exit(77)

    adata_sub = adata[obs_mask]
    print(f"[info] subset: {adata_sub.n_obs:,} cells × {adata_sub.n_vars:,} peaks")

    # ── build triplets ────────────────────────────────────────────────────────
    X = adata_sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    n_obs, n_var = adata_sub.n_obs, adata_sub.n_vars

    if X.shape == (n_obs, n_var):
        cells = np.asarray(adata_sub.obs_names, dtype=str)
        peaks = np.asarray(adata_sub.var_names, dtype=str)
        row_is_cell = True
    elif X.shape == (n_var, n_obs):
        cells = np.asarray(adata_sub.obs_names, dtype=str)
        peaks = np.asarray(adata_sub.var_names, dtype=str)
        row_is_cell = False
    else:
        sys.exit(f"ERROR: unexpected X shape {X.shape} vs n_obs={n_obs}, n_vars={n_var}")

    X = X.tocoo()
    print(f"[info] non-zero entries: {X.nnz:,}")

    # ── write ─────────────────────────────────────────────────────────────────
    print(f"[write] {args.out}")
    with gzip.open(args.out, "wt") as f:
        f.write("cell\tpeak\tcount\n")
        if row_is_cell:
            for i, j, v in zip(X.row, X.col, X.data):
                if v != 0:
                    f.write(f"{cells[i]}\t{peaks[j]}\t{int(v)}\n")
        else:
            for i, j, v in zip(X.row, X.col, X.data):
                if v != 0:
                    f.write(f"{cells[j]}\t{peaks[i]}\t{int(v)}\n")

    print(f"[done] {n_matched} cells written for "
          f"{args.cell_type} × {args.condition}.")


if __name__ == "__main__":
    main()
