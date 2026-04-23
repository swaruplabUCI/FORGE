#!/usr/bin/env python3
"""
make_cicero_triplets_stratified.py

Subset a cell x peak matrix by condition before generating Cicero triplets.
Mirrors Shi et al.'s approach of running Cicero separately per disease condition
to capture condition-specific co-accessibility patterns.

Inputs:
    - Peak matrix h5ad (with obs['condition'] column)
    - Condition key (column name in obs)
    - Condition value (e.g., '5xFAD' or 'SREBF1_OE')

Outputs:
    - cicero_triplets_{condition}.tsv.gz (cell, peak, count)
"""

import argparse
import gzip

import anndata as ad
import numpy as np
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser(
        description="Make condition-stratified Cicero triplets from .h5ad")
    p.add_argument("--peak-matrix", required=True,
                   help="Path to cell x peak matrix .h5ad")
    p.add_argument("--condition-key", required=True,
                   help="Column name in adata.obs for condition (e.g., 'condition')")
    p.add_argument("--condition-value", required=True,
                   help="Condition value to subset (e.g., '5xFAD')")
    p.add_argument("--out", required=True,
                   help="Output triplets file (e.g., cicero_triplets_5xFAD.tsv.gz)")
    return p.parse_args()


def main():
    args = parse_args()

    print(f"Loading peak matrix from {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)

    print(f"  Total cells: {adata.n_obs}")
    print(f"  Total peaks: {adata.n_vars}")

    # Validate condition column exists
    if args.condition_key not in adata.obs.columns:
        raise ValueError(
            f"Column '{args.condition_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}")

    conditions = adata.obs[args.condition_key].astype(str).unique()
    print(f"  Available conditions: {list(conditions)}")

    # Subset to requested condition
    mask = adata.obs[args.condition_key].astype(str) == args.condition_value
    n_cells = int(mask.sum())
    print(f"  Cells matching '{args.condition_value}': {n_cells}")

    if n_cells == 0:
        raise ValueError(
            f"No cells found for {args.condition_key}='{args.condition_value}'. "
            f"Available values: {list(conditions)}")

    adata_sub = adata[mask]

    X = adata_sub.X
    if not sp.isspmatrix(X):
        X = sp.csr_matrix(X)

    n_obs, n_var = adata_sub.n_obs, adata_sub.n_vars

    # Detect orientation (same logic as make_cicero_triplets.py)
    if X.shape == (n_obs, n_var):
        row_is_cell = True
        cells = np.asarray(adata_sub.obs_names, dtype=str)
        peaks = np.asarray(adata_sub.var_names, dtype=str)
    elif X.shape == (n_var, n_obs):
        row_is_cell = False
        cells = np.asarray(adata_sub.obs_names, dtype=str)
        peaks = np.asarray(adata_sub.var_names, dtype=str)
    else:
        raise ValueError(
            f"Unexpected X shape {X.shape} vs n_obs={n_obs}, n_vars={n_var}")

    print(f"  Matrix shape: {X.shape}")
    if row_is_cell:
        print(f"  Rows=cells ({len(cells)}), Cols=peaks ({len(peaks)})")
    else:
        print(f"  Rows=peaks ({len(peaks)}), Cols=cells ({len(cells)})")

    X = X.tocoo()
    print(f"  Non-zero entries: {X.nnz}")

    print(f"Writing triplets to {args.out}")
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

    print(f"Done. Wrote triplets for {n_cells} cells ({args.condition_value}).")


if __name__ == "__main__":
    main()
