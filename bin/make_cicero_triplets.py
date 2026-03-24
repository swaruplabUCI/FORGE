#!/usr/bin/env python3
"""
make_cicero_triplets.py

Convert a cell x peak matrix in .h5ad format into a Cicero-style triplets file:
    cell   peak   count

The script is robust to orientation:
- If rows = cells, cols = peaks, it uses that directly.
- If rows = peaks, cols = cells, it swaps axes accordingly.
"""

import argparse
import gzip

import anndata as ad
import numpy as np
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser(description="Make Cicero triplets from .h5ad cellxpeak matrix")
    p.add_argument(
        "--peak-matrix",
        required=True,
        help="Path to cell x peak matrix .h5ad (e.g. cellxpeak_matrix.h5ad)",
    )
    p.add_argument(
        "--out",
        required=True,
        help="Output triplets file (e.g. cicero_triplets.tsv.gz)",
    )
    return p.parse_args()


def main():
    args = parse_args()

    print(f"Loading peak matrix from {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)

    X = adata.X
    if not sp.isspmatrix(X):
        X = sp.csr_matrix(X)

    n_obs, n_var = adata.n_obs, adata.n_vars

    # Detect orientation
    if X.shape == (n_obs, n_var):
        # rows = obs, cols = var
        # interpret obs as cells, var as peaks
        row_is_cell = True
        cells = np.asarray(adata.obs_names, dtype=str)
        peaks = np.asarray(adata.var_names, dtype=str)
    elif X.shape == (n_var, n_obs):
        # rows = var, cols = obs
        # interpret var as peaks, obs as cells
        row_is_cell = False
        cells = np.asarray(adata.obs_names, dtype=str)
        peaks = np.asarray(adata.var_names, dtype=str)
    else:
        raise ValueError(
            f"Unexpected X shape {X.shape} vs n_obs={n_obs}, n_vars={n_var}. "
            "Cannot infer cells/peaks orientation."
        )

    print("Matrix shape:", X.shape)
    if row_is_cell:
        print(f"Interpreting rows as cells (n={len(cells)}) and columns as peaks (n={len(peaks)})")
    else:
        print(f"Interpreting rows as peaks (n={len(peaks)}) and columns as cells (n={len(cells)})")

    X = X.tocoo()
    print(f"Non-zero entries: {X.nnz}")

    out_path = args.out
    print(f"Writing triplets to {out_path}")
    with gzip.open(out_path, "wt") as f:
        f.write("cell\tpeak\tcount\n")
        if row_is_cell:
            # X[row, col] => cells[row], peaks[col]
            for i, j, v in zip(X.row, X.col, X.data):
                if v != 0:
                    f.write(f"{cells[i]}\t{peaks[j]}\t{int(v)}\n")
        else:
            # X[row, col] => peaks[row], cells[col]
            for i, j, v in zip(X.row, X.col, X.data):
                if v != 0:
                    f.write(f"{cells[j]}\t{peaks[i]}\t{int(v)}\n")

    print("Done.")


if __name__ == "__main__":
    main()
