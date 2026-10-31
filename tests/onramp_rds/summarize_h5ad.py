#!/usr/bin/env python3
"""Converted side of the .rds on-ramp fidelity check.

Reads the h5ad produced by bin/convert_seurat_to_h5ad.R and emits artifacts
with the same filenames and columns as summarize_seurat.R. compare.py diffs
the pair.

Runs in scgpu_extended.sif -- a different container from the R side on
purpose. FORGE hands this h5ad across exactly that boundary, so if a string
came back as bytes or a matrix came back transposed, it surfaces here rather
than three stages downstream.

Every quantity is recomputed from X, never read from a cached obs column, so a
converter that copied a stale column cannot pass by echoing it back.
"""

import argparse
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser()
p.add_argument("--h5ad", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--cell_type_key", default="cell_type")
p.add_argument("--markers", default="")
args = p.parse_args()

os.makedirs(args.outdir, exist_ok=True)

a = ad.read_h5ad(args.h5ad)
X = a.X if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
X = X.tocsr()

cells = a.obs_names.astype(str).tolist()
genes = a.var_names.astype(str).tolist()

ct_key = args.cell_type_key
if ct_key not in a.obs.columns:
    for cand in ("cell_type", "scanvi_prediction", "celltypist_prediction",
                 "cell_type_prediction", "cell_type_marker", "seurat_annotations"):
        if cand in a.obs.columns:
            ct_key = cand
            break
cell_type = a.obs[ct_key].astype(str).to_numpy()

# ---- per-cell QC, recomputed from X ----
total_counts = np.asarray(X.sum(axis=1)).ravel()
n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
pd.DataFrame({
    "cell_id": cells,
    "total_counts": total_counts,
    "n_genes": n_genes,
    "cell_type": cell_type,
}).to_csv(os.path.join(args.outdir, "per_cell.csv"), index=False)

# ---- per-gene totals ----
Xc = X.tocsc()
pd.DataFrame({
    "gene": genes,
    "total_counts": np.asarray(Xc.sum(axis=0)).ravel(),
    "n_cells": np.asarray((Xc > 0).sum(axis=0)).ravel(),
}).to_csv(os.path.join(args.outdir, "per_gene.csv"), index=False)

# ---- embeddings ----
rows = []
for key, mat in a.obsm.items():
    m = np.asarray(mat)
    rows.append(pd.DataFrame({
        "cell_id": cells,
        "embedding": key,
        "dim1": m[:, 0],
        "dim2": m[:, 1] if m.shape[1] > 1 else np.nan,
    }))
emb_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
    columns=["cell_id", "embedding", "dim1", "dim2"])
emb_df.to_csv(os.path.join(args.outdir, "embeddings.csv"), index=False)

# ---- dot-plot quantities ----
if "logcounts" in a.layers:
    E = a.layers["logcounts"]
    E = E.tocsr() if sparse.issparse(E) else sparse.csr_matrix(E)
else:
    E = X

gene_pos = {g: i for i, g in enumerate(genes)}
if args.markers:
    marker_genes = [g.strip() for g in args.markers.split(",") if g.strip() in gene_pos]
else:
    tot = np.asarray(E.sum(axis=0)).ravel()
    # Match R's order(decreasing=TRUE), which breaks ties by first index.
    order = sorted(range(len(genes)), key=lambda i: (-tot[i], i))
    marker_genes = [genes[i] for i in order[:min(20, len(genes))]]

midx = [gene_pos[g] for g in marker_genes]
Esub = E[:, midx].tocsc()

recs = []
for ct in sorted(set(cell_type[(cell_type != "") & (cell_type != "nan")])):
    sel = np.flatnonzero(cell_type == ct)
    block = Esub[sel, :]
    mean_expr = np.asarray(block.mean(axis=0)).ravel()
    pct = np.asarray((block > 0).sum(axis=0)).ravel() / len(sel) * 100.0
    for j, g in enumerate(marker_genes):
        recs.append({"cell_type": ct, "gene": g,
                     "mean_expr": float(mean_expr[j]),
                     "pct_expressed": float(pct[j])})
pd.DataFrame(recs).to_csv(os.path.join(args.outdir, "dotplot.csv"), index=False)

# ---- scalar summary ----
num_cols = [c for c in a.obs.columns if pd.api.types.is_numeric_dtype(a.obs[c])]
summary = {
    "source": "converted_h5ad",
    "n_cells": int(a.n_obs),
    "n_genes": int(a.n_vars),
    "x_sum": float(X.sum()),
    "x_nnz": int(X.nnz),
    "has_lognorm": "logcounts" in a.layers,
    "lognorm_sum": float(a.layers["logcounts"].sum()) if "logcounts" in a.layers else None,
    "cell_type_key": ct_key,
    "cell_type_counts": {k: int(v) for k, v in pd.Series(cell_type).value_counts().items()},
    "embeddings": {k: int(np.asarray(v).shape[0]) for k, v in a.obsm.items()},
    "obs_numeric_means": {c: float(np.nanmean(a.obs[c].to_numpy(dtype=float))) for c in num_cols},
    "first_cells": cells[:5],
    "first_genes": genes[:5],
    "marker_genes": marker_genes,
}
with open(os.path.join(args.outdir, "summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

print(f"Wrote h5ad-side artifacts to {args.outdir}")
print(f"  cells={a.n_obs} genes={a.n_vars} x_sum={X.sum():.0f} cell_type_key={ct_key}")
