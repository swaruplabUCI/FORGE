#!/usr/bin/env python3
"""
nmf_enhancer_programs.py — Shi Fig 2B equivalent.

Decompose pseudobulk distal-enhancer accessibility into NMF programs and
render a heatmap of program activity (rows) × pseudobulk samples (cols).

Inputs:
  --peak-matrix   peak_matrix_annotated.h5ad with cell_type_broad + condition + sample obs
  --enhancers     ccan_enhancer_peaks.bed[.gz] (used to subset to distal CCAN-linked peaks)
  --k-grid        comma-separated NMF component counts to test (default: 6,8,10,12)
  --pick-k        if set, force this k; otherwise pick by elbow on reconstruction error
  --min-cells     skip pseudobulk groups with fewer than this many cells

Outputs (in --outdir):
  nmf_program_heatmap.{pdf,png}
  nmf_elbow.{pdf,png}
  nmf_W.tsv  (programs × peaks loadings)
  nmf_H.tsv  (programs × pseudobulk samples activities)
  nmf_summary.json
"""

import argparse
import gzip
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peak-matrix", required=True)
    p.add_argument("--enhancers", required=True)
    p.add_argument("--cell-type-key", default="cell_type_broad")
    p.add_argument("--condition-key", default="condition")
    p.add_argument("--sample-key", default="sample")
    p.add_argument("--k-grid", default="6,8,10,12")
    p.add_argument("--pick-k", type=int, default=0)
    p.add_argument("--min-cells", type=int, default=200)
    p.add_argument("--top-celltypes", type=int, default=12,
                   help="Restrict to the N largest cell-type-broad groups for legibility")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def load_enhancer_peak_ids(bed_path):
    rows = []
    with _open(bed_path) as f:
        first = f.readline().strip()
        # Detect header
        try:
            int(first.split("\t")[1])
            f.seek(0)
        except (IndexError, ValueError):
            pass
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[1].isdigit():
                continue
            rows.append(f"{parts[0]}:{parts[1]}-{parts[2]}")
    return set(rows)


def pseudobulk_aggregate(adata, peak_mask, cell_type_key, condition_key, sample_key,
                        min_cells, top_celltypes):
    obs = adata.obs[[cell_type_key, condition_key, sample_key]].copy()
    obs.columns = ["ct", "cond", "sample"]
    if top_celltypes:
        keep_cts = obs["ct"].value_counts().head(top_celltypes).index.tolist()
    else:
        keep_cts = sorted(obs["ct"].unique())

    groups = []
    samples = []
    for (ct, cond, samp), idx in obs.groupby(["ct", "cond", "sample"]).groups.items():
        if ct not in keep_cts:
            continue
        if len(idx) < min_cells:
            continue
        groups.append(idx)
        samples.append((ct, cond, samp))

    print(f"[2B] {len(samples)} pseudobulk samples (>={min_cells} cells, "
          f"top {len(keep_cts)} cell types)", flush=True)
    if len(samples) < 4:
        raise SystemExit("[2B] too few pseudobulk samples for NMF")

    peak_idx = np.where(peak_mask)[0]
    print(f"[2B] subsetting to {len(peak_idx)} enhancer peaks", flush=True)

    # Materialize only the enhancer columns once, in CSR. Backed-mode sparse
    # slicing returns one chunk at a time; concatenate per pseudobulk row set.
    obs_rank = {bc: i for i, bc in enumerate(adata.obs_names)}

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    enh_csr = X.tocsr()[:, peak_idx]
    print(f"[2B]   enhancer matrix shape: {enh_csr.shape}", flush=True)

    cols = []
    for idx in groups:
        rows = np.array(sorted(obs_rank[bc] for bc in idx), dtype=np.int64)
        sub = enh_csr[rows]
        col = np.asarray(sub.sum(axis=0)).ravel().astype(np.float32)
        col = col / max(len(rows), 1)
        cols.append(col)
    M = np.vstack(cols).T  # peaks × samples
    sample_labels = [f"{ct}::{cond}::{samp}" for (ct, cond, samp) in samples]
    return M, sample_labels, samples


def pick_k(M, k_grid):
    errs = []
    for k in k_grid:
        model = NMF(n_components=k, init="nndsvd", random_state=0, max_iter=400)
        W = model.fit_transform(M)
        H = model.components_
        recon = W @ H
        err = float(np.linalg.norm(M - recon, "fro") / (np.linalg.norm(M, "fro") + 1e-12))
        errs.append((k, err, W, H))
        print(f"[2B] k={k}  rel_err={err:.4f}", flush=True)
    return errs


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[2B] loading enhancer peak ids from {args.enhancers}", flush=True)
    enh_ids = load_enhancer_peak_ids(args.enhancers)
    print(f"[2B]   {len(enh_ids)} enhancer peaks", flush=True)

    print(f"[2B] reading peak matrix (in-memory) {args.peak_matrix}", flush=True)
    adata = ad.read_h5ad(args.peak_matrix)
    var_ids = pd.Index(adata.var_names.astype(str))
    peak_mask = var_ids.isin(enh_ids)
    n_match = int(peak_mask.sum())
    if n_match == 0:
        raise SystemExit("[2B] no enhancer peaks found in peak matrix var_names")
    print(f"[2B]   matched {n_match}/{len(enh_ids)} enhancer peaks", flush=True)

    M, sample_labels, samples = pseudobulk_aggregate(
        adata, peak_mask, args.cell_type_key, args.condition_key,
        args.sample_key, args.min_cells, args.top_celltypes)
    print(f"[2B] pseudobulk matrix: {M.shape[0]} peaks × {M.shape[1]} samples",
          flush=True)

    # Median-normalize columns to mute library-size effects
    col_sums = M.sum(axis=0)
    col_sums[col_sums == 0] = 1
    M_norm = M / col_sums * np.median(col_sums)

    k_grid = [int(k) for k in args.k_grid.split(",") if k.strip()]
    runs = pick_k(M_norm, k_grid)

    fig_e, ax_e = plt.subplots(figsize=(5, 3.5))
    ks = [r[0] for r in runs]
    errs = [r[1] for r in runs]
    ax_e.plot(ks, errs, "o-", color="#4c72b0")
    ax_e.set_xlabel("k (NMF components)")
    ax_e.set_ylabel("relative reconstruction error")
    ax_e.set_title("NMF elbow")
    fig_e.tight_layout()
    for ext in ("pdf", "png"):
        fig_e.savefig(outdir / f"nmf_elbow.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig_e)

    if args.pick_k:
        chosen = next((r for r in runs if r[0] == args.pick_k), runs[0])
    else:
        # Pick the k with the largest second derivative drop in error.
        if len(errs) >= 3:
            d2 = [errs[i - 1] - 2 * errs[i] + errs[i + 1] for i in range(1, len(errs) - 1)]
            idx = int(np.argmax(d2)) + 1
            chosen = runs[idx]
        else:
            chosen = runs[-1]
    k_chosen, err_chosen, W, H = chosen
    print(f"[2B] picked k={k_chosen} (relative err={err_chosen:.4f})", flush=True)

    var_subset = var_ids[peak_mask]
    pd.DataFrame(W, index=var_subset,
                 columns=[f"P{i+1}" for i in range(k_chosen)]).to_csv(
        outdir / "nmf_W.tsv", sep="\t")
    H_df = pd.DataFrame(H, index=[f"P{i+1}" for i in range(k_chosen)],
                       columns=sample_labels)
    H_df.to_csv(outdir / "nmf_H.tsv", sep="\t")

    # Heatmap rendering (rows = programs, cols = pseudobulk samples) ordered by ct → cond → sample
    order = sorted(range(len(samples)),
                  key=lambda i: (samples[i][0], samples[i][1], samples[i][2]))
    H_z = (H - H.mean(axis=1, keepdims=True)) / (H.std(axis=1, keepdims=True) + 1e-9)
    H_z_ord = H_z[:, order]
    labels_ord = [sample_labels[i] for i in order]

    fig, (ax_top, ax_main) = plt.subplots(
        2, 1, figsize=(max(8, 0.45 * len(labels_ord)), 0.45 * k_chosen + 3),
        gridspec_kw={"height_ratios": [1, k_chosen + 1]}, sharex=True)
    ct_codes = pd.Categorical([samples[i][0] for i in order]).codes
    cond_codes = pd.Categorical([samples[i][1] for i in order]).codes
    ax_top.imshow(np.vstack([ct_codes, cond_codes]),
                 aspect="auto", cmap="tab20")
    ax_top.set_yticks([0, 1]); ax_top.set_yticklabels(["cell type", "condition"])
    ax_top.set_xticks([])

    im = ax_main.imshow(H_z_ord, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax_main.set_yticks(range(k_chosen))
    ax_main.set_yticklabels([f"P{i+1}" for i in range(k_chosen)])
    ax_main.set_xticks(range(len(labels_ord)))
    ax_main.set_xticklabels(labels_ord, rotation=90, fontsize=6)
    plt.colorbar(im, ax=ax_main, fraction=0.02, pad=0.02, label="z-score")
    fig.suptitle(f"NMF programs over distal enhancers (k={k_chosen})", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"nmf_program_heatmap.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "k_chosen": int(k_chosen),
        "k_grid": k_grid,
        "errors": dict(zip([str(k) for k in ks], errs)),
        "n_peaks": int(M.shape[0]),
        "n_samples": int(M.shape[1]),
        "samples": [{"cell_type": ct, "condition": cond, "sample": s}
                    for (ct, cond, s) in samples],
    }
    (outdir / "nmf_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[2B] wrote heatmap, elbow, W/H matrices for k={k_chosen}", flush=True)


if __name__ == "__main__":
    main()
