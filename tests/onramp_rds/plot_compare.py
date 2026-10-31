#!/usr/bin/env python3
"""Render the Seurat-side and h5ad-side artifacts as side-by-side figures.

The numeric gate lives in compare.py; this is the eyeball check that sits on
top of it. Both panels in every figure are drawn from CSVs -- the left one
from the original .rds, the right one from the converted h5ad -- by the same
matplotlib code, so any visible difference is a difference in the data and
not in how two libraries chose to draw it.

Produces, under --outdir:
    umap_compare.png      embedding coloured by cell type
    qc_compare.png        counts/genes distributions and per-cell-type box plots
    dotplot_compare.png   mean expression + percent expressed per cell type
    delta.png             direct scatter of ref vs test for every shared metric
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--ref", required=True)
p.add_argument("--test", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--label", default="")
p.add_argument("--embedding", default="", help="which obsm key to plot; default picks umap then tsne then pca")
args = p.parse_args()

os.makedirs(args.outdir, exist_ok=True)
SIDES = [("original .rds", args.ref), ("converted h5ad", args.test)]


def palette(cats):
    cmap = plt.get_cmap("tab20")
    return {c: cmap(i % 20) for i, c in enumerate(sorted(cats))}


def suptitle(fig, what):
    fig.suptitle(f"{what} — {args.label}" if args.label else what,
                 fontsize=13, fontweight="bold")


# ---------------------------------------------------------------- UMAP ----
emb = {name: pd.read_csv(os.path.join(d, "embeddings.csv")) for name, d in SIDES}
pcs = {name: pd.read_csv(os.path.join(d, "per_cell.csv")) for name, d in SIDES}

available = set(emb["original .rds"].embedding.unique()) & set(emb["converted h5ad"].embedding.unique())
if args.embedding and args.embedding in available:
    key = args.embedding
else:
    key = next((k for k in ("X_umap", "X_tsne", "X_pca") if k in available),
               sorted(available)[0] if available else None)

if key:
    all_ct = set(pcs["original .rds"].cell_type.astype(str)) | set(pcs["converted h5ad"].cell_type.astype(str))
    colors = palette(all_ct)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    for ax, (name, _) in zip(axes, SIDES):
        e = emb[name][emb[name].embedding == key].reset_index(drop=True)
        pc = pcs[name].set_index("cell_id").reindex(e.cell_id.astype(str))
        ct = pc.cell_type.astype(str).to_numpy()
        for c in sorted(set(ct)):
            m = ct == c
            ax.scatter(e.dim1[m], e.dim2[m], s=4, alpha=0.7,
                       color=colors[c], label=c, linewidths=0)
        ax.set_title(f"{name}  (n={len(e)})")
        ax.set_xlabel(f"{key}_1"); ax.set_ylabel(f"{key}_2")
        ax.set_xticks([]); ax.set_yticks([])
    handles, labels = axes[0].get_legend_handles_labels()
    ncol = 1 if len(labels) <= 16 else 2
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=7, markerscale=2.5, frameon=False, ncol=ncol)
    suptitle(fig, f"{key} coloured by cell_type")
    fig.tight_layout(rect=[0, 0, 0.88, 0.95])
    fig.savefig(os.path.join(args.outdir, "umap_compare.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote umap_compare.png")

# ------------------------------------------------------------------ QC ----
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
for col, (name, _) in enumerate(SIDES):
    pc = pcs[name]
    ax = axes[0, col]
    ax.hist(np.log10(pc.total_counts.clip(lower=1)), bins=60, color="#4C72B0", alpha=0.85)
    ax.set_title(f"{name}: total counts/cell")
    ax.set_xlabel("log10(total counts)"); ax.set_ylabel("cells")
    ax.axvline(np.log10(max(pc.total_counts.median(), 1)), color="crimson", ls="--", lw=1,
               label=f"median {pc.total_counts.median():.0f}")
    ax.legend(fontsize=8)

    ax = axes[1, col]
    order = sorted(pc.cell_type.astype(str).unique())
    data = [pc.n_genes[pc.cell_type.astype(str) == c].to_numpy() for c in order]
    ax.boxplot(data, labels=order, showfliers=False)
    ax.set_title(f"{name}: genes detected per cell type")
    ax.set_ylabel("n_genes")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
suptitle(fig, "QC distributions")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(args.outdir, "qc_compare.png"), dpi=130)
plt.close(fig)
print("wrote qc_compare.png")

# ------------------------------------------------------------- dot plot ----
dps = {name: pd.read_csv(os.path.join(d, "dotplot.csv")) for name, d in SIDES}
genes = list(dict.fromkeys(dps["original .rds"].gene))
cts = sorted(dps["original .rds"].cell_type.unique())
vmax = max(dps["original .rds"].mean_expr.max(), dps["converted h5ad"].mean_expr.max(), 1e-12)

fig, axes = plt.subplots(1, 2, figsize=(max(12, len(genes) * 0.85), max(5, len(cts) * 0.40) + 2))
for ax, (name, _) in zip(axes, SIDES):
    d = dps[name]
    piv_m = d.pivot(index="cell_type", columns="gene", values="mean_expr").reindex(index=cts, columns=genes)
    piv_p = d.pivot(index="cell_type", columns="gene", values="pct_expressed").reindex(index=cts, columns=genes)
    gx, cy = np.meshgrid(np.arange(len(genes)), np.arange(len(cts)))
    sc = ax.scatter(gx.ravel(), cy.ravel(),
                    s=(piv_p.to_numpy().ravel() / 100.0) * 150 + 1,
                    c=piv_m.to_numpy().ravel(), cmap="viridis", vmin=0, vmax=vmax,
                    edgecolors="none")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts, fontsize=7)
    ax.set_title(name)
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(-0.6, len(cts) - 0.4)
fig.colorbar(sc, ax=axes, shrink=0.7, label="mean expression")
suptitle(fig, "Dot plot (size = % expressing, colour = mean expression)")
fig.savefig(os.path.join(args.outdir, "dotplot_compare.png"), dpi=130, bbox_inches="tight")
plt.close(fig)
print("wrote dotplot_compare.png")

# ----------------------------------------------------------- ref vs test ----
# If the conversion is faithful every one of these lands exactly on y = x.
pg = {name: pd.read_csv(os.path.join(d, "per_gene.csv")) for name, d in SIDES}
pairs = [
    ("per-cell total counts", pcs["original .rds"].total_counts, pcs["converted h5ad"].total_counts),
    ("per-cell genes detected", pcs["original .rds"].n_genes, pcs["converted h5ad"].n_genes),
    ("per-gene total counts", pg["original .rds"].total_counts, pg["converted h5ad"].total_counts),
    ("dotplot mean expression", dps["original .rds"].mean_expr, dps["converted h5ad"].mean_expr),
]
fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
for ax, (title, x, y) in zip(axes, pairs):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    ax.scatter(x, y, s=5, alpha=0.4, color="#4C72B0", linewidths=0)
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    ax.plot([lo, hi], [lo, hi], color="crimson", lw=1, ls="--")
    md = float(np.abs(x - y).max()) if len(x) == len(y) else float("nan")
    ax.set_title(f"{title}\nmax |Δ| = {md:.2e}", fontsize=9)
    ax.set_xlabel("original .rds"); ax.set_ylabel("converted h5ad")
suptitle(fig, "Original vs converted (identity line in red)")
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig(os.path.join(args.outdir, "delta.png"), dpi=130)
plt.close(fig)
print("wrote delta.png")
print(f"figures in {args.outdir}")
