#!/usr/bin/env python3
"""
plot_scenicplus_results.py — SCENIC+ downstream visualization

Generates publication-quality figures from SCENIC+ outputs:
  1. eRegulon AUCell heatmap-dotplot (RSS + expression per cell type)
  2. eRegulon AUCell UMAP embeddings
  3. eRegulon correlation heatmap (Jaccard overlap)
  4. Top eRegulons per cell type (ranked bar chart)
  5. TF-gene target count summary

Inputs:
  --scplus-mudata     scplusmdata.h5mu from SCENIC+
  --aucell-direct     AUCell_direct.h5mu
  --eregulon-tsv      eRegulon_direct.tsv (TF→Region→Gene triplets)
  --cell-type-key     Column name for cell type annotations (default: cell_type)
  --top-n             Number of top eRegulons per cell type to show (default: 10)
  --dpi               Output DPI (default: 200)
  --outdir            Output directory (default: .)
"""
import argparse
import os
import sys
import warnings

warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scanpy as sc
import mudata as md


def compute_rss(auc_matrix, cell_type_labels):
    """Compute Regulon Specificity Score (RSS) per eRegulon per cell type.

    RSS = 1 - Jensen-Shannon divergence between eRegulon AUC distribution
    and a reference distribution that is 1 for the focal cell type and 0 elsewhere.
    Higher RSS = more specific to that cell type.
    """
    from scipy.spatial.distance import jensenshannon
    from scipy.special import rel_entr

    cell_types = sorted(cell_type_labels.unique())
    regulon_names = auc_matrix.columns if hasattr(auc_matrix, 'columns') else [f'R{i}' for i in range(auc_matrix.shape[1])]

    rss_df = pd.DataFrame(index=regulon_names, columns=cell_types, dtype=float)

    for ct in cell_types:
        # Reference: uniform within cell type, zero outside
        ref = np.zeros(len(cell_type_labels))
        ref[cell_type_labels == ct] = 1.0
        ref = ref / ref.sum()

        for i, reg in enumerate(regulon_names):
            # eRegulon AUC distribution across cells
            auc_vals = auc_matrix[:, i] if hasattr(auc_matrix, 'toarray') else (
                auc_matrix.iloc[:, i].values if hasattr(auc_matrix, 'iloc') else auc_matrix[:, i]
            )
            if hasattr(auc_vals, 'toarray'):
                auc_vals = auc_vals.toarray().flatten()
            auc_dist = np.abs(auc_vals)
            if auc_dist.sum() == 0:
                rss_df.loc[reg, ct] = 0
                continue
            auc_dist = auc_dist / auc_dist.sum()

            # RSS = 1 - JSD
            jsd = jensenshannon(auc_dist, ref)
            rss_df.loc[reg, ct] = 1 - jsd

    return rss_df.astype(float)


def plot_rss_heatmap(rss_df, top_n=10, save_path=None, dpi=200):
    """Plot RSS heatmap showing top N most specific eRegulons per cell type."""
    # Select top N per cell type
    top_regulons = set()
    for ct in rss_df.columns:
        top = rss_df[ct].nlargest(top_n).index.tolist()
        top_regulons.update(top)

    rss_top = rss_df.loc[sorted(top_regulons)]

    # Extract TF name from eRegulon name (e.g., "ATF3_direct_+/+_(1212g)" -> "ATF3")
    tf_labels = [r.split('_')[0] for r in rss_top.index]
    rss_top.index = [f"{tf} ({r.split('(')[1].rstrip(')')}" if '(' in r else tf
                     for tf, r in zip(tf_labels, rss_top.index)]

    fig, ax = plt.subplots(figsize=(max(10, len(rss_top.columns) * 1.2),
                                     max(8, len(rss_top) * 0.35)))
    im = ax.imshow(rss_top.values, aspect='auto', cmap='YlOrRd')
    ax.set_xticks(range(len(rss_top.columns)))
    ax.set_xticklabels(rss_top.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(rss_top)))
    ax.set_yticklabels(rss_top.index, fontsize=8)
    ax.set_title(f'eRegulon Specificity Score (top {top_n} per cell type)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='RSS', shrink=0.5)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        fig.savefig(save_path.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved RSS heatmap: {save_path}")
    plt.close(fig)
    return rss_top


def plot_aucell_umap(scplus_mudata, cell_type_key, top_regulons, save_path=None, dpi=200):
    """Plot UMAP colored by top eRegulon AUCell scores and cell types."""
    rna = scplus_mudata.mod['scRNA_counts']

    # Compute UMAP if not present
    if 'X_umap' not in rna.obsm:
        print("  Computing UMAP on RNA...")
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
        sc.pp.highly_variable_genes(rna, n_top_genes=2000)
        sc.pp.pca(rna, n_comps=30)
        sc.pp.neighbors(rna)
        sc.tl.umap(rna)

    # Cell type UMAP
    fig, axes = plt.subplots(1, 1, figsize=(10, 8))
    sc.pl.umap(rna, color=cell_type_key, ax=axes, show=False, frameon=False,
               title='Cell Types', legend_loc='on data', legend_fontsize=8)

    ct_path = save_path.replace('.png', '_celltypes.png') if save_path else None
    if ct_path:
        fig.savefig(ct_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved cell type UMAP: {ct_path}")
    plt.close(fig)

    # eRegulon AUCell UMAPs (top regulons)
    auc_mod = scplus_mudata.mod.get('direct_gene_based_AUC')
    if auc_mod is None:
        print("  [WARN] No direct_gene_based_AUC modality — skipping AUCell UMAPs")
        return

    # Transfer UMAP coordinates from RNA to AUC
    auc_ad = auc_mod.copy()
    auc_ad.obsm['X_umap'] = rna.obsm['X_umap']

    # Select top regulons that exist in var_names
    avail = [r for r in top_regulons if r in auc_ad.var_names][:12]
    if not avail:
        print("  [WARN] No top regulons found in AUCell var_names")
        return

    n_cols = min(4, len(avail))
    n_rows = (len(avail) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = np.atleast_2d(axes)

    for idx, reg in enumerate(avail):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]
        sc.pl.umap(auc_ad, color=reg, ax=ax, show=False, frameon=False,
                   title=reg.split('_')[0], cmap='viridis', vmin=0)

    # Hide unused axes
    for idx in range(len(avail), n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row, col].axis('off')

    plt.suptitle('eRegulon Activity (AUCell Scores)', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        fig.savefig(save_path.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved AUCell UMAPs: {save_path}")
    plt.close(fig)


def plot_eregulon_correlation(auc_matrix, regulon_names, save_path=None, dpi=200):
    """Plot correlation heatmap between eRegulon AUCell scores."""
    if hasattr(auc_matrix, 'toarray'):
        auc_matrix = auc_matrix.toarray()
    if hasattr(auc_matrix, 'values'):
        auc_matrix = auc_matrix.values

    corr = np.corrcoef(auc_matrix.T)

    # Label with TF names
    tf_labels = [r.split('_')[0] for r in regulon_names]

    fig, ax = plt.subplots(figsize=(max(10, len(tf_labels) * 0.3),
                                     max(8, len(tf_labels) * 0.3)))
    im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    ax.set_xticks(range(len(tf_labels)))
    ax.set_xticklabels(tf_labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(tf_labels)))
    ax.set_yticklabels(tf_labels, fontsize=6)
    ax.set_title('eRegulon AUCell Correlation', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Pearson r', shrink=0.5)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        fig.savefig(save_path.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved correlation heatmap: {save_path}")
    plt.close(fig)


def plot_top_regulons_per_celltype(rss_df, top_n=5, save_path=None, dpi=200):
    """Plot horizontal bar chart of top eRegulons per cell type by RSS."""
    cell_types = rss_df.columns.tolist()
    n_ct = len(cell_types)
    n_cols = min(4, n_ct)
    n_rows = (n_ct + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    for idx, ct in enumerate(cell_types):
        ax = axes[idx]
        top = rss_df[ct].nlargest(top_n)
        tf_labels = [r.split('_')[0] for r in top.index]
        ax.barh(range(top_n), top.values, color='steelblue')
        ax.set_yticks(range(top_n))
        ax.set_yticklabels(tf_labels, fontsize=9)
        ax.set_xlabel('RSS')
        ax.set_title(ct, fontsize=11, fontweight='bold')
        ax.invert_yaxis()

    for idx in range(n_ct, len(axes)):
        axes[idx].axis('off')

    plt.suptitle(f'Top {top_n} eRegulons per Cell Type (by RSS)', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        fig.savefig(save_path.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved top regulons bar chart: {save_path}")
    plt.close(fig)


def plot_tf_target_summary(eregulon_tsv_path, save_path=None, dpi=200):
    """Plot summary of TF-gene target counts from eRegulon TSV."""
    df = pd.read_csv(eregulon_tsv_path, sep='\t')

    # Count unique genes per TF for activating (+/+) regulons
    activating = df[df['eRegulon_name'].str.contains(r'\+/\+', na=False)]
    tf_gene_counts = activating.groupby('TF')['Gene'].nunique().sort_values(ascending=False)

    if len(tf_gene_counts) == 0:
        print("  [WARN] No activating eRegulons found — skipping TF target summary")
        return

    top_tfs = tf_gene_counts.head(30)

    fig, ax = plt.subplots(figsize=(10, max(6, len(top_tfs) * 0.35)))
    ax.barh(range(len(top_tfs)), top_tfs.values, color='darkorange')
    ax.set_yticks(range(len(top_tfs)))
    ax.set_yticklabels(top_tfs.index, fontsize=9)
    ax.set_xlabel('Number of Target Genes')
    ax.set_title('Top TFs by Target Gene Count (Activating eRegulons)', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        fig.savefig(save_path.replace('.png', '.pdf'), dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"  Saved TF target summary: {save_path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description='SCENIC+ visualization')
    parser.add_argument('--scplus-mudata', required=True, help='scplusmdata.h5mu')
    parser.add_argument('--aucell-direct', required=True, help='AUCell_direct.h5mu')
    parser.add_argument('--eregulon-tsv', required=True, help='eRegulon_direct.tsv')
    parser.add_argument('--cell-type-key', default='cell_type')
    parser.add_argument('--top-n', type=int, default=10)
    parser.add_argument('--dpi', type=int, default=200)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=== SCENIC+ Visualization ===")

    # Load data
    print("Loading scplusmdata.h5mu...")
    scplus = md.read_h5mu(args.scplus_mudata)
    rna = scplus.mod['scRNA_counts']

    # Detect cell type column
    ct_key = args.cell_type_key
    if ct_key not in rna.obs.columns:
        for fallback in ['celltypist_prediction', 'scanvi_prediction', 'cell_type']:
            if fallback in rna.obs.columns:
                ct_key = fallback
                break
    print(f"  Cell type key: {ct_key}")
    print(f"  Cell types: {rna.obs[ct_key].nunique()}")
    print(f"  Cells: {rna.n_obs}")

    # Load AUCell scores
    print("Loading AUCell_direct.h5mu...")
    aucell = md.read_h5mu(args.aucell_direct)
    auc_gene = aucell.mod['Gene_based']
    print(f"  eRegulons (gene-based): {auc_gene.n_vars}")

    # Get AUC matrix + cell type labels
    auc_X = auc_gene.X
    if hasattr(auc_X, 'toarray'):
        auc_X = auc_X.toarray()
    cell_types = rna.obs[ct_key].values
    regulon_names = list(auc_gene.var_names)

    # 1. Compute RSS
    print("\n[1/5] Computing Regulon Specificity Scores...")
    rss_df = compute_rss(
        pd.DataFrame(auc_X, columns=regulon_names, index=auc_gene.obs_names),
        pd.Series(cell_types, index=rna.obs_names)
    )
    rss_df.to_csv(os.path.join(args.outdir, 'eregulon_rss.csv'))
    print(f"  RSS matrix: {rss_df.shape}")

    # 2. RSS heatmap
    print("\n[2/5] Plotting RSS heatmap...")
    rss_top = plot_rss_heatmap(
        rss_df, top_n=args.top_n,
        save_path=os.path.join(args.outdir, 'eregulon_rss_heatmap.png'),
        dpi=args.dpi
    )

    # 3. Top regulons per cell type
    print("\n[3/5] Plotting top eRegulons per cell type...")
    plot_top_regulons_per_celltype(
        rss_df, top_n=5,
        save_path=os.path.join(args.outdir, 'eregulon_top_per_celltype.png'),
        dpi=args.dpi
    )

    # 4. eRegulon correlation heatmap
    print("\n[4/5] Plotting eRegulon correlation...")
    plot_eregulon_correlation(
        auc_X, regulon_names,
        save_path=os.path.join(args.outdir, 'eregulon_correlation.png'),
        dpi=args.dpi
    )

    # 5. AUCell UMAPs
    print("\n[5/5] Plotting AUCell UMAPs...")
    # Pick top regulons from RSS (most specific ones overall)
    top_regs = []
    for ct in rss_df.columns:
        top_regs.extend(rss_df[ct].nlargest(3).index.tolist())
    top_regs = list(dict.fromkeys(top_regs))[:12]  # Deduplicate, keep order

    plot_aucell_umap(
        scplus, ct_key, top_regs,
        save_path=os.path.join(args.outdir, 'eregulon_aucell_umap.png'),
        dpi=args.dpi
    )

    # 6. TF target gene summary
    print("\nBonus: TF target gene summary...")
    plot_tf_target_summary(
        args.eregulon_tsv,
        save_path=os.path.join(args.outdir, 'tf_target_summary.png'),
        dpi=args.dpi
    )

    print(f"\n=== SCENIC+ visualization complete. Outputs in {args.outdir} ===")


if __name__ == '__main__':
    main()
