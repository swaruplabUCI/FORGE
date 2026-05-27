#!/usr/bin/env python3
"""
Post-Integration RNA Plotting with Memory Optimization
Generates UMAPs, dotplots, and QC statistics at multiple resolutions
Uses igraph backend for memory efficiency
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import sys
import gc
from pathlib import Path

# Memory optimization
import warnings
warnings.filterwarnings('ignore')

# Allow import of celltypist_broad_map.py from the same bin/ directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from celltypist_broad_map import CELLTYPIST_BROAD_MAP
from h5ad_compat import sanitize_adata



def create_post_integration_plots(adata, output_dir, resolutions=[0.5, 5.0], cell_type_key='scanvi_prediction', tissue_type='pbmc', alt_cell_type_key=None):
    """
    Generate post-integration visualizations with SCANVI embeddings
    """

    print(f"\n{'='*70}")
    print("POST-INTEGRATION PLOTTING")
    print(f"{'='*70}")

    # Find SCANVI embedding
    scanvi_rep = None
    for key in ['X_scANVI', 'X_scanvi_100e', 'X_scanvi']:
        if key in adata.obsm:
            scanvi_rep = key
            break

    if scanvi_rep is None:
        print(' Warning: No SCANVI embedding found, falling back to PCA')
        scanvi_rep = 'X_pca' if 'X_pca' in adata.obsm else None

    print(f"Using embedding: {scanvi_rep}")

    # Recompute neighbors using SCANVI embeddings
    print("Computing neighbors...")
    sc.pp.neighbors(adata, use_rep=scanvi_rep, n_neighbors=30)

    # Leiden clustering at multiple resolutions (CRITICAL: use igraph)
    for res in resolutions:
        res_str = str(res).replace('.', '_')
        cluster_key = f'leiden_scanvi_{res_str}'

        print(f"\nClustering at resolution {res}...")
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=cluster_key,
            flavor='igraph',      # CRITICAL: Use igraph backend
            n_iterations=-1,      # Run to completion
            directed=False        # Required for igraph
        )

        n_clusters = len(adata.obs[cluster_key].unique())
        print(f"  Found {n_clusters} clusters")

        # Force garbage collection
        gc.collect()

    # Compute UMAP
    print("\nComputing UMAP...")
    sc.tl.umap(adata, min_dist=0.3)

    # ========================================
    # INDIVIDUAL UMAP PLOTS
    # ========================================

    print("\nGenerating UMAP visualizations...")

    # Count samples / cell types for dynamic sizing
    n_samples = len(adata.obs['sample'].unique()) if 'sample' in adata.obs.columns else 0
    n_celltypes = len(adata.obs[cell_type_key].unique()) if cell_type_key in adata.obs.columns else 0

    # Plot 1: Leiden 0.5 clustering
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='leiden_scanvi_0_5',
        ax=ax,
        show=False,
        title='Leiden Clustering (Resolution 0.5)',
        legend_loc='on data',
        frameon=False
    )
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'umap_leiden_0.5.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: umap_leiden_0.5.png")

    # Plot 2: Leiden 5.0 clustering
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color='leiden_scanvi_5_0',
        ax=ax,
        show=False,
        title='Leiden Clustering (Resolution 5.0)',
        legend_loc='on data',
        frameon=False
    )
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'umap_leiden_5.0.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: umap_leiden_5.0.png")

    # Plot 3: Cell type predictions — adaptive compression mirroring hdWGCNA logic
    if cell_type_key in adata.obs.columns:
        min_cells_plot = 100
        cell_type_counts = adata.obs[cell_type_key].value_counts()
        n_total_types = len(cell_type_counts)
        n_passing = (cell_type_counts >= min_cells_plot).sum()
        viable_ratio = n_passing / n_total_types if n_total_types > 0 else 1.0

        # Adaptive compression: only compress to broad categories when viable AND the map
        # actually applies.  Guard: if >80% of cells collapse to 'Progenitors/Other' the
        # map is inapplicable for this tissue (e.g. mouse brain labels not in PBMC map).
        use_broad = viable_ratio < 0.5 and len(CELLTYPIST_BROAD_MAP) > 0
        if use_broad:
            plot_key = 'cell_type_broad'
            if plot_key not in adata.obs.columns:
                adata.obs[plot_key] = (
                    adata.obs[cell_type_key]
                    .map(CELLTYPIST_BROAD_MAP)
                    .fillna('Progenitors/Other')
                )
            broad_counts = adata.obs[plot_key].value_counts()
            progenitor_frac = broad_counts.get('Progenitors/Other', 0) / len(adata.obs)
            if progenitor_frac > 0.8:
                print(f"  UMAP: broad map collapsed {progenitor_frac:.0%} → 'Progenitors/Other' "
                      f"(map inapplicable for this tissue); falling back to fine-grained '{cell_type_key}'")
                use_broad = False
            else:
                print(f"  UMAP: {n_passing}/{n_total_types} types >= {min_cells_plot} cells ({viable_ratio:.0%}) — compressing to broad categories")

        if use_broad:
            plot_key = 'cell_type_broad'
            plot_counts = broad_counts
        else:
            plot_key = cell_type_key
            plot_counts = cell_type_counts
            # Lower floor when falling back so small non-PBMC datasets remain informative
            min_cells_plot = 50

        # Filter to types with >= min_cells
        valid_types = plot_counts[plot_counts >= min_cells_plot].index.tolist()
        adata_plot = adata[adata.obs[plot_key].isin(valid_types)].copy()
        n_plot_types = len(valid_types)
        print(f"  UMAP: plotting {n_plot_types} cell types (>= {min_cells_plot} cells) using '{plot_key}'")

        if n_plot_types <= 20:
            fig_width = max(12, 10 + (n_plot_types * 0.3))
            fig, ax = plt.subplots(figsize=(fig_width, 8))
            sc.pl.umap(
                adata_plot,
                color=plot_key,
                ax=ax,
                show=False,
                title=f'Cell Types ({plot_key}, >= {min_cells_plot} cells)',
                legend_loc='right margin',
                legend_fontsize=7,
                frameon=False
            )
            plt.tight_layout()
            fig.savefig(
                os.path.join(output_dir, 'umap_cell_types.png'),
                dpi=300,
                bbox_inches='tight'
            )
            plt.close(fig)
        else:
            # Inject full palette before sc.pl.umap so N>102 types get distinct colors
            _cats_p = (adata_plot.obs[plot_key].cat.categories.tolist()
                       if hasattr(adata_plot.obs[plot_key], 'cat')
                       else sorted(adata_plot.obs[plot_key].unique()))
            if len(adata_plot.uns.get(f'{plot_key}_colors', [])) < len(_cats_p):
                import matplotlib.colors as _mc
                _n_p = len(_cats_p)
                _cmap_p = plt.get_cmap("gist_ncar" if _n_p > 102 else "turbo" if _n_p > 20 else "tab20")
                adata_plot.uns[f'{plot_key}_colors'] = [
                    _mc.to_hex(_cmap_p(i / max(_n_p - 1, 1))) for i in range(_n_p)
                ]
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.umap(
                adata_plot,
                color=plot_key,
                ax=ax,
                show=False,
                title=f'Cell Types ({plot_key}, >= {min_cells_plot} cells)',
                legend_loc='none',
                frameon=False
            )
            handles, labels = ax.get_legend_handles_labels()
            if not handles:
                categories = adata_plot.obs[plot_key].cat.categories if hasattr(adata_plot.obs[plot_key], 'cat') else sorted(adata_plot.obs[plot_key].unique())
                palette = dict(zip(categories, adata_plot.uns.get(f'{plot_key}_colors',
                               sc.pl.palettes.default_102[:len(categories)])))
                handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=palette.get(c, '#888888'), markersize=6, label=c) for c in categories]
                labels = list(categories)
            ncol = max(3, len(labels) // 15 + 1)
            ax.legend(
                handles, labels,
                loc='upper center',
                bbox_to_anchor=(0.5, -0.05),
                ncol=ncol,
                fontsize=6,
                frameon=False,
                columnspacing=1.0,
                handletextpad=0.3
            )
            fig.savefig(
                os.path.join(output_dir, 'umap_cell_types.png'),
                dpi=300,
                bbox_inches='tight'
            )
            plt.close(fig)
        print("  Saved: umap_cell_types.png")

    # Plot 3b: Secondary cell type UMAP (alternate prediction method)
    if alt_cell_type_key and alt_cell_type_key in adata.obs.columns:
        n_alt = len(adata.obs[alt_cell_type_key].unique())
        fig_width = max(12, 10 + (n_alt * 0.3)) if n_alt <= 20 else 12
        # Inject full palette before sc.pl.umap (CellTypist yields 200+ types; default_102 truncates)
        if len(adata.uns.get(f'{alt_cell_type_key}_colors', [])) < n_alt:
            import matplotlib.colors as _mc
            _cats_alt = (adata.obs[alt_cell_type_key].cat.categories.tolist()
                         if hasattr(adata.obs[alt_cell_type_key], 'cat')
                         else sorted(adata.obs[alt_cell_type_key].unique()))
            _cmap_alt = plt.get_cmap("gist_ncar" if n_alt > 102 else "turbo" if n_alt > 20 else "tab20")
            adata.uns[f'{alt_cell_type_key}_colors'] = [
                _mc.to_hex(_cmap_alt(i / max(len(_cats_alt) - 1, 1))) for i in range(len(_cats_alt))
            ]
        fig, ax = plt.subplots(figsize=(fig_width, 8))
        sc.pl.umap(
            adata,
            color=alt_cell_type_key,
            ax=ax,
            show=False,
            title=f'Cell Type Predictions ({alt_cell_type_key})',
            legend_loc='right margin' if n_alt <= 20 else 'none',
            legend_fontsize=7,
            frameon=False
        )
        if n_alt > 20:
            handles, labels = ax.get_legend_handles_labels()
            if not handles:
                categories = adata.obs[alt_cell_type_key].cat.categories
                palette = dict(zip(categories, adata.uns.get(f'{alt_cell_type_key}_colors',
                               sc.pl.palettes.default_102[:len(categories)])))
                handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=palette[c], markersize=6, label=c) for c in categories]
                labels = list(categories)
            ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
                      ncol=max(3, len(labels) // 15 + 1), fontsize=6, frameon=False)
        plt.tight_layout()
        alt_fname = f'umap_cell_types_{alt_cell_type_key}.png'
        fig.savefig(os.path.join(output_dir, alt_fname), dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {alt_fname}")

    # Plot 4: Batch distribution (post-integration)
    if 'batch' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata,
            color='batch',
            ax=ax,
            show=False,
            title='Batch Distribution (Post-Integration)',
            legend_loc='right margin',
            frameon=False
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'umap_batch.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: umap_batch.png")

    # Plot 5: Sample distribution (no legend – too many samples)
    if 'sample' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata,
            color='sample',
            ax=ax,
            show=False,
            title=f'Sample Distribution ({n_samples} samples)',
            legend_loc=None,  # No legend on plot
            size=3,
            frameon=False
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'umap_sample.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: umap_sample.png")

    # ========================================
    # CELL-TYPE MARKER DOTPLOT BY SCANVI PREDICTION (data-derived)
    # ========================================

    if cell_type_key in adata.obs.columns:
        print(f"\nDeriving marker genes by {cell_type_key} via rank_genes_groups...")

        try:
            # Filter to valid groups with >= 10 cells
            group_counts = adata.obs[cell_type_key].value_counts()
            valid_groups = group_counts[group_counts >= 10].index.tolist()
            adata_filtered = adata[adata.obs[cell_type_key].isin(valid_groups)].copy()

            sc.tl.rank_genes_groups(adata_filtered, cell_type_key, method='wilcoxon')

            top_genes = []
            for group in valid_groups:
                genes = adata_filtered.uns['rank_genes_groups']['names'][str(group)][:5]
                top_genes.extend(genes)
            top_genes = list(dict.fromkeys(top_genes))

            if top_genes:
                sc.pl.dotplot(
                    adata_filtered,
                    top_genes,
                    groupby=cell_type_key,
                    standard_scale="var",
                    swap_axes=True,
                    show=False,
                )
                fig = plt.gcf()
                dotplot_fname = f"dotplot_markers_by_{cell_type_key}.png"
                fig.savefig(
                    os.path.join(output_dir, dotplot_fname),
                    dpi=300,
                    bbox_inches="tight",
                )
                plt.close(fig)
                print(f"  Saved: {dotplot_fname} ({len(top_genes)} DE genes, {len(valid_groups)} types)")
            else:
                print("  No DE genes found; skipping dotplot.")

        except Exception as e:
            print(f"  Warning: Could not generate cell-type marker dotplot: {e}")

    # ========================================
    # MARKER GENE DOTPLOTS (CLUSTER-BASED)
    # ========================================

    print("\nComputing marker genes for Leiden clusters...")

    for res in resolutions:
        res_str = str(res).replace('.', '_')
        cluster_key = f'leiden_scanvi_{res_str}'

        try:
            sc.tl.rank_genes_groups(adata, cluster_key, method='wilcoxon')

            # Get top genes
            top_genes = []
            for group in adata.obs[cluster_key].unique()[:10]:  # Top 10 clusters
                genes = adata.uns['rank_genes_groups']['names'][str(group)][:5]
                top_genes.extend(genes)
            top_genes = list(set(top_genes))[:30]

            # Dotplot
            if top_genes:
                sc.pl.dotplot(
                    adata,
                    top_genes,
                    groupby=cluster_key,
                    save=f'_leiden_scanvi_{res_str}.pdf',
                    show=False
                )
                print(f"  Saved: dotplot_leiden_scanvi_{res_str}.pdf")

        except Exception as e:
            print(f"  Warning: Could not generate markers for {cluster_key}: {e}")

        gc.collect()

    # ========================================
    # INTEGRATION METRICS (3 SEPARATE PLOTS)
    # ========================================

    print("\nGenerating integration quality metrics...")

    # Metric 1: Batch mixing entropy
    if 'batch' in adata.obs.columns and len(adata.obs['batch'].unique()) > 1:
        from scipy.stats import entropy

        # Recalculate neighbors for entropy
        sc.pp.neighbors(adata, use_rep=scanvi_rep, n_neighbors=15)

        batch_entropy = []
        for i in range(adata.n_obs):
            neighbors = adata.obsp['connectivities'][i].nonzero()[1][:15]
            if len(neighbors) > 0:
                batch_counts = adata.obs.iloc[neighbors]['batch'].value_counts()
                batch_probs = batch_counts / batch_counts.sum()
                batch_entropy.append(entropy(batch_probs))
            else:
                batch_entropy.append(0)

        adata.obs['batch_entropy_scanvi'] = batch_entropy

        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata,
            color='batch_entropy_scanvi',
            ax=ax,
            show=False,
            title='Batch Mixing Quality (Higher = Better)',
            cmap='viridis',
            frameon=False
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'integration_batch_mixing.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: integration_batch_mixing.png")

    # Metric 2: Cell type composition by sample
    if 'sample' in adata.obs.columns and cell_type_key in adata.obs.columns:
        ct_comp = adata.obs.groupby(['sample', cell_type_key]).size().unstack(fill_value=0)
        ct_comp_pct = ct_comp.div(ct_comp.sum(axis=1), axis=0) * 100

        # Get top cell types
        top_celltypes = adata.obs[cell_type_key].value_counts().head(15).index

        # Dynamic sizing: ~0.35 inches per sample
        fig_height = max(20, n_samples * 0.35)
        fig_width = max(18, len(top_celltypes) * 0.8)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ct_comp_pct[top_celltypes].plot(kind='barh', stacked=True, ax=ax, width=0.85)
        ax.set_xlabel('Percentage', fontsize=12)
        ax.set_ylabel('Sample', fontsize=12)
        ax.set_title(f'Cell Type Composition Across {n_samples} Samples', fontsize=14)
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7,
                  title='Cell Type', title_fontsize=8)
        ax.tick_params(axis='y', labelsize=6)
        ax.tick_params(axis='x', labelsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'integration_celltype_composition.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: integration_celltype_composition.png")

    # Metric 3: QC metrics distribution
    if 'n_genes_by_counts' in adata.obs.columns:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Genes per cell
        axes[0, 0].hist(adata.obs['n_genes_by_counts'], bins=50,
                        alpha=0.7, edgecolor='black')
        median_genes = adata.obs['n_genes_by_counts'].median()
        axes[0, 0].axvline(median_genes, color='r', linestyle='--',
                           label=f'Median: {median_genes:.0f}')
        axes[0, 0].set_xlabel('Number of Genes')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Gene Count Distribution')
        axes[0, 0].legend()

        # UMI per cell
        if 'total_counts' in adata.obs.columns:
            axes[0, 1].hist(adata.obs['total_counts'], bins=50,
                            alpha=0.7, edgecolor='black')
            median_umi = adata.obs['total_counts'].median()
            axes[0, 1].axvline(median_umi, color='r', linestyle='--',
                               label=f'Median: {median_umi:.0f}')
            axes[0, 1].set_xlabel('Total UMI Counts')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('UMI Count Distribution')
            axes[0, 1].legend()
        else:
            axes[0, 1].axis('off')

        # MT percent
        if 'pct_counts_mt' in adata.obs.columns:
            axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50,
                            alpha=0.7, edgecolor='black')
            median_mt = adata.obs['pct_counts_mt'].median()
            axes[1, 0].axvline(median_mt, color='r', linestyle='--',
                               label=f'Median: {median_mt:.2f}%')
            axes[1, 0].set_xlabel('Mitochondrial %')
            axes[1, 0].set_ylabel('Frequency')
            axes[1, 0].set_title('Mitochondrial % Distribution')
            axes[1, 0].legend()
        else:
            axes[1, 0].axis('off')

        # Cells per cluster (leiden 0.5)
        cluster_counts = adata.obs['leiden_scanvi_0_5'].value_counts().sort_index()
        axes[1, 1].bar(range(len(cluster_counts)), cluster_counts.values,
                       alpha=0.7, edgecolor='black')
        axes[1, 1].set_xlabel('Cluster')
        axes[1, 1].set_ylabel('Number of Cells')
        axes[1, 1].set_title('Cells per Cluster (Leiden 0.5)')
        axes[1, 1].set_xticks(range(len(cluster_counts)))
        axes[1, 1].set_xticklabels(cluster_counts.index, rotation=45, ha='right')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'integration_qc_distributions.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("  Saved: integration_qc_distributions.png")

    # ========================================
    # SAVE UPDATED ADATA
    # ========================================

    output_h5ad = os.path.join(output_dir, 'annotated_with_celltype.h5ad')
    sanitize_adata(adata, output_h5ad)
    adata.write(output_h5ad, compression='gzip')
    print(f"\n✓ Updated h5ad saved to: {output_h5ad}")

    print(f"\n{'='*70}")
    print("POST-INTEGRATION PLOTTING COMPLETE")
    print(f"{'='*70}")
    print("\nGenerated plots:")
    print("  UMAPs:")
    print("    - umap_leiden_0.5.png")
    print("    - umap_leiden_5.0.png")
    print("    - umap_cell_types.png")
    print("    - umap_batch.png")
    print("    - umap_sample.png")
    print("  Marker dotplots:")
    print(f"    - dotplot_markers_by_{cell_type_key}.png")
    print("    - dotplot_leiden_scanvi_*.pdf (cluster-based markers)")
    print("  Integration metrics:")
    print("    - integration_batch_mixing.png")
    print("    - integration_celltype_composition.png")
    print("    - integration_qc_distributions.png")
    print(f"{'='*70}")

def extract_cell_types_for_hdwgcna(adata, output_dir, cell_type_key='scanvi_prediction', min_cells=100):
    """
    Extract cell types with sufficient cells for hdWGCNA.
    Adaptively compresses to broad categories when <50% of fine-grained
    types pass the min_cells threshold (e.g. small datasets).
    """
    print(f"\n{'='*70}")
    print("EXTRACTING CELL TYPES FOR hdWGCNA")
    print(f"{'='*70}")

    # Count fine-grained types
    cell_type_counts = adata.obs[cell_type_key].value_counts()
    n_total = len(cell_type_counts)
    n_passing = (cell_type_counts >= min_cells).sum()
    viable_ratio = n_passing / n_total if n_total > 0 else 1.0

    print(f"Fine-grained types: {n_total}, passing min_cells={min_cells}: {n_passing} ({viable_ratio:.0%})")

    # Adaptive compression: if <50% of types are viable, switch to broad categories.
    # Guard: if >80% of cells collapse to 'Progenitors/Other' the map is inapplicable
    # for this tissue (e.g. mouse brain labels not in PBMC map) — revert to fine-grained.
    use_broad = viable_ratio < 0.5 and len(CELLTYPIST_BROAD_MAP) > 0
    if use_broad:
        adata.obs['cell_type_broad'] = (
            adata.obs[cell_type_key]
            .map(CELLTYPIST_BROAD_MAP)
            .fillna('Progenitors/Other')
        )
        broad_counts_ct = adata.obs['cell_type_broad'].value_counts()
        progenitor_frac_ct = broad_counts_ct.get('Progenitors/Other', 0) / len(adata.obs)
        if progenitor_frac_ct > 0.8:
            print(f"Adaptive compression: broad map collapsed {progenitor_frac_ct:.0%} → 'Progenitors/Other' "
                  f"(map inapplicable); reverting to fine-grained '{cell_type_key}'")
            use_broad = False
        else:
            print(f"Adaptive compression: {n_passing}/{n_total} types viable ({viable_ratio:.0%}) — switching to broad categories")

    if use_broad:
        active_key = 'cell_type_broad'
        cell_type_counts = adata.obs[active_key].value_counts()
    else:
        print(f"Using fine-grained types: {n_passing}/{n_total} viable ({viable_ratio:.0%})")
        active_key = cell_type_key

    # Filter by minimum cells
    valid_types = cell_type_counts[cell_type_counts >= min_cells].index.tolist()

    print(f"\nUsing key '{active_key}': {len(valid_types)} types with >= {min_cells} cells:")
    for ct in valid_types:
        print(f"  - {ct}: {cell_type_counts[ct]} cells")

    # Write valid cell types (one per line)
    cell_types_file = os.path.join(output_dir, 'hdwgcna_cell_types.txt')
    with open(cell_types_file, 'w') as f:
        for ct in valid_types:
            f.write(f"{ct}\n")

    # Write counts summary
    counts_df = pd.DataFrame({
        'cell_type': cell_type_counts.index,
        'n_cells': cell_type_counts.values,
        'included': cell_type_counts.index.isin(valid_types)
    })
    counts_file = os.path.join(output_dir, 'hdwgcna_cell_type_counts.csv')
    counts_df.to_csv(counts_file, index=False)

    # Write which key hdWGCNA should use
    key_file = os.path.join(output_dir, 'hdwgcna_cell_type_key.txt')
    with open(key_file, 'w') as f:
        f.write(active_key)

    print(f"\n✓ Cell types written to: {cell_types_file}")
    print(f"✓ Counts summary written to: {counts_file}")
    print(f"✓ Active key: {active_key} (written to {key_file})")
    print(f"{'='*70}\n")

    return cell_types_file, counts_file

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Annotated h5ad file')
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--resolutions', nargs='+', type=float,
                        default=[0.5, 5.0])
    parser.add_argument('--hdwgcna_min_cells', type=int, default=100,
                        help='Minimum cells per type for hdWGCNA')
    parser.add_argument('--cell_type_key', type=str, default=None,
                        help='Column name for cell type predictions (auto-detected if not provided)')
    parser.add_argument('--tissue_type', type=str, default='pbmc',
                        choices=['pbmc', 'brain'],
                        help='Tissue type for marker dotplot (default: pbmc)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading data from {args.input}...")
    adata = sc.read_h5ad(args.input)

    print(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Observations: {list(adata.obs.columns)}")
    print(f"Embeddings: {list(adata.obsm.keys())}")

    # Auto-detect cell type key — prefer explicit arg from pipeline, fallback to auto-detect
    # Validation: skip columns where all values are 'Unknown' (no real annotations)
    def _has_real_labels(col):
        vals = adata.obs[col].dropna().unique()
        return len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown')

    cell_type_key = None
    if args.cell_type_key and args.cell_type_key in adata.obs.columns and _has_real_labels(args.cell_type_key):
        cell_type_key = args.cell_type_key
    else:
        if args.cell_type_key:
            print(f"WARNING: Requested cell_type_key '{args.cell_type_key}' not usable — auto-detecting")
        for candidate in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_prediction']:
            if candidate in adata.obs.columns and _has_real_labels(candidate):
                cell_type_key = candidate
                break
        if cell_type_key is None:
            print("ERROR: No cell type prediction column with real labels found in adata.obs. "
                  "Expected one of: scanvi_prediction, celltypist_prediction, cell_type_prediction. "
                  f"Available columns: {list(adata.obs.columns)}")
            sys.exit(1)
    print(f"Using cell type key: {cell_type_key}")

    # Identify alternate prediction key (for secondary UMAP when both methods ran)
    alt_cell_type_key = None
    all_ct_keys = [k for k in ['scanvi_prediction', 'celltypist_prediction'] if k in adata.obs.columns and k != cell_type_key]
    if all_ct_keys:
        alt_cell_type_key = all_ct_keys[0]
        print(f"Alternate cell type key available: {alt_cell_type_key}")

    create_post_integration_plots(adata, args.output_dir, args.resolutions,
                                   cell_type_key=cell_type_key, tissue_type=args.tissue_type,
                                   alt_cell_type_key=alt_cell_type_key)

    # Extract cell types for hdWGCNA
    extract_cell_types_for_hdwgcna(
        adata,
        args.output_dir,
        cell_type_key=cell_type_key,
        min_cells=args.hdwgcna_min_cells
    )

    # Stamp canonical 'cell_type' + provenance on obs so downstream modules
    # (BUILD_MUDATA, MULTIVI_GAP_FILL) don't have to juggle tool-specific names.
    _CT_CANDIDATES = [
        ('scanvi_prediction', 'scanvi'),
        ('celltypist_prediction', 'celltypist'),
        ('cell_type_prediction', 'celltypist'),
        ('cell_type_marker', 'marker'),
    ]
    chosen_col, chosen_src = None, None
    for col, src in _CT_CANDIDATES:
        if col in adata.obs.columns and _has_real_labels(col):
            chosen_col, chosen_src = col, src
            break
    if chosen_col is not None:
        adata.obs['cell_type'] = adata.obs[chosen_col].values
        adata.obs['cell_type_source'] = chosen_src
        print(f"Stamped canonical 'cell_type' from '{chosen_col}' (source={chosen_src})")

    # Always write the canonical output (renamed from annotated_with_scanvi_clustering.h5ad
    # to reflect that the column is tool-agnostic). Keep the old filename as a symlink
    # only if something downstream still references it.
    output_h5ad = os.path.join(args.output_dir, 'annotated_with_celltype.h5ad')
    sanitize_adata(adata, output_h5ad)
    adata.write(output_h5ad, compression='gzip')
    print(f"Saved annotated h5ad: {output_h5ad}")

if __name__ == "__main__":
    main()
