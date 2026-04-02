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

TISSUE_MARKERS = {
    'pbmc': {
        "T cell (CD4+)": ["CD3D", "CD3E", "CD4", "IL7R"],
        "T cell (CD8+)": ["CD3D", "CD3E", "CD8A", "CD8B"],
        "B cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
        "NK cell": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
        "Monocyte (CD14+)": ["CD14", "LYZ", "S100A9", "S100A8"],
        "Monocyte (CD16+)": ["FCGR3A", "MS4A7", "LILRA5"],
        "Dendritic cell": ["FCER1A", "CST3", "CLEC10A"],
        "Platelet": ["PPBP", "PF4"],
    },
    'brain': {
        "Oligodendrocyte": ["MBP", "PLP1", "MAG", "MOG"],
        "Astrocyte": ["GFAP", "S100B", "ALDH1L1", "SLC1A3"],
        "Microglia": ["P2RY12", "CX3CR1", "TMEM119"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "KDR"],
        "Pericyte/VLMC": ["PDGFRB", "PDGFRA", "COL1A1", "COL1A2"],
        "Pan-neuronal": ["SYT1", "SNAP25"],
        "Excitatory (IT/ET)": ["SLC17A7", "SATB2", "RORB", "FEZF2", "BCL11B"],
        "Inhibitory (Pvalb/Sst/Vip)": ["GAD1", "GAD2", "PVALB", "SST", "VIP", "LAMP5"],
        "OPC": ["PDGFRA", "CSPG4", "OLIG2"],
    },
}


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

    # Plot 3: Cell type predictions – syntax aligned with interactive tests
    if cell_type_key in adata.obs.columns:
        if n_celltypes <= 20:
            # Few enough cell types: legend beside the plot
            fig_width = max(12, 10 + (n_celltypes * 0.3))
            fig, ax = plt.subplots(figsize=(fig_width, 8))
            sc.pl.umap(
                adata,
                color=cell_type_key,
                ax=ax,
                show=False,
                title=f'Cell Type Predictions ({cell_type_key})',
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
            # Many cell types: plot UMAP on top, legend below
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.umap(
                adata,
                color=cell_type_key,
                ax=ax,
                show=False,
                title=f'Cell Type Predictions ({cell_type_key})',
                legend_loc='none',
                frameon=False
            )
            # Build a standalone legend underneath
            handles, labels = ax.get_legend_handles_labels()
            if not handles:
                # scanpy sometimes puts legend on the figure, not the axes
                categories = adata.obs[cell_type_key].cat.categories
                palette = dict(zip(categories, adata.uns.get(f'{cell_type_key}_colors',
                               sc.pl.palettes.default_102[:len(categories)])))
                handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=palette[c], markersize=6, label=c) for c in categories]
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
    # CELL-TYPE MARKER DOTPLOT BY SCANVI PREDICTION
    # ========================================

    if cell_type_key in adata.obs.columns:
        print(f"\nGenerating marker dotplot by {cell_type_key}...")

        # Tissue-aware markers with case-insensitive matching
        key_markers = TISSUE_MARKERS.get(tissue_type, TISSUE_MARKERS.get('pbmc', {}))

        # Case-insensitive map: lower -> true name
        var_lower_to_true = {v.lower(): v for v in adata.var_names}
        markers_in_data = {}
        for cell_type, markers in key_markers.items():
            present_true = [
                var_lower_to_true[m.lower()]
                for m in markers
                if m.lower() in var_lower_to_true
            ]
            if present_true:
                markers_in_data[cell_type] = present_true

        if markers_in_data:
            sc.pl.dotplot(
                adata,
                markers_in_data,
                groupby=cell_type_key,
                standard_scale="var",
                swap_axes=True,
                show=False,
            )
            fig = plt.gcf()
            fig.savefig(
                os.path.join(output_dir, "dotplot_markers_by_scanvi_prediction.png"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)
            print("  Saved: dotplot_markers_by_scanvi_prediction.png")
        else:
            print("  No specified markers found in this dataset; skipping dotplot.")

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

    output_h5ad = os.path.join(output_dir, 'annotated_with_scanvi_clustering.h5ad')
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
    print("    - dotplot_markers_by_scanvi_prediction.png")
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

    # Adaptive compression: if <50% of types are viable, switch to broad categories
    use_broad = viable_ratio < 0.5 and len(CELLTYPIST_BROAD_MAP) > 0
    if use_broad:
        print(f"Adaptive compression: {n_passing}/{n_total} types viable ({viable_ratio:.0%}) — switching to broad categories")
        adata.obs['cell_type_broad'] = (
            adata.obs[cell_type_key]
            .map(CELLTYPIST_BROAD_MAP)
            .fillna('Progenitors/Other')
        )
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
    if args.cell_type_key and args.cell_type_key in adata.obs.columns:
        cell_type_key = args.cell_type_key
    else:
        if args.cell_type_key:
            print(f"WARNING: Requested cell_type_key '{args.cell_type_key}' not found in adata.obs — auto-detecting")
        cell_type_key = None
        for candidate in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_prediction']:
            if candidate in adata.obs.columns:
                cell_type_key = candidate
                break
        if cell_type_key is None:
            print("ERROR: No cell type prediction column found in adata.obs. "
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

    # Re-save h5ad if cell_type_broad was added (needed by downstream CONVERT_H5AD_TO_SEURAT)
    if 'cell_type_broad' in adata.obs.columns:
        output_h5ad = os.path.join(args.output_dir, 'annotated_with_scanvi_clustering.h5ad')
        adata.write(output_h5ad, compression='gzip')
        print(f"Re-saved h5ad with cell_type_broad column: {output_h5ad}")

if __name__ == "__main__":
    main()
