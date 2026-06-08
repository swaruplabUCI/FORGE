#!/usr/bin/env python3
"""
merge_annotations.py

Merge cell type annotations into the peak matrix while preserving original barcodes.
Supports two annotation modes:
  - Marker mode (--annotations JSON): cluster-level marker-based → 'cell_type' column
  - CellTypist mode (--celltypist-h5ad): per-cell CellTypist predictions → 'celltypist_prediction' column
Metadata merge from sample CSV always runs in both modes.
"""

import argparse
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
import json


def build_palette(cats):
    """{category: hex} covering any N without silent gray overflow.

    N <= 102 → scanpy default_102 (hand-tuned contrast).
    N  > 102 → evenly sampled gist_ncar (full hue+value range).
    """
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import scanpy as sc
    n = len(cats)
    if n <= 102:
        colors = list(sc.pl.palettes.default_102[:n])
    else:
        cmap = plt.get_cmap("gist_ncar")
        colors = [mcolors.to_hex(cmap(i / max(n - 1, 1))) for i in range(n)]
    return dict(zip(cats, colors))

# Allow import of celltypist_broad_map.py from the same bin/ directory
sys.path.insert(0, str(Path(__file__).resolve().parent))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peak-matrix", required=True, help="Input peak matrix h5ad")
    p.add_argument("--metadata", required=True, help="Sample metadata CSV")
    p.add_argument("--output", default="peak_matrix_annotated.h5ad", help="Output h5ad")

    # Marker mode args (mutually exclusive with celltypist)
    p.add_argument("--annotations", default=None, help="Cell type annotations JSON (marker mode)")
    p.add_argument("--resolution", default=None, help="Leiden resolution key (marker mode, e.g. leiden_0_5)")

    # CellTypist mode arg
    p.add_argument("--celltypist-h5ad", default=None,
                   help="CellTypist-annotated gene_matrix h5ad (celltypist mode)")
    p.add_argument("--plot-dir", default=".",
                   help="Directory for output plots")

    # scATAnno mode arg
    p.add_argument("--scatanno-h5ad", default=None,
                   help="scATAnno-annotated h5ad (scatanno mode)")

    args = p.parse_args()

    # Validate: exactly one mode must be specified
    modes = sum([bool(args.annotations), bool(args.celltypist_h5ad), bool(args.scatanno_h5ad)])
    if modes != 1:
        p.error("Must specify exactly one of --annotations, --celltypist-h5ad, or --scatanno-h5ad")
    if args.annotations and not args.resolution:
        p.error("--resolution is required when using --annotations (marker mode)")

    return args


def extract_base_sample(sample_str, metadata_samples=None):
    """
    Extract base sample name from compound IDs by matching against known
    metadata sample names with progressively shorter prefixes.
    """
    s = str(sample_str)

    # Try progressively shorter prefixes to find a match in metadata
    if metadata_samples:
        parts = s.split('_')
        for i in range(len(parts), 0, -1):
            candidate = '_'.join(parts[:i])
            if candidate in metadata_samples:
                return candidate

    # No match found — return as-is
    return s


def apply_marker_annotations(peak_matrix, annotations_path, resolution_key):
    """Marker mode: map cluster IDs to cell types from JSON → 'cell_type' column."""
    print(f"\n{'='*70}")
    print("ANNOTATION MODE: Marker-based (super-user)")
    print(f"{'='*70}")

    print(f"Loading cell type annotations from {annotations_path}")
    with open(annotations_path, 'r') as f:
        annotations_dict = json.load(f)

    if resolution_key in annotations_dict:
        cluster_to_celltype = annotations_dict[resolution_key]
        print(f"Using nested structure with resolution key: {resolution_key}")
    else:
        print(f"Using flat structure (no resolution key found)")
        cluster_to_celltype = annotations_dict

    print(f"Cluster to cell type mapping (first 5):")
    for i, (k, v) in enumerate(list(cluster_to_celltype.items())[:5]):
        print(f"  {k} -> {v}")

    if resolution_key not in peak_matrix.obs.columns:
        raise ValueError(
            f"Resolution column '{resolution_key}' not found in peak_matrix.obs. "
            f"Available columns: {list(peak_matrix.obs.columns)}"
        )

    print("\nMapping clusters to cell types...")
    sample_value = list(cluster_to_celltype.values())[0]
    if isinstance(sample_value, dict):
        print("  Detected nested dictionary format with confidence scores")
        cluster_to_celltype_simple = {
            cluster: info['cell_type']
            for cluster, info in cluster_to_celltype.items()
        }
    else:
        print("  Detected simple string format")
        cluster_to_celltype_simple = cluster_to_celltype

    peak_matrix.obs['cell_type'] = (
        peak_matrix.obs[resolution_key]
        .astype(str)
        .map(cluster_to_celltype_simple)
    )

    unmapped = peak_matrix.obs['cell_type'].isna().sum()
    if unmapped > 0:
        print(f"  WARNING: {unmapped} cells could not be mapped to cell types")
        unique_clusters = peak_matrix.obs[resolution_key].unique()
        print(f"  Unique clusters in data: {sorted([str(c) for c in unique_clusters])}")
        print(f"  Clusters in mapping: {sorted(cluster_to_celltype_simple.keys())}")

    cell_type_counts = peak_matrix.obs['cell_type'].value_counts()
    print(f"  Cell type distribution:")
    print(cell_type_counts)


def apply_celltypist_annotations(peak_matrix, celltypist_h5ad_path):
    """CellTypist mode: transfer per-cell predictions from gene_matrix → peak_matrix."""
    from celltypist_broad_map import CELLTYPIST_BROAD_MAP

    print(f"\n{'='*70}")
    print("ANNOTATION MODE: CellTypist (default)")
    print(f"{'='*70}")

    print(f"Loading CellTypist-annotated gene matrix from {celltypist_h5ad_path}")
    gene_matrix = ad.read_h5ad(celltypist_h5ad_path)
    print(f"  Gene matrix: {gene_matrix.n_obs} cells")
    print(f"  Gene matrix obs columns: {list(gene_matrix.obs.columns)}")

    # Barcode alignment: gene_matrix and peak_matrix share barcodes from the same pipeline run
    peak_barcodes = set(peak_matrix.obs_names)
    gene_barcodes = set(gene_matrix.obs_names)
    common = peak_barcodes & gene_barcodes
    print(f"\n  Peak matrix barcodes: {len(peak_barcodes)}")
    print(f"  Gene matrix barcodes: {len(gene_barcodes)}")
    print(f"  Common barcodes: {len(common)}")

    if len(common) == 0:
        raise ValueError(
            "No common barcodes between peak_matrix and CellTypist gene_matrix. "
            "These should share barcodes from the same ATAC pipeline run."
        )

    if len(common) < len(peak_barcodes):
        missing = len(peak_barcodes) - len(common)
        print(f"  WARNING: {missing} peak_matrix cells not in gene_matrix — will be labeled 'Unknown'")

    # Transfer CellTypist columns by index alignment
    ct_columns = ['celltypist_prediction', 'cell_type_prediction', 'celltypist_confidence',
                  'celltypist_majority_voting']
    transferred = []
    for col in ct_columns:
        if col in gene_matrix.obs.columns:
            # Use .reindex to safely align by barcode, filling missing with appropriate defaults
            peak_matrix.obs[col] = gene_matrix.obs[col].reindex(peak_matrix.obs_names)
            transferred.append(col)

    print(f"\n  Transferred columns: {transferred}")

    if 'celltypist_prediction' not in transferred:
        raise ValueError(
            "CellTypist gene_matrix missing 'celltypist_prediction' column. "
            f"Available: {list(gene_matrix.obs.columns)}"
        )

    # Fill NaN predictions for cells not in gene_matrix
    col = peak_matrix.obs['celltypist_prediction']
    if hasattr(col, 'cat'):
        col = col.cat.add_categories('Unknown')
    peak_matrix.obs['celltypist_prediction'] = col.fillna('Unknown')

    # Apply hierarchical condensation: CELLTYPIST_BROAD_MAP → 'cell_type_broad'
    broad = peak_matrix.obs['celltypist_prediction'].map(CELLTYPIST_BROAD_MAP)
    if hasattr(broad, 'cat') and 'Other' not in broad.cat.categories:
        broad = broad.cat.add_categories('Other')
    peak_matrix.obs['cell_type_broad'] = broad.fillna('Other')

    # Report
    print(f"\n  Fine-grained predictions (celltypist_prediction):")
    fine_counts = peak_matrix.obs['celltypist_prediction'].value_counts()
    for ct in fine_counts.head(15).index:
        print(f"    {ct}: {fine_counts[ct]}")
    if len(fine_counts) > 15:
        print(f"    ... and {len(fine_counts) - 15} more types")

    print(f"\n  Broad categories (cell_type_broad):")
    broad_counts = peak_matrix.obs['cell_type_broad'].value_counts()
    for ct in broad_counts.index:
        print(f"    {ct}: {broad_counts[ct]}")


def apply_scatanno_annotations(peak_matrix, scatanno_h5ad_path):
    """scATAnno mode: transfer per-cell predictions from scATAnno h5ad → peak_matrix."""
    print(f"\n{'='*70}")
    print("ANNOTATION MODE: scATAnno (reference peak atlas)")
    print(f"{'='*70}")

    print(f"Loading scATAnno-annotated data from {scatanno_h5ad_path}")
    scatanno = ad.read_h5ad(scatanno_h5ad_path)
    print(f"  scATAnno data: {scatanno.n_obs} cells")
    print(f"  scATAnno obs columns: {list(scatanno.obs.columns)}")

    # Barcode alignment
    peak_barcodes = set(peak_matrix.obs_names)
    sa_barcodes = set(scatanno.obs_names)
    common = peak_barcodes & sa_barcodes
    print(f"\n  Peak matrix barcodes: {len(peak_barcodes)}")
    print(f"  scATAnno barcodes: {len(sa_barcodes)}")
    print(f"  Common barcodes: {len(common)}")

    if len(common) == 0:
        raise ValueError(
            "No common barcodes between peak_matrix and scATAnno data. "
            "Check barcode format alignment."
        )

    if len(common) < len(peak_barcodes):
        missing = len(peak_barcodes) - len(common)
        print(f"  WARNING: {missing} peak_matrix cells not in scATAnno — will be labeled 'Unknown'")

    # Transfer scATAnno columns by index alignment.
    # cell_type_broad is written by run_scatanno.py when the reference atlas
    # contains a class-level hierarchy column (class_label, cell_type_class, etc.).
    sa_columns = ['cell_type_prediction', 'pred_y', 'cluster_annotation',
                  'cell_type_broad',
                  'Uncertainty_Score', 'uncertainty_score_step1', 'uncertainty_score_step2',
                  '1.knn-based_celltype', '2.corrected_celltype']
    transferred = []
    for col in sa_columns:
        if col in scatanno.obs.columns:
            peak_matrix.obs[col] = scatanno.obs[col].reindex(peak_matrix.obs_names)
            transferred.append(col)

    print(f"\n  Transferred columns: {transferred}")

    # Ensure cell_type_prediction exists (primary column for downstream)
    if 'cell_type_prediction' not in transferred:
        if 'pred_y' in transferred:
            peak_matrix.obs['cell_type_prediction'] = peak_matrix.obs['pred_y']
            print("  Created cell_type_prediction from pred_y")
        else:
            raise ValueError(
                "scATAnno h5ad missing both 'cell_type_prediction' and 'pred_y'. "
                f"Available: {list(scatanno.obs.columns)}"
            )

    # Fill NaN for cells not in scATAnno
    for col in ['cell_type_prediction', 'pred_y', 'cluster_annotation']:
        if col in peak_matrix.obs.columns:
            series = peak_matrix.obs[col]
            if hasattr(series, 'cat'):
                if 'Unknown' not in series.cat.categories:
                    series = series.cat.add_categories('Unknown')
            peak_matrix.obs[col] = series.fillna('Unknown')

    # Populate cell_type_broad.
    # Priority 1: reference-derived broad classes written by run_scatanno.py when
    #             the reference atlas has a class-level hierarchy column.
    # Priority 2: identity — use cell_type_prediction directly so downstream
    #             modules always have a usable broad label.
    if ('cell_type_broad' in peak_matrix.obs.columns and
            peak_matrix.obs['cell_type_broad'].notna().any() and
            (peak_matrix.obs['cell_type_broad'].astype(str) != 'Unknown').any()):
        print("  Using reference-derived cell_type_broad from scATAnno output.")
    else:
        print("  No reference-derived cell_type_broad found — using cell_type_prediction as cell_type_broad.")
        peak_matrix.obs['cell_type_broad'] = (
            peak_matrix.obs['cell_type_prediction'].astype(str).astype('category')
        )

    # Report
    print(f"\n  Cell type predictions (cell_type_prediction):")
    ct_counts = peak_matrix.obs['cell_type_prediction'].value_counts()
    for ct in ct_counts.head(15).index:
        print(f"    {ct}: {ct_counts[ct]}")
    if len(ct_counts) > 15:
        print(f"    ... and {len(ct_counts) - 15} more types")

    if 'cell_type_broad' in peak_matrix.obs.columns:
        print(f"\n  Broad classes (cell_type_broad):")
        for ct, n in peak_matrix.obs['cell_type_broad'].value_counts().items():
            print(f"    {ct}: {n}")

    n_unknown = (peak_matrix.obs['cell_type_prediction'].str.lower() == 'unknown').sum()
    print(f"\n  Unknown fraction: {n_unknown}/{peak_matrix.n_obs} ({100*n_unknown/peak_matrix.n_obs:.1f}%)")

    if 'Uncertainty_Score' in peak_matrix.obs.columns:
        u = peak_matrix.obs['Uncertainty_Score'].astype(float)
        print(f"  Uncertainty: mean={u.mean():.3f}, median={u.median():.3f}")


def generate_celltype_plots(celltypist_h5ad_path, peak_matrix, plot_dir):
    """Generate UMAPs and dotplots using data-derived markers after CellTypist annotation."""
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from celltypist_broad_map import CELLTYPIST_BROAD_MAP

    print(f"\n{'='*70}")
    print("GENERATING CELL-TYPE VISUALIZATIONS (data-derived markers)")
    print(f"{'='*70}")

    plot_dir = Path(plot_dir)
    gene_matrix = ad.read_h5ad(celltypist_h5ad_path)

    # Transfer cell type labels from peak_matrix to gene_matrix by barcode
    common = set(gene_matrix.obs_names) & set(peak_matrix.obs_names)
    if len(common) == 0:
        print("  WARNING: No common barcodes — skipping plots")
        return

    for col in ['celltypist_prediction', 'cell_type_broad']:
        if col in peak_matrix.obs.columns:
            gene_matrix.obs[col] = peak_matrix.obs[col].reindex(gene_matrix.obs_names)
            if hasattr(gene_matrix.obs[col], 'cat') and 'Unknown' not in gene_matrix.obs[col].cat.categories:
                gene_matrix.obs[col] = gene_matrix.obs[col].cat.add_categories('Unknown')
            gene_matrix.obs[col] = gene_matrix.obs[col].fillna('Unknown')

    # Filter to cells with valid labels and sufficient group size
    groupby_col = 'cell_type_broad'
    if groupby_col not in gene_matrix.obs.columns:
        print(f"  WARNING: '{groupby_col}' not found — skipping plots")
        return

    valid_mask = ~gene_matrix.obs[groupby_col].isin(['Unknown', 'Other', ''])
    group_counts = gene_matrix.obs.loc[valid_mask, groupby_col].value_counts()
    valid_groups = group_counts[group_counts >= 10].index.tolist()

    if len(valid_groups) < 2:
        print(f"  WARNING: Only {len(valid_groups)} cell type groups with >= 10 cells — skipping DE")
        return

    gm_filtered = gene_matrix[valid_mask & gene_matrix.obs[groupby_col].isin(valid_groups)].copy()
    print(f"  Using {gm_filtered.n_obs} cells across {len(valid_groups)} broad cell types for DE")

    # Derive markers via rank_genes_groups
    try:
        sc.tl.rank_genes_groups(gm_filtered, groupby_col, method='wilcoxon')

        top_genes = []
        for group in valid_groups:
            genes = gm_filtered.uns['rank_genes_groups']['names'][str(group)][:5]
            top_genes.extend(genes)
        top_genes = list(dict.fromkeys(top_genes))

        print(f"  Derived {len(top_genes)} DE marker genes across {len(valid_groups)} cell types")
    except Exception as e:
        print(f"  WARNING: rank_genes_groups failed: {e}")
        return

    # Dotplot: DE markers grouped by broad cell type
    if top_genes:
        try:
            sc.pl.dotplot(
                gm_filtered,
                top_genes,
                groupby=groupby_col,
                standard_scale='var',
                swap_axes=True,
                show=False,
            )
            fig = plt.gcf()
            fig.savefig(plot_dir / "atac_dotplot_markers_by_cell_type.pdf",
                        dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved: atac_dotplot_markers_by_cell_type.pdf")
        except Exception as e:
            print(f"  WARNING: Could not generate dotplot: {e}")

    # UMAPs colored by cell type (if UMAP coords exist)
    # Adaptive compression: mirror plot_postscanvi.py logic
    if 'X_umap' in gene_matrix.obsm:
        min_cells_plot = 100

        # Determine best column via viability check on celltypist_prediction
        ct_col = 'celltypist_prediction'
        if ct_col in gene_matrix.obs.columns:
            ct_counts = gene_matrix.obs[ct_col].value_counts()
            n_total_types = len(ct_counts)
            n_passing = (ct_counts >= min_cells_plot).sum()
            viable_ratio = n_passing / n_total_types if n_total_types > 0 else 1.0

            if viable_ratio < 0.5 and 'cell_type_broad' in gene_matrix.obs.columns:
                print(f"  UMAP: {n_passing}/{n_total_types} types >= {min_cells_plot} cells "
                      f"({viable_ratio:.0%}) — compressing to broad categories")
                plot_col = 'cell_type_broad'
            else:
                plot_col = ct_col
        elif 'cell_type_broad' in gene_matrix.obs.columns:
            plot_col = 'cell_type_broad'
        else:
            plot_col = None

        if plot_col is not None:
            plot_counts = gene_matrix.obs[plot_col].value_counts()
            # Exclude 'Other'/'Unknown' from viability filter but keep them plotted
            valid_types = plot_counts[
                (plot_counts >= min_cells_plot) &
                (~plot_counts.index.isin(['Unknown']))
            ].index.tolist()
            gm_plot = gene_matrix[gene_matrix.obs[plot_col].isin(valid_types)].copy()
            n_plot_types = len(valid_types)
            categories = sorted([str(c) for c in valid_types])
            gm_plot.obs[plot_col] = pd.Categorical(
                gm_plot.obs[plot_col].astype(str), categories=categories, ordered=False,
            )
            palette = build_palette(categories)
            print(f"  UMAP: plotting {n_plot_types} cell types (>= {min_cells_plot} cells) using '{plot_col}'")

            try:
                if n_plot_types <= 20:
                    fig_width = max(12, 10 + (n_plot_types * 0.3))
                    fig, ax = plt.subplots(figsize=(fig_width, 8))
                    sc.pl.umap(
                        gm_plot, color=plot_col, ax=ax, show=False, palette=palette,
                        frameon=False, legend_loc='right margin', legend_fontsize=7,
                        title=f'ATAC Cell Types ({plot_col}, >= {min_cells_plot} cells)',
                    )
                else:
                    fig, ax = plt.subplots(figsize=(12, 10))
                    sc.pl.umap(
                        gm_plot, color=plot_col, ax=ax, show=False, palette=palette,
                        frameon=False, legend_loc='none',
                        title=f'ATAC Cell Types ({plot_col}, >= {min_cells_plot} cells)',
                    )
                    handles = [plt.Line2D([0], [0], marker='o', color='w',
                               markerfacecolor=palette[c], markersize=6, label=c)
                               for c in categories]
                    ncol = max(3, len(categories) // 15 + 1)
                    ax.legend(handles, categories, loc='upper center',
                              bbox_to_anchor=(0.5, -0.05), ncol=ncol, fontsize=6,
                              frameon=False, columnspacing=1.0, handletextpad=0.3)

                plt.tight_layout()
                fig.savefig(plot_dir / "atac_umap_cell_types.pdf", dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved: atac_umap_cell_types.pdf")
            except Exception as e:
                print(f"  WARNING: Could not generate UMAP: {e}")
        else:
            print("  No cell type columns available — skipping UMAP plots")
    else:
        print("  No X_umap in gene_matrix — skipping UMAP plots")


def generate_scatanno_plots(peak_matrix, plot_dir):
    """Generate ATAC UMAP colored by scATAnno cell type predictions."""
    import scanpy as sc
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print(f"\n{'='*70}")
    print("GENERATING CELL-TYPE VISUALIZATIONS (scATAnno)")
    print(f"{'='*70}")

    plot_dir = Path(plot_dir)

    if 'X_umap' not in peak_matrix.obsm:
        print("  No X_umap in peak_matrix — skipping UMAP plots")
        return

    plot_col = 'cell_type_prediction'
    if plot_col not in peak_matrix.obs.columns:
        print(f"  '{plot_col}' not found — skipping UMAP plots")
        return

    min_cells_plot = 100
    ct_counts = peak_matrix.obs[plot_col].value_counts()
    valid_types = ct_counts[
        (ct_counts >= min_cells_plot) &
        (~ct_counts.index.isin(['Unknown']))
    ].index.tolist()

    if len(valid_types) < 2:
        print(f"  Only {len(valid_types)} cell type groups with >= {min_cells_plot} cells — skipping")
        return

    pm_plot = peak_matrix[peak_matrix.obs[plot_col].isin(valid_types)].copy()
    n_plot_types = len(valid_types)
    categories = sorted([str(c) for c in valid_types])
    pm_plot.obs[plot_col] = pd.Categorical(
        pm_plot.obs[plot_col].astype(str), categories=categories, ordered=False,
    )
    palette = build_palette(categories)
    print(f"  UMAP: plotting {n_plot_types} cell types (>= {min_cells_plot} cells)")

    try:
        if n_plot_types <= 20:
            fig_width = max(12, 10 + (n_plot_types * 0.3))
            fig, ax = plt.subplots(figsize=(fig_width, 8))
            sc.pl.umap(
                pm_plot, color=plot_col, ax=ax, show=False, palette=palette,
                frameon=False, legend_loc='right margin', legend_fontsize=7,
                title=f'ATAC Cell Types (scATAnno, >= {min_cells_plot} cells)',
            )
        else:
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.umap(
                pm_plot, color=plot_col, ax=ax, show=False, palette=palette,
                frameon=False, legend_loc='none',
                title=f'ATAC Cell Types (scATAnno, >= {min_cells_plot} cells)',
            )
            handles = [plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=palette[c], markersize=6, label=c)
                       for c in categories]
            ncol = max(3, len(categories) // 15 + 1)
            ax.legend(handles, categories, loc='upper center',
                      bbox_to_anchor=(0.5, -0.05), ncol=ncol, fontsize=6,
                      frameon=False, columnspacing=1.0, handletextpad=0.3)

        plt.tight_layout()
        fig.savefig(plot_dir / "atac_umap_cell_types.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: atac_umap_cell_types.pdf")
    except Exception as e:
        print(f"  WARNING: Could not generate UMAP: {e}")


def merge_sample_metadata(peak_matrix, metadata_path, original_barcodes):
    """Merge sample-level metadata onto peak_matrix. Shared by both modes."""
    print(f"\n{'='*70}")
    print("MERGING SAMPLE METADATA")
    print(f"{'='*70}")

    print(f"Loading metadata from {metadata_path}")
    metadata = pd.read_csv(metadata_path)
    print(f"Metadata shape: {metadata.shape}")
    print(f"Metadata columns: {list(metadata.columns)}")

    # Identify the correct sample column in metadata
    sample_col_metadata = None
    for col_name in ['sample_id', 'sample', 'Sample', 'sample_name']:
        if col_name in metadata.columns:
            sample_col_metadata = col_name
            print(f"Found sample column in metadata: '{sample_col_metadata}'")
            print(f"Sample values: {sorted(metadata[sample_col_metadata].unique())}")
            break

    if sample_col_metadata is None:
        raise ValueError(f"Could not find sample column in metadata. Available: {list(metadata.columns)}")

    # Check if peak_matrix has 'sample' column
    if 'sample' not in peak_matrix.obs.columns:
        raise ValueError(f"'sample' column not found in peak_matrix.obs. Available: {list(peak_matrix.obs.columns)}")

    print(f"\nPeak matrix sample values (raw): {sorted(peak_matrix.obs['sample'].unique())}")

    # CRITICAL FIX: Extract base sample name from compound IDs
    print("\nExtracting base sample names from peak matrix...")
    meta_sample_set = set(metadata[sample_col_metadata].unique())
    peak_matrix.obs['base_sample'] = peak_matrix.obs['sample'].apply(
        lambda x: extract_base_sample(x, metadata_samples=meta_sample_set)
    )
    print(f"Extracted base samples: {sorted(peak_matrix.obs['base_sample'].unique())}")

    # Check overlap between peak matrix and metadata
    pm_samples = set(peak_matrix.obs['base_sample'].dropna().unique())
    meta_samples = set(metadata[sample_col_metadata].unique())
    overlap = pm_samples.intersection(meta_samples)

    print(f"\n Merge compatibility check:")
    print(f"  Samples in peak matrix: {len(pm_samples)}")
    print(f"  Samples in metadata: {len(meta_samples)}")
    print(f"   Overlapping samples: {len(overlap)}")
    print(f"   Missing in metadata: {sorted(pm_samples - meta_samples)}")
    print(f"   Missing in peak matrix: {sorted(meta_samples - pm_samples)}")

    if len(overlap) == 0:
        raise ValueError(
            "\n MERGE FAILURE: No overlapping samples!\n"
            f"  Peak matrix: {sorted(pm_samples)}\n"
            f"  Metadata: {sorted(meta_samples)}\n"
            "  Check if sample naming conventions match."
        )

    print(f"   Found {len(overlap)} overlapping samples")

    # Rename metadata sample column to 'base_sample' for merging
    if sample_col_metadata != 'base_sample':
        print(f"\nRenaming metadata column '{sample_col_metadata}' to 'base_sample' for merge")
        metadata = metadata.rename(columns={sample_col_metadata: 'base_sample'})

    # Check for duplicate samples in metadata (should not happen after filtering)
    duplicates = metadata['base_sample'].duplicated()
    if duplicates.any():
        print(f"\nWARNING: Found {duplicates.sum()} duplicate samples in metadata!")
        dup_samples = metadata[duplicates]['base_sample'].unique()
        print(f"Duplicate samples: {dup_samples}")
        # Keep only first occurrence
        metadata = metadata.drop_duplicates(subset=['base_sample'], keep='first')
        print(f"Kept first occurrence of duplicates. New shape: {metadata.shape}")

    print("\nMerging metadata on 'base_sample' column...")

    # Store original index to restore later
    original_index = peak_matrix.obs.index.copy()

    # Reset index to avoid issues during merge
    peak_matrix.obs = peak_matrix.obs.reset_index(drop=True)

    # Merge metadata (left join to keep all cells)
    peak_matrix.obs = peak_matrix.obs.merge(
        metadata,
        on='base_sample',
        how='left',
        suffixes=('', '_meta')
    )

    # Restore original barcodes as index
    peak_matrix.obs.index = original_index
    peak_matrix.obs.index.name = None

    # Drop the temporary base_sample column
    if 'base_sample' in peak_matrix.obs.columns:
        peak_matrix.obs.drop(['base_sample'], axis=1, inplace=True)

    # Verify the shape is preserved
    if len(peak_matrix.obs) != len(original_barcodes):
        raise RuntimeError(
            f"CRITICAL: Number of cells changed during merge! "
            f"Original: {len(original_barcodes)}, Current: {len(peak_matrix.obs)}"
        )

    # Verify condition column was added
    if 'condition' not in peak_matrix.obs.columns:
        print(f"\nWARNING: 'condition' column missing after merge!")
        print(f"  Available columns: {list(peak_matrix.obs.columns)}")
    else:
        print(f"\n Successfully added 'condition' column")
        print(f"  Condition values: {peak_matrix.obs['condition'].unique()}")
        print(f"  Condition counts: {peak_matrix.obs['condition'].value_counts().to_dict()}")

    # Final verification
    print("\nVerification:")
    print(f"  Final shape: {peak_matrix.shape}")
    print(f"  Original shape cells: {len(original_barcodes)}")
    print(f"  First 5 obs_names: {list(peak_matrix.obs_names[:5])}")
    print(f"  Original first 5: {list(original_barcodes[:5])}")

    # Check that barcodes are preserved
    if not np.array_equal(peak_matrix.obs_names, original_barcodes):
        raise RuntimeError("CRITICAL: Barcodes were not preserved correctly!")

    print(" Barcodes successfully preserved")


def main():
    args = parse_args()

    print(f"Loading peak matrix from {args.peak_matrix}")
    peak_matrix = ad.read_h5ad(args.peak_matrix)

    print(f"Original peak matrix:")
    print(f"  Shape: {peak_matrix.shape}")
    print(f"  First 5 obs_names: {list(peak_matrix.obs_names[:5])}")
    print(f"  obs columns: {list(peak_matrix.obs.columns)}")

    # CRITICAL: Store original barcodes BEFORE any manipulation
    original_barcodes = peak_matrix.obs_names.to_numpy().copy()
    print(f"\nStored {len(original_barcodes)} original barcodes")

    # --- Apply cell type annotations (mode-dependent) ---
    if args.celltypist_h5ad:
        apply_celltypist_annotations(peak_matrix, args.celltypist_h5ad)
    elif args.scatanno_h5ad:
        apply_scatanno_annotations(peak_matrix, args.scatanno_h5ad)
    else:
        apply_marker_annotations(peak_matrix, args.annotations, args.resolution)

    # --- Generate cell-type plots ---
    if args.celltypist_h5ad:
        generate_celltype_plots(args.celltypist_h5ad, peak_matrix, args.plot_dir)
    elif args.scatanno_h5ad and args.plot_dir:
        generate_scatanno_plots(peak_matrix, args.plot_dir)

    # --- Merge sample metadata (always) ---
    merge_sample_metadata(peak_matrix, args.metadata, original_barcodes)

    # Final summary
    print(f"\nFinal peak matrix:")
    print(f"  Shape: {peak_matrix.shape}")
    print(f"  obs columns: {list(peak_matrix.obs.columns)}")
    if 'cell_type' in peak_matrix.obs.columns:
        print(f"  Cell types (marker): {list(peak_matrix.obs['cell_type'].unique())}")
    if 'celltypist_prediction' in peak_matrix.obs.columns:
        n_types = peak_matrix.obs['celltypist_prediction'].nunique()
        print(f"  CellTypist predictions: {n_types} types")
    if 'cell_type_broad' in peak_matrix.obs.columns:
        print(f"  Broad categories: {list(peak_matrix.obs['cell_type_broad'].unique())}")
    if 'condition' in peak_matrix.obs.columns:
        print(f"  Conditions: {list(peak_matrix.obs['condition'].unique())}")

    # FIX-21: Convert NaN in object columns to empty strings for h5ad compat.
    for col in peak_matrix.obs.columns:
        if peak_matrix.obs[col].dtype == object:
            peak_matrix.obs[col] = peak_matrix.obs[col].fillna('').astype(str)

    # Save
    print(f"\nSaving annotated peak matrix to {args.output}")
    peak_matrix.write_h5ad(args.output, compression='gzip')
    print("Done!")


if __name__ == "__main__":
    main()
