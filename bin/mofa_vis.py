#!/usr/bin/env python3
"""
MOFA+ Visualization Pipeline - Multiome-Compatible
Properly extracts cell type annotations from RNA modality
"""

import pandas as pd
import scanpy as sc
import mudata as md
import numpy as np
from pathlib import Path
import mofax as mfx
from h5ad_compat import sanitize_mudata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import f_oneway
import warnings
import json
import argparse
import sys
from datetime import datetime

warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(
        description='MOFA+ visualization and downstream analysis'
    )
    
    parser.add_argument('--mofa_model', type=str, required=True,
                       help='MOFA model file (.hdf5)')
    parser.add_argument('--mudata_file', type=str, required=True,
                       help='MuData file (.h5mu)')
    parser.add_argument('--metadata_file', type=str, default=None,
                       help='MOFA metadata JSON (optional)')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory')
    
    return parser.parse_args()

def extract_lane_from_sample(sample_str):
    """Extract lane from sample_id like 'L1_Donor_0_july'"""
    try:
        parts = str(sample_str).split('_')
        if len(parts) > 0 and parts[0].startswith('L'):
            return parts[0]  # Returns 'L1', 'L2', etc.
        return 'Unknown'
    except:
        return 'Unknown'

def extract_batch_from_sample(sample_str):
    """Extract batch/timepoint from sample_id like 'L1_Donor_0_july'"""
    try:
        parts = str(sample_str).split('_')
        if len(parts) >= 4:
            return parts[3]  # Returns 'july', 'nov', etc.
        return 'Unknown'
    except:
        return 'Unknown'

def main():
    args = parse_args()
    
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    FIGURES_DIR = OUTPUT_DIR / "figures"
    FIGURES_DIR.mkdir(exist_ok=True)
    
    print("="*80)
    print("MOFA+ VISUALIZATION - COMPLETE PIPELINE")
    print("="*80)
    
    # ==========================================================================
    # STEP 1: LOAD MOFA MODEL
    # ==========================================================================
    print(f"\nLoading MOFA model: {args.mofa_model}")
    model = mfx.mofa_model(str(args.mofa_model))
    
    print(f"""
Model dimensions:
  Cells: {model.shape[0]}
  Features: {model.shape[1]}
  Groups: {', '.join(model.groups)}
  Views: {', '.join(model.views)}
""")
    
    # Get factors and weights
    factors_array = model.get_factors(df=False)
    factors_df = model.get_factors(df=True)
    weights_array = model.get_weights(df=False)
    weights_df = model.get_weights(df=True)
    
    n_factors = factors_array.shape[1]
    print(f"  Factors shape: {factors_array.shape}")
    print(f"  Number of factors: {n_factors}")
    
    # ==========================================================================
    # STEP 2: COMPREHENSIVE WEIGHT VISUALIZATIONS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 2: WEIGHT VISUALIZATIONS")
    print("="*80)
    
    # 2.1 Weights for all factors
    print("\nGenerating weights plot (all factors)...")
    fig = plt.figure(figsize=(12, 8))
    mfx.plot_weights(model, n_features=10)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_weights_all_factors.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_weights_all_factors.png")
    
    # 2.2 Ranked weights for individual factors
    print(f"\nGenerating ranked weights for {min(3, n_factors)} factors...")
    for factor_idx in range(min(3, n_factors)):
        plt.figure(figsize=(6, 8))
        mfx.plot_weights_ranked(model, factor=factor_idx, n_features=10,
                               y_repel_coef=0.04, x_rank_offset=-150)
        plt.title(f"Factor {factor_idx+1} - Top Weighted Features")
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / f'mofa_weights_ranked_factor{factor_idx+1}.png',
                    dpi=300, bbox_inches='tight')
        plt.close()
    print(f"✓ Saved ranked weights for {min(3, n_factors)} factors")
    
    # 2.3 Weights heatmap
    print("\nGenerating weights heatmap...")
    plt.figure(figsize=(12, 8))
    mfx.plot_weights_heatmap(model, n_features=10, 
                             factors=list(range(min(7, n_factors))),
                             xticklabels_size=6, w_abs=True, 
                             cmap="viridis", cluster_factors=False)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_weights_heatmap.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_weights_heatmap.png")
    
    # 2.4 Weights correlation
    print("\nGenerating weights correlation...")
    plt.figure(figsize=(8, 7))
    mfx.plot_weights_correlation(model)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_weights_correlation.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_weights_correlation.png")
    
    # ==========================================================================
    # STEP 3: VARIANCE EXPLAINED
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 3: VARIANCE EXPLAINED")
    print("="*80)
    
    r2_df = model.get_r2()
    view_totals = r2_df.groupby('View')['R2'].sum()
    
    print("\nVariance explained by view:")
    for view, total in view_totals.items():
        view_name = 'RNA (view0)' if 'view0' in view else 'ATAC (view1)'
        print(f"  {view_name}: {total:.2f}%")
    
    # Variance plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Total per view
    view_labels = ['RNA (view0)', 'ATAC (view1)']
    axes[0].bar(range(len(view_totals)), view_totals.values, 
                color=['#E64B35', '#4DBBD5'])
    axes[0].set_xticks(range(len(view_totals)))
    axes[0].set_xticklabels(view_labels)
    axes[0].set_ylabel('Total Variance Explained (%)')
    axes[0].set_title('MOFA+ Total Variance Explained')
    axes[0].set_ylim([0, 100])
    
    # Add value labels
    for i, v in enumerate(view_totals.values):
        axes[0].text(i, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontsize=12)
    
    # Per factor
    r2_pivot = r2_df.pivot(index='Factor', columns='View', values='R2')
    r2_pivot.columns = ['RNA (view0)', 'ATAC (view1)']
    r2_pivot.plot(kind='bar', ax=axes[1], color=['#E64B35', '#4DBBD5'])
    axes[1].set_ylabel('Variance Explained (%)')
    axes[1].set_xlabel('MOFA Factor')
    axes[1].legend(title='View')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_variance_explained.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_variance_explained.png")
    
    # Save variance table
    variance_csv = OUTPUT_DIR / "mofa_variance_explained.csv"
    r2_df.to_csv(variance_csv, index=False)
    print(f"✓ Saved: {variance_csv}")
    
    # ==========================================================================
    # STEP 4: FACTOR VISUALIZATIONS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 4: FACTOR VISUALIZATIONS")
    print("="*80)
    
    # 4.1 Factors correlation
    print("\nGenerating factors correlation heatmap...")
    plt.figure(figsize=(8, 7))
    mfx.plot_factors_correlation(model)
    plt.title("Factor Correlation (Pearson r)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_factors_correlation.png',
                dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_factors_correlation.png")
    
    # 4.2 Factors scatter plots (first 6 factors)
    print("\nGenerating factor scatter plots...")
    n_scatter = min(6, n_factors)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    plot_idx = 0
    for i in range(n_scatter - 1):
        for j in range(i + 1, n_scatter):
            if plot_idx >= 6:
                break
            axes[plot_idx].scatter(factors_array[:, i], factors_array[:, j],
                                  alpha=0.5, s=10)
            axes[plot_idx].set_xlabel(f'Factor {i+1}')
            axes[plot_idx].set_ylabel(f'Factor {j+1}')
            axes[plot_idx].set_title(f'Factor {i+1} vs {j+1}')
            plot_idx += 1
            if plot_idx >= 6:
                break
    
    # Hide unused subplots
    for idx in range(plot_idx, 6):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'mofa_factors_scatter.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: mofa_factors_scatter.png")
    
    # ==========================================================================
    # STEP 5: LOAD MUDATA AND EXTRACT METADATA
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 5: LOADING MUDATA & EXTRACTING METADATA")
    print("="*80)
    
    print(f"\nLoading MuData: {args.mudata_file}")
    mdata = md.read_h5mu(args.mudata_file)
    
    # Match cells (downsampling was done in integration)
    if mdata.n_obs > factors_array.shape[0]:
        print(f"  Subsetting to {factors_array.shape[0]} cells")
        mdata = mdata[:factors_array.shape[0], :].copy()
    
    mdata.obsm['X_mofa'] = factors_array
    print(f"✓ Added MOFA factors to MuData")
    print(f"  Working with {mdata.n_obs:,} cells")
    
    # ==========================================================================
    # CRITICAL: Extract cell type from RNA modality
    # ==========================================================================
    print("\n" + "="*80)
    print("EXTRACTING METADATA FROM MODALITIES")
    print("="*80)
    
    # Cell type from RNA predictions (CellTypist or scANVI)
    cell_type_source = None
    for key in ['rna:celltypist_prediction', 'rna:cell_type_prediction',
                'rna:scanvi_prediction', 'celltypist_prediction',
                'cell_type_prediction', 'scanvi_prediction',
                'atac:scanvi_prediction']:
        if key in mdata.obs.columns:
            vals = mdata.obs[key].dropna().unique()
            if len(vals) == 1 and vals[0] == 'Unknown':
                continue  # skip all-Unknown columns
            mdata.obs['cell_type'] = mdata.obs[key]
            cell_type_source = key
            print(f"Extracted cell types from {key}")
            print(f"  Unique cell types: {mdata.obs['cell_type'].nunique()}")
            break
    if cell_type_source is None:
        print("\n" + "!"*80)
        print("WARNING: No cell type prediction column found in any modality.")
        print("  Searched: rna:celltypist_prediction, rna:cell_type_prediction,")
        print("            rna:scanvi_prediction, celltypist_prediction,")
        print("            cell_type_prediction, scanvi_prediction, atac:scanvi_prediction")
        print(f"  Available obs columns: {[c for c in mdata.obs.columns if 'type' in c.lower() or 'predict' in c.lower() or 'label' in c.lower()]}")
        print("  Cell type UMAP will be SKIPPED.")
        print("!"*80)
    
    # Extract lane from sample_id
    if 'sample_id' in mdata.obs.columns:
        mdata.obs['lane'] = mdata.obs['sample_id'].apply(extract_lane_from_sample)
        print(f"✓ Extracted lane from sample_id")
        print(f"  Unique lanes: {mdata.obs['lane'].nunique()}")
        
        # Extract batch/timepoint
        mdata.obs['batch'] = mdata.obs['sample_id'].apply(extract_batch_from_sample)
        print(f"✓ Extracted batch from sample_id")
        print(f"  Unique batches: {mdata.obs['batch'].nunique()}")
    
    # Print summary
    print("\nFinal metadata summary:")
    for col in ['cell_type', 'batch', 'lane', 'sample_id']:
        if col in mdata.obs.columns:
            n_unique = mdata.obs[col].nunique()
            print(f"  {col}: {n_unique} unique values")
    
    # ==========================================================================
    # STEP 6: COMPUTE UMAP AND CLUSTERING
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 6: DOWNSTREAM ANALYSIS")
    print("="*80)
    
    print("\nComputing neighbors on MOFA factors...")
    sc.pp.neighbors(mdata, use_rep='X_mofa', n_neighbors=15)
    
    print("Computing UMAP...")
    sc.tl.umap(mdata, min_dist=0.2, spread=1.0, random_state=42)
    
    print("Leiden clustering...")
    sc.tl.leiden(mdata, key_added='leiden_mofa', resolution=0.5)
    print(f"✓ Identified {mdata.obs['leiden_mofa'].nunique()} clusters")
    
    # ==========================================================================
    # STEP 7: COMPREHENSIVE UMAP VISUALIZATIONS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 7: GENERATING UMAP VISUALIZATIONS")
    print("="*80)
    
    # Get UMAP coordinates
    umap_coords = mdata.obsm['X_umap']
    
    # 7.1 UMAP by cell type (RNA-based annotations)
    if 'cell_type' in mdata.obs.columns:
        print(f"\nPlotting UMAP by cell type...")

        cell_types = mdata.obs['cell_type'].astype(str)
        unique_types = sorted(cell_types.unique())
        n_types = len(unique_types)

        # Use tab20 for <=20 types, turbo colormap for more
        if n_types <= 20:
            colors = plt.cm.tab20(np.linspace(0, 1, n_types))
        else:
            colors = plt.cm.turbo(np.linspace(0.05, 0.95, n_types))
        color_map = dict(zip(unique_types, colors))

        # Dynamic layout: legend below for many types, right for few
        if n_types <= 20:
            fig, ax = plt.subplots(figsize=(12, 8))
        else:
            fig, ax = plt.subplots(figsize=(12, 10))

        for ct in unique_types:
            mask = cell_types == ct
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                      c=[color_map[ct]], label=ct, s=5, alpha=0.7)

        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        source_desc = cell_type_source if cell_type_source else 'unknown'
        ax.set_title(f'MOFA UMAP — Cell types from RNA ({source_desc})')

        if n_types <= 20:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        else:
            ncol = max(3, n_types // 10)
            ax.legend(bbox_to_anchor=(0.5, -0.12), loc='upper center',
                     fontsize=6, ncol=ncol, frameon=False)

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'umap_mofa_rna_celltype.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: umap_mofa_rna_celltype.png")
    else:
        print("\nSKIPPED: Cell type UMAP — no cell_type column in mdata.obs "
              "(cell_type_source was not resolved during metadata extraction)")

    # 7.1b UMAP by ATAC cell type (FIX-67)
    atac_ct_key = None
    for candidate in ['atac:cell_type', 'atac:scanvi_prediction']:
        if candidate in mdata.obs.columns:
            atac_ct_key = candidate
            break
    if atac_ct_key is None and 'atac' in mdata.mod:
        if 'cell_type' in mdata.mod['atac'].obs.columns:
            mdata.obs['atac_cell_type'] = mdata.mod['atac'].obs['cell_type']
            atac_ct_key = 'atac_cell_type'

    if atac_ct_key:
        print(f"\nPlotting UMAP by ATAC cell type ({atac_ct_key})...")
        fig, ax = plt.subplots(figsize=(10, 8))
        cell_types = mdata.obs[atac_ct_key].astype(str)
        unique_types = cell_types.unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_types)))
        color_map = dict(zip(unique_types, colors))
        for ct_label in unique_types:
            mask = cell_types == ct_label
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                      c=[color_map[ct_label]], label=ct_label, s=5, alpha=0.7)
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title('UMAP colored by ATAC Cell Type')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'umap_mofa_atac_celltype.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: umap_mofa_atac_celltype.png")

    # 7.2 UMAP by clusters
    print("\nPlotting UMAP by MOFA clusters...")
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    clusters = mdata.obs['leiden_mofa'].astype(str)
    unique_clusters = sorted(clusters.unique(), key=lambda x: int(x) if x.isdigit() else 0)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
    color_map = dict(zip(unique_clusters, colors))
    
    for cluster in unique_clusters:
        mask = clusters == cluster
        ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                  c=[color_map[cluster]], label=cluster, s=5, alpha=0.7)
    
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP colored by MOFA Leiden clusters')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'umap_mofa_clusters.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Saved: umap_mofa_clusters.png")
    
    # 7.3 UMAP by batch
    if 'batch' in mdata.obs.columns:
        print(f"\nPlotting UMAP by batch...")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        batches = mdata.obs['batch'].astype(str)
        unique_batches = batches.unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(unique_batches)))
        color_map = dict(zip(unique_batches, colors))
        
        for batch in unique_batches:
            mask = batches == batch
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                      c=[color_map[batch]], label=batch, s=5, alpha=0.7)
        
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title('UMAP colored by Batch')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'umap_mofa_batch.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: umap_mofa_batch.png")
    
    # 7.4 UMAP by lane
    if 'lane' in mdata.obs.columns:
        print(f"\nPlotting UMAP by lane...")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        lanes = mdata.obs['lane'].astype(str)
        unique_lanes = lanes.unique()
        colors = plt.cm.tab20b(np.linspace(0, 1, len(unique_lanes)))
        color_map = dict(zip(unique_lanes, colors))
        
        for lane in unique_lanes:
            mask = lanes == lane
            ax.scatter(umap_coords[mask, 0], umap_coords[mask, 1],
                      c=[color_map[lane]], label=lane, s=5, alpha=0.7)
        
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title('UMAP colored by Lane')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=6)
        
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'umap_mofa_lane.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: umap_mofa_lane.png")
    
    # ==========================================================================
    # STEP 8: CELL TYPE - FACTOR CORRELATION ANALYSIS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 8: CELL TYPE - FACTOR ANALYSIS")
    print("="*80)
    
    if 'cell_type' in mdata.obs.columns:
        # Create factor-celltype DataFrame
        factor_celltype_df = pd.DataFrame({
            'Cell_Type': mdata.obs['cell_type'].values,
            **{f'Factor_{i+1}': mdata.obsm['X_mofa'][:, i] for i in range(n_factors)}
        })
        
        # Compute mean factor values per cell type
        factor_means = factor_celltype_df.groupby('Cell_Type').mean()
        
        print("\nCell types with highest loading per factor:")
        for col in factor_means.columns:
            max_celltype = factor_means[col].idxmax()
            max_val = factor_means[col].max()
            print(f"  {col}: {max_celltype} ({max_val:.3f})")
        
        # Save table
        factor_celltype_csv = OUTPUT_DIR / 'mofa_factors_by_celltype.csv'
        factor_means.to_csv(factor_celltype_csv)
        print(f"\n✓ Saved: {factor_celltype_csv}")
        
        # 8.1 Heatmap (dynamic sizing, annotations only when legible)
        print("\nGenerating factor-celltype heatmap...")
        n_types = len(factor_means.index)
        show_annot = n_types <= 15
        fig_w = max(10, n_types * 0.6)
        fig_h = max(8, n_factors * 0.5)
        plt.figure(figsize=(fig_w, fig_h))
        sns.heatmap(factor_means.T, cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'Mean Factor Value'},
                    annot=show_annot, fmt='.2f' if show_annot else '')
        plt.xlabel('Cell Type')
        plt.ylabel('MOFA Factor')
        plt.title('Mean MOFA Factor Values by Cell Type')
        plt.xticks(rotation=45, ha='right', fontsize=max(6, 10 - n_types // 10))
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'mofa_factors_celltype_heatmap.png',
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: mofa_factors_celltype_heatmap.png")
        
        # 8.2 Boxplots
        print("\nGenerating factor-celltype boxplots...")
        fig, axes = plt.subplots(3, 3, figsize=(15, 12))
        axes = axes.flatten()
        
        for i in range(min(9, n_factors)):
            factor_celltype_df.boxplot(column=f'Factor_{i+1}', by='Cell_Type', 
                                       ax=axes[i], rot=45)
            axes[i].set_title(f'Factor {i+1}')
            axes[i].set_xlabel('')
            axes[i].get_figure().suptitle('')
        
        # Hide unused subplots
        for i in range(min(9, n_factors), 9):
            axes[i].axis('off')
        
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'mofa_factors_celltype_boxplots.png',
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: mofa_factors_celltype_boxplots.png")
        
        # 8.3 ANOVA testing
        n_cell_types = factor_celltype_df['Cell_Type'].nunique()
        if n_cell_types < 2:
            print(f"\nSkipping ANOVA — only {n_cell_types} cell type group(s) (need >= 2)")
            anova_results = []
        else:
            print(f"\nPerforming ANOVA across {n_cell_types} cell types...")
            anova_results = []

            for i in range(n_factors):
                factor_name = f'Factor_{i+1}'

                groups = [factor_celltype_df[factor_celltype_df['Cell_Type'] == ct][factor_name].values
                          for ct in factor_celltype_df['Cell_Type'].unique()]

                f_stat, p_val = f_oneway(*groups)

                anova_results.append({
                    'Factor': factor_name,
                    'F_statistic': f_stat,
                    'P_value': p_val,
                    'Significant': 'Yes' if p_val < 0.05 else 'No'
                })

                sig_marker = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
                print(f"  {factor_name}: F={f_stat:.2f}, p={p_val:.2e} {sig_marker}")

        anova_df = pd.DataFrame(anova_results)
        anova_csv = OUTPUT_DIR / 'mofa_factors_anova_celltype.csv'
        anova_df.to_csv(anova_csv, index=False)
        print(f"\n✓ Saved: {anova_csv}")
    
    # ==========================================================================
    # STEP 9: BATCH EFFECTS IN FACTORS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 9: BATCH EFFECTS IN FACTORS")
    print("="*80)
    
    if 'batch' in mdata.obs.columns:
        print("\nGenerating factors-by-batch boxplots...")
        
        n_plot_factors = min(8, n_factors)
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        axes = axes.flatten()
        
        for i in range(n_plot_factors):
            plot_df = pd.DataFrame({
                'batch': mdata.obs['batch'].values,
                f'Factor_{i+1}': mdata.obsm['X_mofa'][:, i]
            })
            
            sns.boxplot(data=plot_df, x='batch', y=f'Factor_{i+1}', ax=axes[i])
            axes[i].set_title(f'Factor {i+1} by Batch')
            axes[i].set_xlabel('Batch')
            axes[i].set_ylabel(f'Factor {i+1} Value')
            axes[i].tick_params(axis='x', rotation=45)
        
        # Hide unused subplot
        for i in range(n_plot_factors, 8):
            axes[i].axis('off')
        
        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'mofa_factors_by_batch.png', 
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("✓ Saved: mofa_factors_by_batch.png")
    
    # ==========================================================================
    # STEP 10: SAVE INTEGRATED DATA AND OUTPUTS
    # ==========================================================================
    print("\n" + "="*80)
    print("STEP 10: SAVING OUTPUTS")
    print("="*80)
    
    # Generate timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Save integrated MuData
    output_file = OUTPUT_DIR / f"integrated_{timestamp}_mofa_complete.h5mu"
    sanitize_mudata(mdata, output_file)
    mdata.write_h5mu(output_file)
    print(f"\n✓ Saved integrated MuData: {output_file}")
    
    # Export factor loadings
    factors_export_df = pd.DataFrame(
        factors_array,
        index=mdata.obs_names,
        columns=[f"Factor_{i+1}" for i in range(n_factors)]
    )
    factors_export_df['leiden_mofa'] = mdata.obs['leiden_mofa'].values
    
    if 'batch' in mdata.obs.columns:
        factors_export_df['batch'] = mdata.obs['batch'].values
    if 'cell_type' in mdata.obs.columns:
        factors_export_df['cell_type'] = mdata.obs['cell_type'].values
    if 'lane' in mdata.obs.columns:
        factors_export_df['lane'] = mdata.obs['lane'].values
    
    factors_csv = OUTPUT_DIR / "mofa_factors.csv"
    factors_export_df.to_csv(factors_csv)
    print(f"✓ Saved factor loadings: {factors_csv}")
    
    # Save summary
    summary = {
        'timestamp': datetime.now().isoformat(),
        'n_cells': int(mdata.n_obs),
        'n_factors': int(n_factors),
        'n_clusters': int(mdata.obs['leiden_mofa'].nunique()),
        'variance_explained': {
            view: float(total) for view, total in view_totals.items()
        },
        'metadata_columns_extracted': {
            'cell_type': 'cell_type' in mdata.obs.columns,
            'batch': 'batch' in mdata.obs.columns,
            'lane': 'lane' in mdata.obs.columns
        },
        'output_file': str(output_file.name),
        'factors_csv': 'mofa_factors.csv',
        'variance_csv': 'mofa_variance_explained.csv'
    }
    
    summary_file = OUTPUT_DIR / "mofa_integration_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✓ Saved summary: {summary_file}")
    
    # ==========================================================================
    # FINAL SUMMARY
    # ==========================================================================
    print("\n" + "="*80)
    print("✅ MOFA+ VISUALIZATION COMPLETE!")
    print("="*80)
    
    print(f"\nSummary:")
    print(f"  {mdata.n_obs:,} cells analyzed")
    print(f"  {n_factors} MOFA factors")
    print(f"  {mdata.obs['leiden_mofa'].nunique()} clusters")
    
    print(f"\nVariance explained:")
    for view, var in summary['variance_explained'].items():
        view_name = 'RNA' if 'view0' in view else 'ATAC'
        print(f"  {view_name}: {var:.2f}%")
    
    # List all generated plots
    plot_files = list(FIGURES_DIR.glob('*.png'))
    print(f"\n✓ Figures generated: {len(plot_files)}")
    for plot_file in sorted(plot_files):
        print(f"  - {plot_file.name}")
    
    print(f"\nOutputs saved to: {OUTPUT_DIR}")
    print(f"  - Integrated MuData: {output_file.name}")
    print(f"  - Factor loadings: {factors_csv.name}")
    
    print("\n" + "="*80)
    
    model.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
