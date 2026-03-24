#!/usr/bin/env python3
"""
Enhanced concatenation with comprehensive QC metrics recalculation
Aligns with template workflow including HVG selection and initial embeddings
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import argparse
import json
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re

def concatenate_and_qc(h5ad_files, output_file, metrics_file, plots_dir):
    """Concatenate multiple h5ad files with comprehensive QC at experiment level"""
    
    # ============================================
    # SETUP AND FILE LOADING
    # ============================================
    
    # Create plots directory
    os.makedirs(plots_dir, exist_ok=True)
    
    # Load all samples
    adatas = {}
    sample_info = []
    
    for file in h5ad_files:
        print(f"Loading {file}...")
        adata = sc.read_h5ad(file)
        
        # Ensure counts layer exists
        if 'counts' not in adata.layers:
            print(f"  WARNING: No counts layer in {file}, creating from X")
            adata.layers['counts'] = adata.X.copy()
        
        # Extract sample name
        base_name = os.path.basename(file).replace('_filtered.h5ad', '').replace('_qc_filtered.h5ad', '')
        base_name = re.sub(r'^_\d+_', '', base_name)
        # Store in dictionary for named concatenation
        adatas[base_name] = adata
        
        # Collect sample info
        sample_info.append({
            'sample': base_name,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'batch': (
                'june' if 'june' in base_name.lower()
                else 'july' if 'july' in base_name.lower()
                else 'nov'  if 'nov'  in base_name.lower()
                else 'unknown'
            )
        })
    
    # ============================================
    # CONCATENATION
    # ============================================
    
    print("\nStarting concatenation of adata objects...")
    adata = ad.concat(adatas, join='outer', label="sample", index_unique=None)
    print("Concatenation done.")
    
    # Make names unique
    adata.obs_names_make_unique()
    print("Observation names made unique.")
    adata.var_names_make_unique()
    print("Variable names made unique.")
    
    # Print sample distribution
    print("\nSummary of observations per sample:")
    sample_counts = adata.obs["sample"].value_counts()
    n_samples = len(sample_counts)
    print(sample_counts)
    
    # Add batch information
    adata.obs['batch'] = adata.obs['sample'].apply(
        lambda x: (
            'june' if 'june' in x.lower()
            else 'july' if 'july' in x.lower()
            else 'nov'  if 'nov'  in x.lower()
            else 'unknown'
        )
    )

    # ============================================
    # EXPERIMENT-LEVEL QC (following template)
    # ============================================
    
    print("\nApplying experiment-level QC metrics...")
    
    # Define gene categories based on species
    # Assuming human based on your data
    adata.var["mt"] = adata.var_names.str.startswith("MT-")  # MT for human
    print("Identified mitochondrial genes.")
    
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))  # uppercase for human
    print("Identified ribosomal genes.")
    
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^P]")
    print("Identified hemoglobin genes.")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    
    # Preserve counts layer
    adata.layers["counts"] = adata.X.copy()
    
    # ============================================
    # NORMALIZATION AND HVG SELECTION
    # ============================================
    
    print("\nStarting normalization and log transformation...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Store raw counts
    adata.raw = adata
    print("Experiment-level normalization and log transformation done.")
    
    print("\nSelecting highly variable genes...")
    # Use batch_key to account for batch effects
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=2000,
            batch_key='batch',
            subset=False,
            flavor='seurat_v3'
        )
    except (ValueError, IndexError) as e:
        print(f"Warning: Batch-wise HVG selection failed ({e})")
        print("Falling back to non-batch HVG selection...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=2000,
            subset=False,
            flavor='seurat_v3'
        )
    
    # Plot HVG
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig(os.path.join(plots_dir, "highly_variable_genes.png"), dpi=300)
    plt.close()
    
    # ============================================
    # INITIAL DIMENSIONALITY REDUCTION
    # ============================================
    
    print("\nRunning PCA, neighbor calculation, and UMAP...")
    
    # PCA on HVGs
    sc.tl.pca(adata, n_comps=50, mask_var="highly_variable")
    
    # Plot PCA variance
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig(os.path.join(plots_dir, "pca_variance_ratio.png"), dpi=300)
    plt.close()
    
    # Compute neighbors and UMAP
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    # Initial clustering at multiple resolutions
    resolutions = [0.5, 5.0]
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}',
                 flavor='igraph', n_iterations=-1)

    # Save cluster counts for a coarse resolution (e.g., 0.5)
    cluster_key = 'leiden_0.5'
    if cluster_key in adata.obs.columns:
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        cluster_counts.to_csv(
            os.path.join(plots_dir, f'{cluster_key}_cluster_counts.csv'),
            header=['cell_count']
        )
        print(f"Cluster counts saved to {cluster_key}_cluster_counts.csv")
    
    # ============================================
    # QC VISUALIZATIONS
    # ============================================
    
    print("\nGenerating QC visualizations...")
    
    # UMAP colored by sample
    n_samples = len(sample_counts)

    fig_width = max(10, 0.12 * n_samples)  # wider for many samples
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    sc.pl.umap(
        adata,
        color="sample",
        size=2,
        show=False,
        ax=ax,
        legend_loc='right margin',   # use the extra width
        legend_fontsize=5,
    )

    plt.savefig(os.path.join(plots_dir, "umap_by_sample.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # UMAP colored by batch
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color="batch", show=False, ax=ax)
    plt.savefig(os.path.join(plots_dir, "umap_by_batch.png"), dpi=300)
    plt.close()
    
    # QC violin plots
    n_samples = len(sample_counts)
    fig_width = max(18, 0.25 * n_samples)  # 0.25–0.3 inches per sample works well
    fig_height = 10

    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_height))

    sc.pl.violin(
        adata, 'pct_counts_mt', groupby='sample',
        rotation=90, ax=axes[0, 0], show=False
    )
    axes[0, 0].set_title('MT% by Sample')

    sc.pl.violin(
        adata, 'total_counts', groupby='sample',
        rotation=90, ax=axes[0, 1], show=False
    )
    axes[0, 1].set_title('Total Counts by Sample')

    sc.pl.violin(
        adata, 'n_genes_by_counts', groupby='sample',
        rotation=90, ax=axes[1, 0], show=False
    )
    axes[1, 0].set_title('Genes by Counts by Sample')

    sc.pl.violin(
        adata, 'pct_counts_ribo', groupby='sample',
        rotation=90, ax=axes[1, 1], show=False
    )
    axes[1, 1].set_title('Ribosomal% by Sample')

    for ax in axes.flat:
        ax.tick_params(axis='x', labelsize=6)  # shrink font a bit

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "qc_metrics_by_sample.png"), dpi=300)
    plt.close()
    
    # ============================================
    # PER-SAMPLE METRICS (ROW FORMAT)
    # ============================================
    
    print("\nGenerating comprehensive per-sample metrics...")
    
    # Create per-sample summary with each sample as a row
    sample_metrics = []
    for sample_name in sample_counts.index:
        sample_mask = adata.obs['sample'] == sample_name
        sample_adata = adata[sample_mask]
        
        metrics = {
            'sample': sample_name,
            'cells': sample_adata.n_obs,
            'genes': sample_adata.n_vars,
            'median_genes_per_cell': sample_adata.obs['n_genes_by_counts'].median(),
            'median_umi_per_cell': sample_adata.obs['total_counts'].median(),
            'median_mt_percent': sample_adata.obs['pct_counts_mt'].median(),
            'median_ribo_percent': sample_adata.obs['pct_counts_ribo'].median(),
            'median_hb_percent': sample_adata.obs['pct_counts_hb'].median(),
            'mean_genes_per_cell': sample_adata.obs['n_genes_by_counts'].mean(),
            'mean_umi_per_cell': sample_adata.obs['total_counts'].mean(),
            'mean_mt_percent': sample_adata.obs['pct_counts_mt'].mean(),
            'mean_ribo_percent': sample_adata.obs['pct_counts_ribo'].mean(),
            'pct_of_total_cells': (sample_adata.n_obs / adata.n_obs) * 100,
            'batch': adata.obs.loc[sample_mask, 'batch'].iloc[0]
        }
        sample_metrics.append(metrics)
    
    # Convert to DataFrame (rows = samples)
    sample_df = pd.DataFrame(sample_metrics)
    
    # Add overall summary row
    overall_metrics = {
        'sample': 'OVERALL_TOTAL',
        'cells': adata.n_obs,
        'genes': adata.n_vars,
        'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
        'median_umi_per_cell': adata.obs['total_counts'].median(),
        'median_mt_percent': adata.obs['pct_counts_mt'].median(),
        'median_ribo_percent': adata.obs['pct_counts_ribo'].median(),
        'median_hb_percent': adata.obs['pct_counts_hb'].median(),
        'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
        'mean_umi_per_cell': adata.obs['total_counts'].mean(),
        'mean_mt_percent': adata.obs['pct_counts_mt'].mean(),
        'mean_ribo_percent': adata.obs['pct_counts_ribo'].mean(),
        'pct_of_total_cells': 100.0,
        'batch': 'all'
    }
    sample_df = pd.concat([sample_df, pd.DataFrame([overall_metrics])], ignore_index=True)
    
    # Save per-sample metrics (PRIMARY OUTPUT - ROW FORMAT)
    sample_df.to_csv(metrics_file, index=False)
    print(f"Per-sample metrics (row format) saved to: {metrics_file}")

    # ------------------------------------------------------------------
    # Neural marker dotplot (template-style, pre-integration)
    # ------------------------------------------------------------------
    key_markers = {
        'Oligodendrocyte': ['MBP', 'PLP1', 'MAG', 'MOG'],
        'Astrocyte': ['GFAP', 'S100B', 'ALDH1L1', 'SLC1A3'],
        'Microglia': ['P2RY12', 'CX3CR1', 'TMEM119'],
        'Vasculature': ['PDGFRB', 'FLT1'],
        'Pan-neuronal': ['SYT1'],
        'Excitatory Neuron': ['SLC17A7', 'SATB2', 'RORB', 'FEZF2', 'THEMIS'],
        'Inhibitory Neuron': ['GAD1', 'GAD2', 'PVALB', 'SST', 'VIP', 'LAMP5'],
        'OPC': ['PDGFRA', 'CSPG4', 'OLIG2'],
        'Endothelial cells': ['PECAM1', 'VWF', 'VEGFA', 'CLDN5'],
        'Pericytes': ['ACTA2', 'S100A4', 'NRP1'],
        'Immune cells': ['CD4', 'CD68'],
    }
    # Restrict to markers actually present in the data
    markers_in_data = {}
    for cell_type, markers in key_markers.items():
        present = [m for m in markers if m in adata.var_names]
        if present:
            markers_in_data[cell_type] = present

    if markers_in_data and 'leiden_0.5' in adata.obs.columns:
        try:
            dp = sc.pl.dotplot(
                adata,
                markers_in_data,
                groupby="leiden_0.5",
                standard_scale="var",
                swap_axes=True,
                show=False,
            )

            fig = plt.gcf()
            fig.savefig(
                os.path.join(plots_dir, "neural_markers_dotplot_leiden_0_5.png"),
                dpi=300,
                bbox_inches="tight",
            )
            plt.close(fig)
            print("Neural marker dotplot saved.")

        except Exception as e:
            print(f"Warning: could not generate neural marker dotplot: {e}")
    
    # ============================================
    # CONSOLIDATE LANE-LEVEL QC METRICS
    # ============================================
    
    print("\nConsolidating lane-level QC metrics from RNA_QC step...")
    
    # Find all QC metrics CSVs
    qc_dir = os.path.dirname(os.path.dirname(metrics_file))  # Go up to results/rna/
    qc_metrics_dir = os.path.join(qc_dir, 'qc')
    
    all_qc_dfs = []
    if os.path.exists(qc_metrics_dir):
        qc_files = glob.glob(os.path.join(qc_metrics_dir, '*_qc_metrics.csv'))
        print(f"Found {len(qc_files)} QC metrics files to consolidate...")
        
        for qc_file in qc_files:
            try:
                df = pd.read_csv(qc_file)
                all_qc_dfs.append(df)
            except Exception as e:
                print(f"  Warning: Could not read {qc_file}: {e}")
    
    if all_qc_dfs:
        # Concatenate all lane-level QC metrics
        consolidated_qc = pd.concat(all_qc_dfs, ignore_index=True)
        
        # Save consolidated QC metrics
        consolidated_file = metrics_file.replace('.csv', '_consolidated_qc.csv')
        consolidated_qc.to_csv(consolidated_file, index=False)
        print(f"Consolidated QC metrics saved to: {consolidated_file}")
        print(f"  Total rows (lanes + donors): {len(consolidated_qc)}")

        # ------------------------------------------------------------------
        # Template-style plots: filtered vs unfiltered, diffs, heatmap, stacked bar
        # Use ONLY lane-level rows (where unfiltered_cells is present)
        # ------------------------------------------------------------------
        qc = consolidated_qc.copy()
        qc = qc[pd.to_numeric(qc['unfiltered_cells'], errors='coerce').notna()].copy()
        qc['unfiltered_cells'] = pd.to_numeric(qc['unfiltered_cells'], errors='coerce')
        qc['unfiltered_genes'] = pd.to_numeric(qc['unfiltered_genes'], errors='coerce')
        qc['final_filtered_total'] = pd.to_numeric(qc['final_filtered_total'], errors='coerce')
        qc['final_filtered_total_genes'] = pd.to_numeric(qc['final_filtered_total_genes'], errors='coerce')

        # Differences
        qc['Diff_Cells'] = qc['unfiltered_cells'] - qc['final_filtered_total']
        qc['Diff_Genes'] = qc['unfiltered_genes'] - qc['final_filtered_total_genes']

        n_lanes = len(qc)

        # 4-panel plot: filtered vs unfiltered + differences
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))

        # Plot 1: Filtered vs Unfiltered Cells
        axes[0, 0].plot(qc['sample'], qc['unfiltered_cells'], label='Unfiltered Cells', marker='o')
        axes[0, 0].plot(qc['sample'], qc['final_filtered_total'], label='Filtered Cells', marker='x')
        axes[0, 0].set_xlabel('Sample (lane)')
        axes[0, 0].set_ylabel('Number of Cells')
        axes[0, 0].set_title('Filtered vs Unfiltered Cells')
        axes[0, 0].legend()
        axes[0, 0].tick_params(axis='x', rotation=90)

        avg_unfiltered_cells = qc['unfiltered_cells'].mean()
        sem_unfiltered_cells = qc['unfiltered_cells'].std() / np.sqrt(n_lanes)
        avg_filtered_cells = qc['final_filtered_total'].mean()
        sem_filtered_cells = qc['final_filtered_total'].std() / np.sqrt(n_lanes)
        axes[0, 0].text(
            0.95, 0.95,
            f'Avg Unfiltered: {avg_unfiltered_cells:.2f} ± {sem_unfiltered_cells:.2f}\n'
            f'Avg Filtered: {avg_filtered_cells:.2f} ± {sem_filtered_cells:.2f}',
            ha='right', va='top', transform=axes[0, 0].transAxes,
            bbox=dict(facecolor='white', alpha=0.7)
        )

        # Plot 2: Filtered vs Unfiltered Genes
        axes[0, 1].plot(qc['sample'], qc['unfiltered_genes'], label='Unfiltered Genes', marker='o')
        axes[0, 1].plot(qc['sample'], qc['final_filtered_total_genes'], label='Filtered Genes', marker='x')
        axes[0, 1].set_xlabel('Sample (lane)')
        axes[0, 1].set_ylabel('Number of Genes')
        axes[0, 1].set_title('Filtered vs Unfiltered Genes')
        axes[0, 1].legend()
        axes[0, 1].tick_params(axis='x', rotation=90)

        avg_unfiltered_genes = qc['unfiltered_genes'].mean()
        sem_unfiltered_genes = qc['unfiltered_genes'].std() / np.sqrt(n_lanes)
        avg_filtered_genes = qc['final_filtered_total_genes'].mean()
        sem_filtered_genes = qc['final_filtered_total_genes'].std() / np.sqrt(n_lanes)
        axes[0, 1].text(
            0.95, 0.95,
            f'Avg Unfiltered: {avg_unfiltered_genes:.2f} ± {sem_unfiltered_genes:.2f}\n'
            f'Avg Filtered: {avg_filtered_genes:.2f} ± {sem_filtered_genes:.2f}',
            ha='right', va='top', transform=axes[0, 1].transAxes,
            bbox=dict(facecolor='white', alpha=0.7)
        )

        # Plot 3: Difference in Cells
        axes[1, 0].bar(qc['sample'], qc['Diff_Cells'], color='blue')
        axes[1, 0].set_xlabel('Sample (lane)')
        axes[1, 0].set_ylabel('Difference in Number of Cells')
        axes[1, 0].set_title('Difference (Unfiltered - Filtered) in Cells')
        axes[1, 0].tick_params(axis='x', rotation=90)
        avg_diff_cells = qc['Diff_Cells'].mean()
        sem_diff_cells = qc['Diff_Cells'].std() / np.sqrt(n_lanes)
        axes[1, 0].text(
            0.95, 0.95,
            f'Avg Diff: {avg_diff_cells:.2f}\nSEM: {sem_diff_cells:.2f}',
            ha='right', va='top', transform=axes[1, 0].transAxes,
            bbox=dict(facecolor='white', alpha=0.7)
        )

        # Plot 4: Difference in Genes
        axes[1, 1].bar(qc['sample'], qc['Diff_Genes'], color='green')
        axes[1, 1].set_xlabel('Sample (lane)')
        axes[1, 1].set_ylabel('Difference in Number of Genes')
        axes[1, 1].set_title('Difference (Unfiltered - Filtered) in Genes')
        axes[1, 1].tick_params(axis='x', rotation=90)
        avg_diff_genes = qc['Diff_Genes'].mean()
        sem_diff_genes = qc['Diff_Genes'].std() / np.sqrt(n_lanes)
        axes[1, 1].text(
            0.95, 0.95,
            f'Avg Diff: {avg_diff_genes:.2f}\nSEM: {sem_diff_genes:.2f}',
            ha='right', va='top', transform=axes[1, 1].transAxes,
            bbox=dict(facecolor='white', alpha=0.7)
        )

        fig.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'pre_post_filter_and_differences_cells_genes.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # Filtering efficiency heatmap
        filtering_steps = [
            'frac_filt_mingenes', 'frac_filt_3cellthresh', 'cellfracfiltScr',
            'ngbc_frac_filt', 'pctmt_frac_filt', 'tc_frac_filt'
        ]
        # keep only existing columns
        filtering_steps = [c for c in filtering_steps if c in qc.columns]
        if filtering_steps:
            heatmap_data = qc.set_index('sample')[filtering_steps]
            plt.figure(figsize=(12, 8))
            sns.heatmap(heatmap_data, annot=True, cmap='viridis',
                        cbar_kws={'label': 'Fraction Filtered'})
            plt.title('Filtering Efficiency Heatmap')
            plt.xlabel('Filtering Step')
            plt.ylabel('Sample (lane)')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'filtering_efficiency_heatmap.png'),
                        dpi=300, bbox_inches='tight')
            plt.close()

        # Stacked bar plot of cells remaining after each filtering step
        stacked_cols = [
            'post_mingenes_cells', 'post_3cellthresh_cells', 'total_cells_afterScr',
            'post_ngbc_filter_cells', 'post_pctcountsmt_filter_cells',
            'post_totalcounts_filter_cells'
        ]
        stacked_cols = [c for c in stacked_cols if c in qc.columns]
        if stacked_cols:
            stacked_data = qc.set_index('sample')[stacked_cols]
            stacked_data.plot(kind='bar', stacked=True, figsize=(14, 8), colormap='tab20')
            plt.title('Cells Remaining After Each Filtering Step (Lane Level)')
            plt.xlabel('Sample (lane)')
            plt.ylabel('Number of Cells')
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Filtering Step', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'cells_remaining_after_each_filter_step_stacked.png'),
                        dpi=300, bbox_inches='tight')
            plt.close()
    else:
        print("  No QC metrics files found to consolidate")
    
    # ============================================
    # ADDITIONAL COMPREHENSIVE PLOTS
    # ============================================
    
    print("\nGenerating additional comprehensive plots...")
    
    # Plot 1: Cells per sample bar chart
    fig, ax = plt.subplots(figsize=(max(14, 0.3 * n_samples), 6))
    sample_df_plot = sample_df[sample_df['sample'] != 'OVERALL_TOTAL'].copy()
    sample_df_plot = sample_df_plot.sort_values('cells', ascending=False)
    
    bars = ax.bar(range(len(sample_df_plot)), sample_df_plot['cells'])
    ax.set_xticks(range(len(sample_df_plot)))
    ax.set_xticklabels(sample_df_plot['sample'], rotation=45, ha='right', fontsize=6)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Number of Cells')
    ax.set_title('Cells per Sample After QC Filtering')
    ax.axhline(y=sample_df_plot['cells'].mean(), color='r', linestyle='--', 
               label=f"Mean: {sample_df_plot['cells'].mean():.0f}")
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'cells_per_sample_barplot.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: QC metrics distributions (4-panel)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Genes per cell
    axes[0, 0].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
    median_genes = adata.obs['n_genes_by_counts'].median()
    axes[0, 0].axvline(median_genes, color='r', linestyle='--', 
                   label=f"Median: {median_genes:.0f}")
    axes[0, 0].set_xlabel('Genes per Cell')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Genes per Cell')
    axes[0, 0].legend()

    # UMI per cell
    axes[0, 1].hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7)
    median_umi = adata.obs['total_counts'].median()
    axes[0, 1].axvline(median_umi, color='r', linestyle='--', 
                   label=f"Median: {median_umi:.0f}")
    axes[0, 1].set_xlabel('UMI per Cell')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Distribution of UMI per Cell')
    axes[0, 1].legend()

    # MT percent
    axes[1, 0].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7)
    median_mt = adata.obs['pct_counts_mt'].median()
    axes[1, 0].axvline(median_mt, color='r', linestyle='--', 
                   label=f"Median: {median_mt:.2f}%")
    axes[1, 0].set_xlabel('Mitochondrial %')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Mitochondrial %')
    axes[1, 0].legend()

    # ALL samples by cell count (sorted, horizontal bar)
    all_samples = sample_df_plot.sort_values('cells', ascending=True)  # ascending for bottom-to-top
    fig_height_samples = max(12, len(all_samples) * 0.3)  # 0.3 inches per sample

    # Create subplot with dynamic height
    plt.close(fig)  # Close the 2x2 subplot
    fig = plt.figure(figsize=(15, fig_height_samples))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])

    # Recreate first 3 plots
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
    ax1.axvline(median_genes, color='r', linestyle='--', label=f"Median: {median_genes:.0f}")
    ax1.set_xlabel('Genes per Cell')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Genes per Cell')
    ax1.legend()

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7)
    ax2.axvline(median_umi, color='r', linestyle='--', label=f"Median: {median_umi:.0f}")
    ax2.set_xlabel('UMI per Cell')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of UMI per Cell')
    ax2.legend()

    ax3 = fig.add_subplot(gs[1, 0])
    ax3.hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7)
    ax3.axvline(median_mt, color='r', linestyle='--', label=f"Median: {median_mt:.2f}%")
    ax3.set_xlabel('Mitochondrial %')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Mitochondrial %')
    ax3.legend()

    # ALL samples horizontal bar chart (bottom-right, spanning full height)
    ax4 = fig.add_subplot(gs[1, 1])
    y_pos = np.arange(len(all_samples))
    ax4.barh(y_pos, all_samples['cells'], color='steelblue', alpha=0.7)
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(all_samples['sample'], fontsize=6)
    ax4.set_xlabel('Number of Cells')
    ax4.set_title(f'All {len(all_samples)} Samples by Cell Count')
    ax4.invert_yaxis()  # Highest at top

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'qc_distributions_4panel.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()    
    
    # Plot 3: Batch comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    batch_summary = sample_df_plot.groupby('batch').agg({
        'cells': 'sum',
        'median_genes_per_cell': 'mean',
        'median_umi_per_cell': 'mean'
    }).reset_index()
    
    axes[0].bar(batch_summary['batch'], batch_summary['cells'])
    axes[0].set_ylabel('Total Cells')
    axes[0].set_title('Cells per Batch')
    axes[0].tick_params(axis='x', rotation=45)
    
    axes[1].bar(batch_summary['batch'], batch_summary['median_genes_per_cell'])
    axes[1].set_ylabel('Mean Median Genes/Cell')
    axes[1].set_title('Genes per Cell by Batch')
    axes[1].tick_params(axis='x', rotation=45)
    
    axes[2].bar(batch_summary['batch'], batch_summary['median_umi_per_cell'])
    axes[2].set_ylabel('Mean Median UMI/Cell')
    axes[2].set_title('UMI per Cell by Batch')
    axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'batch_comparison.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"All plots saved to: {plots_dir}")
    
    # ============================================
    # SUMMARY METRICS FOR LOGGING
    # ============================================
    
    # Calculate comprehensive summary metrics
    concat_summary = {
        'total_cells': int(adata.n_obs),
        'total_genes': int(adata.n_vars),
        'n_samples': n_samples,
        'n_batches': int(adata.obs['batch'].nunique()),
        'n_hvgs': int(adata.var['highly_variable'].sum()),
        'median_genes_per_cell': float(adata.obs['n_genes_by_counts'].median()),
        'median_umi_per_cell': float(adata.obs['total_counts'].median()),
        'median_mt_percent': float(adata.obs['pct_counts_mt'].median()),
        'median_ribo_percent': float(adata.obs['pct_counts_ribo'].median()),
        'n_mt_genes': int(adata.var['mt'].sum()),
        'n_ribo_genes': int(adata.var['ribo'].sum()),
        'n_hb_genes': int(adata.var['hb'].sum())
    }
    
    # Save JSON summary for easy loading
    with open(metrics_file.replace('.csv', '_summary.json'), 'w') as f:
        json.dump(concat_summary, f, indent=2)
    
    # ============================================
    # SAVE OUTPUT
    # ============================================
    
    print(f"\nSaving concatenated data to {output_file}...")
    adata.write(output_file)
    
    print("\n" + "="*70)
    print("CONCATENATION COMPLETE - SUMMARY")
    print("="*70)
    print(f"Total cells:             {concat_summary['total_cells']:>10,}")
    print(f"Total genes:             {concat_summary['total_genes']:>10,}")
    print(f"Number of samples:       {concat_summary['n_samples']:>10}")
    print(f"Number of batches:       {concat_summary['n_batches']:>10}")
    print(f"Highly variable genes:   {concat_summary['n_hvgs']:>10,}")
    print(f"Median genes/cell:       {concat_summary['median_genes_per_cell']:>10.1f}")
    print(f"Median UMI/cell:         {concat_summary['median_umi_per_cell']:>10.1f}")
    print(f"Median MT%:              {concat_summary['median_mt_percent']:>10.2f}")
    print(f"Median Ribo%:            {concat_summary['median_ribo_percent']:>10.2f}")
    print("="*70)
    print(f"\nOutputs saved:")
    print(f"  - Concatenated h5ad:     {output_file}")
    print(f"  - Per-sample metrics:    {metrics_file}")
    print(f"  - Summary JSON:          {metrics_file.replace('.csv', '_summary.json')}")
    if all_qc_dfs:
        print(f"  - Consolidated QC:       {metrics_file.replace('.csv', '_consolidated_qc.csv')}")
    print(f"  - Plots directory:       {plots_dir}/")
    print("="*70)
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='Enhanced concatenation with QC')
    parser.add_argument('--inputs', nargs='+', required=True, 
                        help='Input h5ad files')
    parser.add_argument('--output', required=True, 
                        help='Output concatenated h5ad file')
    parser.add_argument('--metrics', required=True, 
                        help='Output metrics CSV file')
    parser.add_argument('--plots_dir', default='concat_plots',
                        help='Directory for QC plots')
    
    args = parser.parse_args()
    
    # Run enhanced concatenation
    concatenate_and_qc(args.inputs, args.output, args.metrics, args.plots_dir)

if __name__ == "__main__":
    main()
