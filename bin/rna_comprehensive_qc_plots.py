#!/usr/bin/env python3
"""
Comprehensive RNA QC and Post-Integration Plotting
Generates UMAPs, dotplots, and QC statistics at multiple resolutions
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from pathlib import Path

def create_qc_summary_table(adata_dict, output_dir):
    """Create comprehensive QC statistics table"""
    
    sample_summary_data = []
    
    for sample_name, adata in adata_dict.items():
        stats = {
            'sample': sample_name,
            'unfiltered_cells': adata.n_obs,
            'unfiltered_genes': adata.n_vars,
        }
        
        # Add filtering statistics if available
        if 'qc_stats' in adata.uns:
            stats.update(adata.uns['qc_stats'])
        
        # Calculate derived metrics
        if 'n_genes_by_counts' in adata.obs.columns:
            stats['pre_mingenes_cells'] = adata.n_obs
            stats['pre_mingenes_genes'] = adata.n_vars
            
            mask_mingenes = adata.obs['n_genes_by_counts'] >= 100
            stats['post_mingenes_cells'] = mask_mingenes.sum()
            stats['frac_filt_mingenes'] = 1 - (stats['post_mingenes_cells'] / stats['pre_mingenes_cells'])
        
        # MT percentage filtering
        if 'pct_counts_mt' in adata.obs.columns:
            mask_mt = adata.obs['pct_counts_mt'] < 10
            stats['post_pctcountsmt_filter_cells'] = mask_mt.sum()
            stats['pctmt_frac_filt'] = 1 - (stats['post_pctcountsmt_filter_cells'] / adata.n_obs)
        
        stats['final_filtered_total'] = adata.n_obs
        stats['final_filtered_total_genes'] = adata.n_vars
        
        sample_summary_data.append(stats)
    
    df = pd.DataFrame(sample_summary_data)
    df.to_csv(os.path.join(output_dir, 'sample_QC_descriptive_statistics.csv'), index=False)
    
    return df

def create_filtering_plots(qc_df, output_dir):
    """Create 4-panel filtering visualization"""
    
    qc_df['Diff_Cells'] = qc_df['unfiltered_cells'] - qc_df['final_filtered_total']
    qc_df['Diff_Genes'] = qc_df['unfiltered_genes'] - qc_df['final_filtered_total_genes']
    
    n_samples = len(qc_df)
    avg_diff_cells = qc_df['Diff_Cells'].mean()
    sem_diff_cells = qc_df['Diff_Cells'].std() / np.sqrt(n_samples)
    avg_diff_genes = qc_df['Diff_Genes'].mean()
    sem_diff_genes = qc_df['Diff_Genes'].std() / np.sqrt(n_samples)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Filtered vs Unfiltered Cells
    axes[0, 0].plot(qc_df['sample'], qc_df['unfiltered_cells'], 
                    label='Unfiltered Cells', marker='o')
    axes[0, 0].plot(qc_df['sample'], qc_df['final_filtered_total'], 
                    label='Filtered Cells', marker='x')
    axes[0, 0].set_xlabel('Sample Name')
    axes[0, 0].set_ylabel('Number of Cells')
    axes[0, 0].set_title('Filtered vs Unfiltered Cells')
    axes[0, 0].legend()
    axes[0, 0].tick_params(axis='x', rotation=90)
    
    avg_unfiltered_cells = qc_df['unfiltered_cells'].mean()
    sem_unfiltered_cells = qc_df['unfiltered_cells'].std() / np.sqrt(n_samples)
    avg_filtered_cells = qc_df['final_filtered_total'].mean()
    sem_filtered_cells = qc_df['final_filtered_total'].std() / np.sqrt(n_samples)
    axes[0, 0].text(0.95, 0.95, 
        f'Avg Unfiltered: {avg_unfiltered_cells:.2f} ± {sem_unfiltered_cells:.2f}\n'
        f'Avg Filtered: {avg_filtered_cells:.2f} ± {sem_filtered_cells:.2f}', 
        horizontalalignment='right', verticalalignment='top', 
        transform=axes[0, 0].transAxes, 
        bbox=dict(facecolor='white', alpha=0.7))
    
    # Plot 2: Filtered vs Unfiltered Genes
    axes[0, 1].plot(qc_df['sample'], qc_df['unfiltered_genes'], 
                    label='Unfiltered Genes', marker='o')
    axes[0, 1].plot(qc_df['sample'], qc_df['final_filtered_total_genes'], 
                    label='Filtered Genes', marker='x')
    axes[0, 1].set_xlabel('Sample Name')
    axes[0, 1].set_ylabel('Number of Genes')
    axes[0, 1].set_title('Filtered vs Unfiltered Genes')
    axes[0, 1].legend()
    axes[0, 1].tick_params(axis='x', rotation=90)
    
    # Plot 3: Difference in Cells
    axes[1, 0].bar(qc_df['sample'], qc_df['Diff_Cells'], color='blue')
    axes[1, 0].set_xlabel('Sample Name')
    axes[1, 0].set_ylabel('Difference in Number of Cells')
    axes[1, 0].set_title('Difference (Unfiltered - Filtered) in Cells')
    axes[1, 0].tick_params(axis='x', rotation=90)
    axes[1, 0].text(0.95, 0.95, 
                    f'Avg Diff: {avg_diff_cells:.2f}\nSEM: {sem_diff_cells:.2f}', 
                    horizontalalignment='right', verticalalignment='top', 
                    transform=axes[1, 0].transAxes, 
                    bbox=dict(facecolor='white', alpha=0.7))
    
    # Plot 4: Difference in Genes
    axes[1, 1].bar(qc_df['sample'], qc_df['Diff_Genes'], color='green')
    axes[1, 1].set_xlabel('Sample Name')
    axes[1, 1].set_ylabel('Difference in Number of Genes')
    axes[1, 1].set_title('Difference (Unfiltered - Filtered) in Genes')
    axes[1, 1].tick_params(axis='x', rotation=90)
    axes[1, 1].text(0.95, 0.95, 
                    f'Avg Diff: {avg_diff_genes:.2f}\nSEM: {sem_diff_genes:.2f}', 
                    horizontalalignment='right', verticalalignment='top', 
                    transform=axes[1, 1].transAxes, 
                    bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pre_post_filter_and_differences.png'), dpi=300)
    plt.close()

def create_leiden_umaps_dotplots(adata, output_dir, resolutions=[0.5, 5.0], n_hvg=2000):
    """
    Create UMAPs and dotplots at different Leiden resolutions
    Uses mouse brain marker genes
    """
    
    # Mouse brain marker genes (from your specification)
    key_markers = {
        'Oligodendrocyte': ['Mbp', 'Plp1', 'Mag', 'Mog'],
        'Astrocyte': ['Gfap', 'S100b', 'Aldh1l1', 'Slc1a3'],
        'Microglia': ['P2ry12', 'Cx3cr1', 'Tmem119'],
        'Vasculature': ['Pdgfrb', 'Flt1'],
        'Pan-neuronal': ['Syt1'],
        'Excitatory Neuron': ['Slc17a7', 'Satb2', 'Rorb', 'Themis', 'Fezf2'],
        'Inhibitory Neuron': ['Gad1', 'Gad2', 'Pvalb', 'Sst', 'Vip', 'Lamp5'],
        'Ependyma': ['Foxj1'],
        'OPC': ['Pdgfra', 'Cspg4', 'Olig2'],
        'Immune cells': ['Cd4', 'Cd68'],
        'Pericytes': ['Acta2', 'S100a4', 'Nrp1'],
        'Endothelial cells': ['Pecam1', 'Vwf', 'Vegfa', 'Cldn5']
    }
    
    # Identify HVG and create filtered marker dict
    print(f"Identifying {n_hvg} highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=False, batch_key="sample")
    sc.pl.highly_variable_genes(adata)
    
    highly_variable_genes = adata.var_names[adata.var['highly_variable']].tolist()
    
    # Create HVG_markers_in_data
    HVG_markers_in_data = {}
    for cell_type, markers in key_markers.items():
        present_markers = [marker for marker in markers if marker in highly_variable_genes]
        if present_markers:
            HVG_markers_in_data[cell_type] = present_markers
    
    print(f"Found {len(HVG_markers_in_data)} cell types with markers in HVG set")
    
    # PCA, neighbors, UMAP
    print("Running PCA, neighbor calculation, and UMAP...")
    sc.tl.pca(adata, n_comps=50, mask_var="highly_variable")
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    # Iterate through resolutions
    for res in resolutions:
        res_str = str(res).replace('.', '_')
        cluster_key = f'hvg{n_hvg//1000}k_igraph_res{res}'
        
        print(f"\nProcessing resolution {res}...")
        
        # Leiden clustering with igraph
        sc.tl.leiden(adata, key_added=cluster_key, resolution=res, 
                     n_iterations=-1, flavor='igraph')
        
        # UMAP colored by clusters
        sc.pl.umap(adata, color=[cluster_key], 
                   save=f'_{cluster_key}.pdf')
        
        # Dotplot with markers
        if HVG_markers_in_data:
            sc.pl.dotplot(adata, HVG_markers_in_data, groupby=cluster_key, 
                         standard_scale="var", swap_axes=True,
                         save=f'_{cluster_key}_markers.pdf')
        
        # Quantify cell counts per cluster
        cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
        summary_df = pd.DataFrame({
            "Cluster": cluster_counts.index.astype(str),
            "Cell Count": cluster_counts.values
        })
        
        summary_df.to_csv(
            os.path.join(output_dir, f'{cluster_key}_cluster_counts.csv'),
            index=False
        )
        
        print(f"  Resolution {res}: {len(cluster_counts)} clusters")
        print(summary_df)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Concatenated h5ad file')
    parser.add_argument('--output_dir', default='rna_qc_plots')
    parser.add_argument('--n_hvg', type=int, default=2000, 
                       help='Number of highly variable genes')
    parser.add_argument('--resolutions', nargs='+', type=float, 
                       default=[0.2, 0.5, 1.0, 5.0],
                       help='Leiden clustering resolutions')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    # Create sample-wise QC summary if samples are tracked
    if 'sample' in adata.obs.columns:
        print("Creating sample-wise QC summaries...")
        sample_dict = {}
        for sample in adata.obs['sample'].unique():
            sample_dict[sample] = adata[adata.obs['sample'] == sample].copy()
        
        qc_df = create_qc_summary_table(sample_dict, args.output_dir)
        create_filtering_plots(qc_df, args.output_dir)
    
    # Create Leiden UMAPs and dotplots
    print("\nGenerating clustering and visualization...")
    create_leiden_umaps_dotplots(adata, args.output_dir, 
                                 resolutions=args.resolutions,
                                 n_hvg=args.n_hvg)
    
    print(f"\n✅ All visualizations saved to {args.output_dir}")

if __name__ == "__main__":
    main()
