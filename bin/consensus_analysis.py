#!/usr/bin/env python3
"""
Consensus Analysis for Bootstrap MOFA Results
Compares factors, clusters, and integration quality across iterations
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import adjusted_rand_score
import argparse
import sys

def load_bootstrap_results(bootstrap_dir):
    """Load all bootstrap iteration results"""
    bootstrap_dir = Path(bootstrap_dir)
    
    results = {}
    n_iterations = len(list(bootstrap_dir.glob("iteration_*")))
    
    for i in range(n_iterations):
        iter_dir = bootstrap_dir / f"iteration_{i:02d}"
        
        if not iter_dir.exists():
            continue
        
        results[i] = {
            'factors': pd.read_csv(iter_dir / "factors.csv", index_col=0),
            'variance': pd.read_csv(iter_dir / "variance_explained.csv"),
            'weights': pd.read_csv(iter_dir / "weights.csv", index_col=0),
            'cell_ids': open(iter_dir / "sampled_cell_ids.txt").read().splitlines()
        }
    
    return results

def compute_factor_stability(results):
    """
    Compute factor correlation across bootstrap iterations
    """
    print("\n" + "="*80)
    print("FACTOR STABILITY ANALYSIS")
    print("="*80)
    
    n_iterations = len(results)
    n_factors = results[0]['factors'].shape[1]
    
    # Pairwise factor correlations
    correlation_matrix = np.zeros((n_iterations, n_iterations))
    
    for i in range(n_iterations):
        for j in range(i, n_iterations):
            # Find overlapping cells
            cells_i = set(results[i]['cell_ids'])
            cells_j = set(results[j]['cell_ids'])
            overlap = list(cells_i & cells_j)
            
            if len(overlap) > 10:
                # Subset factors to overlapping cells
                factors_i = results[i]['factors'].loc[overlap]
                factors_j = results[j]['factors'].loc[overlap]
                
                # Average correlation across all factors
                corrs = []
                for f in range(n_factors):
                    corr, _ = spearmanr(factors_i.iloc[:, f], factors_j.iloc[:, f])
                    corrs.append(abs(corr))  # Use absolute correlation
                
                avg_corr = np.mean(corrs)
                correlation_matrix[i, j] = avg_corr
                correlation_matrix[j, i] = avg_corr
    
    # Diagonal is 1
    np.fill_diagonal(correlation_matrix, 1.0)
    
    mean_corr = correlation_matrix[np.triu_indices(n_iterations, k=1)].mean()
    std_corr = correlation_matrix[np.triu_indices(n_iterations, k=1)].std()
    
    print(f"\nFactor stability (Spearman correlation):")
    print(f"  Mean: {mean_corr:.3f} ± {std_corr:.3f}")
    print(f"  Range: [{correlation_matrix[np.triu_indices(n_iterations, k=1)].min():.3f}, "
          f"{correlation_matrix[np.triu_indices(n_iterations, k=1)].max():.3f}]")
    
    return correlation_matrix

def compute_variance_stability(results):
    """
    Compute variance explained stability
    """
    print("\n" + "="*80)
    print("VARIANCE EXPLAINED STABILITY")
    print("="*80)
    
    var_dfs = [r['variance'] for r in results.values()]
    
    # Aggregate by view
    rna_vars = []
    atac_vars = []
    
    for var_df in var_dfs:
        rna_var = var_df[var_df['View'] == 'view0']['R2'].sum()
        atac_var = var_df[var_df['View'] == 'view1']['R2'].sum()
        rna_vars.append(rna_var)
        atac_vars.append(atac_var)
    
    print(f"\nRNA variance explained:")
    print(f"  Mean: {np.mean(rna_vars):.2f}% ± {np.std(rna_vars):.2f}%")
    print(f"  Range: [{np.min(rna_vars):.2f}%, {np.max(rna_vars):.2f}%]")
    
    print(f"\nATAC variance explained:")
    print(f"  Mean: {np.mean(atac_vars):.2f}% ± {np.std(atac_vars):.2f}%")
    print(f"  Range: [{np.min(atac_vars):.2f}%, {np.max(atac_vars):.2f}%]")
    
    return {
        'rna': {'mean': np.mean(rna_vars), 'std': np.std(rna_vars)},
        'atac': {'mean': np.mean(atac_vars), 'std': np.std(atac_vars)}
    }

def plot_stability_results(correlation_matrix, variance_stats, output_dir):
    """
    Generate stability plots
    """
    output_dir = Path(output_dir)
    
    # Factor correlation heatmap
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Correlation matrix
    sns.heatmap(correlation_matrix, annot=True, fmt='.2f', 
                cmap='RdYlGn', vmin=0, vmax=1,
                ax=axes[0], cbar_kws={'label': 'Spearman |r|'})
    axes[0].set_title('Factor Stability Across Iterations')
    axes[0].set_xlabel('Iteration')
    axes[0].set_ylabel('Iteration')
    
    # Variance boxplot
    var_data = pd.DataFrame({
        'RNA': [variance_stats['rna']['mean']] * 2,
        'ATAC': [variance_stats['atac']['mean']] * 2
    })
    
    var_df = pd.melt(var_data, var_name='Modality', value_name='Variance (%)')
    axes[1].bar(['RNA', 'ATAC'], 
                [variance_stats['rna']['mean'], variance_stats['atac']['mean']],
                yerr=[variance_stats['rna']['std'], variance_stats['atac']['std']],
                capsize=5, color=['#E64B35', '#4DBBD5'])
    axes[1].set_ylabel('Variance Explained (%)')
    axes[1].set_title('Variance Stability Across Iterations')
    axes[1].set_ylim([0, 100])
    
    plt.tight_layout()
    plt.savefig(output_dir / 'stability_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\n✓ Saved: {output_dir / 'stability_analysis.png'}")

def parse_args():
    parser = argparse.ArgumentParser(description='Consensus Analysis')
    
    parser.add_argument('--bootstrap_dir', type=str, required=True,
                       help='Bootstrap results directory')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("BOOTSTRAP CONSENSUS ANALYSIS")
    print("="*80)
    
    # Load results
    print(f"\nLoading bootstrap results from: {args.bootstrap_dir}")
    results = load_bootstrap_results(args.bootstrap_dir)
    print(f"✓ Loaded {len(results)} iterations")
    
    # Compute stability metrics
    correlation_matrix = compute_factor_stability(results)
    variance_stats = compute_variance_stability(results)
    
    # Plot results
    plot_stability_results(correlation_matrix, variance_stats, OUTPUT_DIR)
    
    # Save summary
    summary = {
        'n_iterations': len(results),
        'factor_correlation_mean': float(correlation_matrix[np.triu_indices(len(results), k=1)].mean()),
        'factor_correlation_std': float(correlation_matrix[np.triu_indices(len(results), k=1)].std()),
        'variance_stats': variance_stats
    }
    
    with open(OUTPUT_DIR / 'consensus_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✅ Consensus analysis complete!")
    print(f"Results saved to: {OUTPUT_DIR}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
