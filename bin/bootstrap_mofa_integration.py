#!/usr/bin/env python3
"""
Bootstrap MOFA+ Integration for Stability Analysis
Based on: https://github.com/bioFAM/MOFA2/blob/master/mofapy2/notebooks/getting_started_python.ipynb
"""

import numpy as np
import pandas as pd
import mudata as md
from pathlib import Path
from mofapy2.run.entry_point import entry_point
import warnings
import argparse
import sys
import json
import psutil
import time

warnings.filterwarnings('ignore')

def log_memory(stage):
    """Log current memory usage"""
    process = psutil.Process()
    mem_info = process.memory_info()
    mem_gb = mem_info.rss / (1024 ** 3)
    
    system_mem = psutil.virtual_memory()
    system_used_gb = system_mem.used / (1024 ** 3)
    system_total_gb = system_mem.total / (1024 ** 3)
    system_pct = system_mem.percent
    
    log_entry = {
        'stage': stage,
        'timestamp': time.time(),
        'process_memory_gb': mem_gb,
        'system_used_gb': system_used_gb,
        'system_total_gb': system_total_gb,
        'system_percent': system_pct
    }
    
    # Append to memory log
    with open('memory_log.jsonl', 'a') as f:
        f.write(json.dumps(log_entry) + '\n')
    
    print(f"[{stage}] Memory: Process={mem_gb:.2f}GB, System={system_used_gb:.2f}/{system_total_gb:.2f}GB ({system_pct:.1f}%)")

def parse_args():
    parser = argparse.ArgumentParser(description='Bootstrap MOFA+ integration')
    parser.add_argument('--mudata_file', type=str, required=True)
    parser.add_argument('--n_iterations', type=int, default=10)
    parser.add_argument('--sample_fraction', type=float, default=0.1)
    parser.add_argument('--n_factors', type=int, default=10)
    parser.add_argument('--output_dir', type=str, default='bootstrap_results')
    parser.add_argument('--seed', type=int, default=42)
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("BOOTSTRAP MOFA INTEGRATION")
    print("="*80)
    print(f"Iterations: {args.n_iterations}")
    print(f"Sample fraction: {args.sample_fraction * 100}%")
    print(f"Random seed: {args.seed}")
    
    # Load full dataset
    log_memory("start_loading_data")
    print(f"\nLoading MuData: {args.mudata_file}")
    mdata = md.read_h5mu(args.mudata_file)
    log_memory("loaded_full_data")
    
    n_cells_total = mdata.n_obs
    n_cells_sample = int(n_cells_total * args.sample_fraction)
    
    print(f"Full dataset: {n_cells_total:,} cells")
    print(f"Will sample: {n_cells_sample:,} cells per iteration")
    
    # Bootstrap iterations
    np.random.seed(args.seed)
    
    for i in range(args.n_iterations):
        print(f"\n{'='*80}")
        print(f"ITERATION {i+1}/{args.n_iterations}")
        print(f"{'='*80}")
        
        log_memory(f"iter_{i}_start")
        
        # Sample cells (stratified by cell type if available)
        print("\nSampling cells...")
        if 'cell_type' in mdata.obs.columns or 'scanvi_prediction' in mdata.obs.columns:
            cell_type_col = 'scanvi_prediction' if 'scanvi_prediction' in mdata.obs.columns else 'cell_type'
            
            # Stratified sampling
            sampled_indices = []
            for cell_type in mdata.obs[cell_type_col].unique():
                ct_indices = np.where(mdata.obs[cell_type_col] == cell_type)[0]
                n_sample_ct = max(1, int(len(ct_indices) * args.sample_fraction))
                if len(ct_indices) > n_sample_ct:
                    sampled = np.random.choice(ct_indices, size=n_sample_ct, replace=False)
                else:
                    sampled = ct_indices
                sampled_indices.extend(sampled)
            
            sampled_indices = np.array(sampled_indices)
            print(f"  Stratified sampling: {len(sampled_indices)} cells across {mdata.obs[cell_type_col].nunique()} cell types")
        else:
            # Random sampling
            sampled_indices = np.random.choice(n_cells_total, size=n_cells_sample, replace=False)
            print(f"  Random sampling: {len(sampled_indices)} cells")
        
        mdata_subset = mdata[sampled_indices, :].copy()
        log_memory(f"iter_{i}_sampled")
        
        # Prepare data for MOFA
        print("\nPreparing MOFA input...")
        rna = mdata_subset.mod['rna']
        atac = mdata_subset.mod['atac']
        
        # Select features (top 2000 genes, top 5000 peaks by expression)
        if hasattr(rna.X, 'toarray'):
            gene_means = np.array(rna.X.mean(axis=0)).flatten()
        else:
            gene_means = rna.X.mean(axis=0)
        
        top_genes_idx = np.argsort(gene_means)[::-1][:2000]
        rna_subset = rna[:, top_genes_idx].copy()
        log_memory(f"iter_{i}_rna_extracted")
        
        if hasattr(atac.X, 'toarray'):
            peak_means = np.array(atac.X.mean(axis=0)).flatten()
        else:
            peak_means = atac.X.mean(axis=0)
        
        top_peaks_idx = np.argsort(peak_means)[::-1][:5000]
        atac_subset = atac[:, top_peaks_idx].copy()
        log_memory(f"iter_{i}_atac_extracted")
        
        # Convert to dense
        rna_dense = rna_subset.X.toarray() if hasattr(rna_subset.X, 'toarray') else rna_subset.X
        atac_dense = atac_subset.X.toarray() if hasattr(atac_subset.X, 'toarray') else atac_subset.X
        
        # CRITICAL: Nested list format for mofapy2
        data_mat = [[rna_dense], [atac_dense]]
        
        # Train MOFA model
        print("\nRunning MOFA+ integration...")
        
        ent = entry_point()
        
        # Set data options BEFORE loading data
        ent.set_data_options(scale_views=False)
        
        # Load data
        ent.set_data_matrix(data_mat, likelihoods=["gaussian", "gaussian"])
        
        # Model options
        ent.set_model_options(
            factors=args.n_factors,
            spikeslab_weights=True,
            ard_weights=True
        )
        
        # Training options
        ent.set_train_options(
            convergence_mode="fast",
            dropR2=0.001,
            gpu_mode=False,
            seed=args.seed + i,  # Different seed per iteration
            verbose=False
        )
        
        # Build and train
        ent.build()
        ent.run()
        
        # Save model
        model_file = output_dir / f"mofa_model_iter_{i}.hdf5"
        ent.save(outfile=str(model_file))
        log_memory(f"iter_{i}_complete")
        
        print(f"✓ Iteration {i+1} complete. Model saved: {model_file.name}")
        
        # Clean up
        del mdata_subset, rna_subset, atac_subset, rna_dense, atac_dense, data_mat, ent
    
    # Save summary
    summary = {
        'n_iterations': args.n_iterations,
        'n_cells_total': n_cells_total,
        'n_cells_per_iteration': n_cells_sample,
        'sample_fraction': args.sample_fraction,
        'n_factors': args.n_factors,
        'seed': args.seed,
        'timestamp': pd.Timestamp.now().isoformat()
    }
    
    summary_file = output_dir / "bootstrap_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n{'='*80}")
    print("✅ BOOTSTRAP COMPLETE!")
    print(f"{'='*80}")
    print(f"Models saved in: {output_dir}")
    print(f"Summary: {summary_file}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
