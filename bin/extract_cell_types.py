#!/usr/bin/env python3
"""
Extract cell types from SCANVI h5ad with minimum cell threshold
"""
import argparse
import scanpy as sc
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad', required=True, help='Input h5ad file')
    parser.add_argument('--cell_type_key', required=True, help='Column name for cell types')
    parser.add_argument('--min_cells', type=int, default=100, help='Minimum cells per type')
    parser.add_argument('--output', required=True, help='Output text file with cell types')
    parser.add_argument('--counts_output', required=True, help='Output CSV with counts')
    args = parser.parse_args()
    
    # Load h5ad
    print(f"Loading {args.h5ad}...")
    adata = sc.read_h5ad(args.h5ad)
    
    # Count cells per type
    cell_type_counts = adata.obs[args.cell_type_key].value_counts()
    
    # Filter by minimum cells
    valid_types = cell_type_counts[cell_type_counts >= args.min_cells].index.tolist()
    
    print(f"Found {len(valid_types)} cell types with >= {args.min_cells} cells")
    
    # Write valid cell types
    with open(args.output, 'w') as f:
        for ct in valid_types:
            f.write(f"{ct}\n")
    
    # Write counts summary
    counts_df = pd.DataFrame({
        'cell_type': cell_type_counts.index,
        'n_cells': cell_type_counts.values,
        'included': cell_type_counts.index.isin(valid_types)
    })
    counts_df.to_csv(args.counts_output, index=False)
    
    print(f"Cell types written to {args.output}")
    print(f"Counts summary written to {args.counts_output}")

if __name__ == '__main__':
    main()
