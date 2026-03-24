#!/usr/bin/env python3
"""
Assign samples to arbitrary test groups for feasibility testing.
Will be replaced with real metadata when available.
"""

import scanpy as sc
import argparse
import json
from pathlib import Path

def assign_groups(h5ad_path, output_path, group_mapping_json):
    """
    Assign test groups to samples in h5ad based on JSON mapping.
    
    Args:
        h5ad_path: Path to annotated h5ad from scANVI
        output_path: Path to output h5ad with groups
        group_mapping_json: Path to JSON with sample->group mapping
    """
    print(f"Loading AnnData from {h5ad_path}...")
    adata = sc.read_h5ad(h5ad_path)
    
    # Load group mapping
    with open(group_mapping_json, 'r') as f:
        sample_to_group = json.load(f)
    
    print(f"Loaded mapping for {len(sample_to_group)} samples")
    
    # Add group column
    if 'sample' not in adata.obs.columns:
        raise ValueError("No 'sample' column found in adata.obs")
    
    adata.obs['condition_group'] = adata.obs['sample'].map(sample_to_group)
    
    # Check for unmapped samples
    unmapped = adata.obs['condition_group'].isna().sum()
    if unmapped > 0:
        print(f"WARNING: {unmapped} cells from unmapped samples")
        # Assign to 'Unknown' group
        adata.obs['condition_group'].fillna('Unknown', inplace=True)
    
    # Print summary
    print("\nGroup assignment summary:")
    print(adata.obs['condition_group'].value_counts())
    
    # Save
    print(f"\nSaving to {output_path}...")
    adata.write_h5ad(output_path)
    
    # Save summary statistics
    summary = {
        'total_cells': int(adata.n_obs),
        'total_samples': int(adata.obs['sample'].nunique()),
        'groups': adata.obs['condition_group'].value_counts().to_dict()
    }
    
    summary_path = Path(output_path).parent / "group_assignment_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to {summary_path}")

def main():
    parser = argparse.ArgumentParser(description='Assign test groups for DE testing')
    parser.add_argument('--input', required=True, help='Input h5ad file')
    parser.add_argument('--output', required=True, help='Output h5ad file')
    parser.add_argument('--mapping', required=True, help='JSON file with sample->group mapping')
    
    args = parser.parse_args()
    
    assign_groups(args.input, args.output, args.mapping)

if __name__ == '__main__':
    main()
