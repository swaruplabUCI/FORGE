#!/usr/bin/env python3
"""
merge_annotations.py 

Merge cell type annotations into the peak matrix while preserving original barcodes.
"""

import argparse
import anndata as ad
import pandas as pd
import numpy as np
import json

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peak-matrix", required=True, help="Input peak matrix h5ad")
    p.add_argument("--annotations", required=True, help="Cell type annotations JSON")
    p.add_argument("--metadata", required=True, help="Sample metadata CSV")
    p.add_argument("--resolution", required=True, help="Leiden resolution to use (e.g., leiden_0_5)")
    p.add_argument("--output", default="peak_matrix_annotated.h5ad", help="Output h5ad")
    return p.parse_args()

def extract_base_sample(sample_str, metadata_samples=None):
    """
    Extract base sample name from compound IDs.
    Examples:
      - "sample1_batch1_batch1" -> "sample1"
      - "L10_Donor_0_nov" -> "L10_Donor_0"
      - "Donor_A_july" -> "Donor_A"
      - "10k_PBMC_10k_batch" -> "10k_PBMC"  (when metadata has "10k_PBMC")
    """
    s = str(sample_str)

    # If we have metadata samples, try progressively shorter prefixes to find a match
    if metadata_samples:
        parts = s.split('_')
        # Try removing trailing parts one at a time (longest match first)
        for i in range(len(parts), 0, -1):
            candidate = '_'.join(parts[:i])
            if candidate in metadata_samples:
                return candidate

    # Fallback: strip known batch keywords
    parts = s.split('_')
    batch_keywords = ['batch', 'june', 'july', 'nov', 'november']
    filtered = [p for p in parts if p.lower() not in batch_keywords]

    if filtered:
        return '_'.join(filtered)

    return parts[0]

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
    
    # Load annotations JSON
    print(f"\nLoading cell type annotations from {args.annotations}")
    with open(args.annotations, 'r') as f:
        annotations_dict = json.load(f)
    
    # Check if the JSON has the resolution as a key or if it's flat (cluster IDs only)
    resolution_key = args.resolution
    
    if resolution_key in annotations_dict:
        # Nested structure: {"leiden_0_5": {"0": "CellType1", ...}}
        cluster_to_celltype = annotations_dict[resolution_key]
        print(f"Using nested structure with resolution key: {resolution_key}")
    else:
        # Flat structure: {"0": "CellType1", "1": "CellType2", ...}
        print(f"Using flat structure (no resolution key found)")
        cluster_to_celltype = annotations_dict
    
    print(f"Cluster to cell type mapping (first 5):")
    for i, (k, v) in enumerate(list(cluster_to_celltype.items())[:5]):
        print(f"  {k} -> {v}")
    
    # Check if the resolution column exists in peak_matrix.obs
    if resolution_key not in peak_matrix.obs.columns:
        raise ValueError(
            f"Resolution column '{resolution_key}' not found in peak_matrix.obs. "
            f"Available columns: {list(peak_matrix.obs.columns)}"
        )
    
    # Map clusters to cell types
    print("\nMapping clusters to cell types...")
    
    # Check if values are dictionaries (nested structure) or simple strings
    sample_value = list(cluster_to_celltype.values())[0]
    if isinstance(sample_value, dict):
        # Extract 'cell_type' field from nested dictionaries
        print("  Detected nested dictionary format with confidence scores")
        cluster_to_celltype_simple = {
            cluster: info['cell_type'] 
            for cluster, info in cluster_to_celltype.items()
        }
    else:
        # Values are already simple strings
        print("  Detected simple string format")
        cluster_to_celltype_simple = cluster_to_celltype
    
    peak_matrix.obs['cell_type'] = (
        peak_matrix.obs[resolution_key]
        .astype(str)
        .map(cluster_to_celltype_simple)
    )
    
    # Check for unmapped clusters
    unmapped = peak_matrix.obs['cell_type'].isna().sum()
    if unmapped > 0:
        print(f"  WARNING: {unmapped} cells could not be mapped to cell types")
        unique_clusters = peak_matrix.obs[resolution_key].unique()
        print(f"  Unique clusters in data: {sorted([str(c) for c in unique_clusters])}")
        print(f"  Clusters in mapping: {sorted(cluster_to_celltype_simple.keys())}")
    
    cell_type_counts = peak_matrix.obs['cell_type'].value_counts()
    print(f"  Cell type distribution:")
    print(cell_type_counts)
    
    # Load and merge metadata
    print(f"\nLoading metadata from {args.metadata}")
    metadata = pd.read_csv(args.metadata)
    print(f"Metadata shape: {metadata.shape}")
    print(f"Metadata columns: {list(metadata.columns)}")
    
    # FIX: Filter metadata to only include individual sample entries (not lane entries)
    if 'sample_type' in metadata.columns:
        print(f"\nFiltering metadata to only 'demux' sample types...")
        metadata = metadata[metadata['sample_type'] == 'demux'].copy()
        print(f"Filtered metadata shape: {metadata.shape}")
    else:
        print("\nWARNING: 'sample_type' column not found in metadata - using all rows")
    
    # Identify the correct sample column in metadata
    sample_col_metadata = None
    for col_name in ['demux_sample', 'sample', 'sample_id', 'Sample', 'sample_name']:
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
    
    # Final check
    print(f"\nFinal peak matrix:")
    print(f"  Shape: {peak_matrix.shape}")
    print(f"  obs columns: {list(peak_matrix.obs.columns)}")
    if 'cell_type' in peak_matrix.obs.columns:
        print(f"  Cell types: {list(peak_matrix.obs['cell_type'].unique())}")
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
