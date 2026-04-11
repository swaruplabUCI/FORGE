#!/usr/bin/env python3
"""
Prepare reference and query for SCVI/SCANVI integration  v2.2.0

FIX-33b: Changed HVG batch_key from 'source' to 'sample'.
  When run on the reference alone, 'source' is a constant ('reference'),
  so batch_key had no effect.  Using 'sample' (Donor ID) gives
  batch-aware HVG selection across donors.
"""

import argparse
import scanpy as sc
import os
import numpy as np
from h5ad_compat import sanitize_adata

def prepare_reference_and_query(ref_path, query_path, species, results_dir, n_top_genes=12000):
    """
    Prepare reference and query by:
    1. Loading both datasets
    2. Subsetting to common genes
    3. Normalizing and selecting HVGs on reference
    4. Saving prepared datasets
    """
    print("="*60)
    print("PREPARING REFERENCE AND QUERY FOR SCVI/SCANVI")
    print("="*60)

    # FIX-29b: Handle ref_path being a directory (config uses ref_dir_*)
    #          Auto-detect the correct .h5ad file inside the directory.
    if os.path.isdir(ref_path):
        print(f"ref_path is a directory: {ref_path}")
        h5ad_files = sorted([f for f in os.listdir(ref_path) if f.endswith('.h5ad')])
        print(f"  Found {len(h5ad_files)} h5ad files: {h5ad_files}")

        if not h5ad_files:
            raise FileNotFoundError(f"No .h5ad files found in reference directory: {ref_path}")

        # Priority order for file selection:
        #   1. SEAAD *final-nuclei* (preferred: QC-filtered, latest)
        #   2. SEAAD *all-nuclei*
        #   3. Reference *final-nuclei*
        #   4. Reference *all-nuclei*
        #   5. Any .h5ad file (fallback)
        priority_patterns = [
            'SEAAD_MTG_RNAseq_final-nuclei',
            'SEAAD_MTG_RNAseq_all-nuclei',
            'Reference_MTG_RNAseq_final-nuclei',
            'Reference_MTG_RNAseq_all-nuclei',
        ]

        selected = None
        for pattern in priority_patterns:
            matches = [f for f in h5ad_files if pattern in f]
            if matches:
                # If multiple matches for same pattern, pick newest (last alphabetically by date)
                selected = sorted(matches)[-1]
                break

        if selected is None:
            selected = h5ad_files[0]

        ref_path = os.path.join(ref_path, selected)
        print(f"  Selected reference file: {selected}")

    # Load reference
    print(f"Loading reference from {ref_path}")
    ref = sc.read_h5ad(ref_path)

    # Load query
    print(f"Loading query from {query_path}")
    adata_query = sc.read_h5ad(query_path)

    # ============================================
    # CRITICAL: SET UP BATCH KEYS FOR SCVI
    # ============================================

    # Label source/batch column for tracking
    adata_query.obs['source'] = 'query'
    ref.obs['source'] = 'reference'

    # CRITICAL: Preserve donor-level batch information for scVI
    # SEAAD uses 'Donor ID' (with space) as the donor identifier
    if 'Donor ID' in ref.obs.columns:
        ref.obs['sample'] = ref.obs['Donor ID'].astype(str)
        n_donors = ref.obs['sample'].nunique()
        print(f"✓ Using SEAAD 'Donor ID' as batch key: {n_donors} unique donors")
    elif 'donor_id' in ref.obs.columns:
        ref.obs['sample'] = ref.obs['donor_id'].astype(str)
        n_donors = ref.obs['sample'].nunique()
        print(f"✓ Using 'donor_id' as batch key: {n_donors} unique donors")
    else:
        print("⚠ WARNING: No donor column found in reference, using 'reference' as single batch")
        ref.obs['sample'] = 'reference'

    # For query, preserve original sample information
    if 'Sample' in adata_query.obs.columns:
        adata_query.obs['sample'] = adata_query.obs['Sample'].astype(str)
        print(f"✓ Query samples: {adata_query.obs['sample'].nunique()} unique samples")
    elif 'sample' not in adata_query.obs.columns:
        print("⚠ WARNING: No sample column in query, using 'query' as single batch")
        adata_query.obs['sample'] = 'query'

    print(f"\n{'='*60}")
    print(f"Batch correction summary:")
    print(f"  Reference batches: {ref.obs['sample'].nunique()}")
    print(f"  Query batches: {adata_query.obs['sample'].nunique()}")
    print(f"  Total batches for scVI: {ref.obs['sample'].nunique() + adata_query.obs['sample'].nunique()}")
    print(f"{'='*60}\n")

    # Species-specific filtering and label key
    if species == 'human':
        # Drop rows with NaN in SEAAD annotation columns
        print("Filtering SEAAD reference...")
        ref = ref[
            ~ref.obs['Subclass'].isna() &
            ~ref.obs['Class'].isna() &
            ~ref.obs['Supertype'].isna() &
            ~ref.obs['Supertype (non-expanded)'].isna()
        ]
        SCANVI_LABELS_KEY = 'Subclass'
    else:  # mouse
        print("Filtering Allen reference...")
        # Support both column naming conventions: 'subclass_label'/'class_label' (ABC Atlas)
        # and 'subclass'/'class' (older Allen VISp reference)
        sub_col = 'subclass_label' if 'subclass_label' in ref.obs.columns else 'subclass'
        cls_col = 'class_label' if 'class_label' in ref.obs.columns else 'class'
        ref = ref[
            ~ref.obs[sub_col].isna() &
            ~ref.obs[cls_col].isna()
        ]
        SCANVI_LABELS_KEY = sub_col

    print(f"Reference after filtering: {ref.n_obs} cells, {ref.n_vars} genes")

    # Subsample large references to keep training tractable
    max_ref_cells = 50000
    if ref.n_obs > max_ref_cells:
        print(f"\nSubsampling reference from {ref.n_obs} to {max_ref_cells} cells (stratified by {SCANVI_LABELS_KEY})...")
        sc.pp.subsample(ref, n_obs=max_ref_cells)
        print(f"Reference after subsampling: {ref.n_obs} cells")

    # ============================================
    # CRITICAL: SUBSET TO COMMON GENES
    # ============================================
    print("\nSubsetting to common genes...")

    # Detect if reference has raw counts or is already normalized/scaled.
    # Raw counts: non-negative integers. Pre-processed: floats with negatives or small decimals.
    import scipy.sparse as sp
    ref_x_sample = ref.X[:min(100, ref.n_obs), :min(100, ref.n_vars)]
    if sp.issparse(ref_x_sample):
        ref_x_sample = ref_x_sample.toarray()
    ref_is_raw = (ref_x_sample >= 0).all() and (ref_x_sample == ref_x_sample.astype(int)).all()
    print(f"Reference X appears {'raw counts' if ref_is_raw else 'pre-processed (normalized/scaled)'}")

    # Save raw counts BEFORE subsetting
    if ref_is_raw:
        ref.layers["counts"] = ref.X.copy()
    if 'counts' not in adata_query.layers:
        adata_query.layers["counts"] = adata_query.X.copy()

    # Find common genes
    common_var_names = ref.var_names.intersection(adata_query.var_names)
    print(f"Common genes: {len(common_var_names)}")

    if len(common_var_names) == 0:
        raise ValueError("No common genes found between reference and query!")

    # Subset both datasets
    ref = ref[:, common_var_names].copy()
    adata_query = adata_query[:, common_var_names].copy()

    # ============================================
    # NORMALIZE AND SELECT HVGs ON REFERENCE
    # ============================================
    if ref_is_raw:
        print("\nNormalizing reference and selecting HVGs...")
        sc.pp.normalize_total(ref, target_sum=1e4)
        sc.pp.log1p(ref)
    else:
        print("\nReference already normalized — skipping normalize_total + log1p")

    # FIX-33b: Select HVGs with donor-level batch awareness.
    if ref_is_raw:
        sc.pp.highly_variable_genes(
            ref,
            n_top_genes=n_top_genes,
            batch_key="sample",
            subset=False
        )
    else:
        # Pre-processed reference: already HVG-selected, mark all common genes as HVG
        ref.var['highly_variable'] = True
        print(f"Pre-processed reference: marking all {ref.n_vars} common genes as HVG")

    n_hvgs = ref.var['highly_variable'].sum()
    print(f"Selected {n_hvgs} highly variable genes")

    # Save prepared datasets
    ref_out = os.path.join(results_dir, f'Refmapping_{n_top_genes}HVG_{species}.h5ad')
    query_out = os.path.join(results_dir, 'Refmapping_query.h5ad')
    label_out = os.path.join(results_dir, 'label_key.txt')

    sanitize_adata(ref, ref_out)
    ref.write(ref_out)
    sanitize_adata(adata_query, query_out)
    adata_query.write(query_out)

    # Save label key for downstream use
    with open(label_out, 'w') as f:
        f.write(SCANVI_LABELS_KEY)

    print(f"\nSaved prepared reference: {ref_out}")
    print(f"Saved prepared query: {query_out}")
    print(f"Saved label key: {label_out}")

    return ref, adata_query, SCANVI_LABELS_KEY

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', required=True, help='Query h5ad file (concatenated)')
    parser.add_argument('--species', required=True, choices=['human', 'mouse'])
    parser.add_argument('--ref_path', required=True, help='Path to reference h5ad file or directory')
    parser.add_argument('--results_dir', required=True)
    parser.add_argument('--n_top_genes', type=int, default=12000)

    args = parser.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)

    # Prepare reference and query
    prepare_reference_and_query(
        args.ref_path,
        args.query,
        args.species,
        args.results_dir,
        args.n_top_genes
    )

    print("\n" + "="*60)
    print("REFERENCE PREPARATION COMPLETE!")
    print("="*60)

if __name__ == '__main__':
    main()
