#!/usr/bin/env python3
"""
MOFA+ High-Memory Integration
Single-shot integration using ALL cells
"""

import numpy as np
import pandas as pd
import mudata as md
from pathlib import Path
import argparse
import sys
import json
import warnings
from h5ad_compat import sanitize_mudata

warnings.filterwarnings('ignore')

def run_mofa_high_memory(mudata_path, output_dir, n_factors=10,
                         convergence_mode='fast', seed=42):
    """
    Run MOFA+ integration in HIGH-MEMORY mode (no bootstrapping)
    Uses ALL cells in a single model training
    """
    print(f"\n{'='*70}")
    print("MOFA+ HIGH-MEMORY MODE (Single-Shot Integration)")
    print(f"{'='*70}")

    # -------------------------------------------------------------------------
    # STEP 1: LOAD MUDATA
    # -------------------------------------------------------------------------
    print("\nLoading MuData...")
    mdata = md.read_h5mu(mudata_path)

    print(f"Loaded: {mdata.n_obs:,} cells")
    print(f"  RNA genes: {mdata.mod['rna'].n_vars:,}")
    print(f"  ATAC peaks: {mdata.mod['atac'].n_vars:,}")

    # -------------------------------------------------------------------------
    # STEP 2: FEATURE SELECTION
    # -------------------------------------------------------------------------
    print("\nSelecting features...")
    import scanpy as sc

    # RNA: top 5000 HVGs
    rna = mdata.mod['rna']
    # Ensure RNA is log1p-normalized (seurat HVG flavor requires it; raw counts overflow expm1)
    import numpy as np
    from scipy import sparse
    x_max = rna.X.data.max() if sparse.issparse(rna.X) else rna.X.max()
    if x_max > 50:  # raw counts, not yet normalized
        print("  Normalizing RNA (raw counts detected)...")
        if 'counts' not in rna.layers:
            rna.layers['counts'] = rna.X.copy()
        sc.pp.normalize_total(rna, target_sum=1e4)
        sc.pp.log1p(rna)
    if 'highly_variable' not in rna.var.columns:
        sc.pp.highly_variable_genes(rna, n_top_genes=5000)
    rna_features = rna.var_names[rna.var['highly_variable']]

    # ATAC: top 10000 variable peaks (reuse HVG machinery on the peak matrix)
    atac = mdata.mod['atac']
    if 'highly_variable' not in atac.var.columns:
        sc.pp.highly_variable_genes(atac, n_top_genes=10000)
    atac_features = atac.var_names[atac.var['highly_variable']]

    print(f"  RNA:  {len(rna_features):,} HVGs")
    print(f"  ATAC: {len(atac_features):,} peaks")

    # Subset MuData to selected features
    mdata_subset = md.MuData({
        'rna': rna[:, rna_features].copy(),
        'atac': atac[:, atac_features].copy()
    })

    # Copy cell-level metadata from original MuData (cell type annotations, etc.)
    for col in mdata.obs.columns:
        if col not in mdata_subset.obs.columns:
            mdata_subset.obs[col] = mdata.obs[col].values
    # Also copy modality-level obs (e.g., rna:celltypist_prediction)
    for mod_key in ['rna', 'atac']:
        if mod_key in mdata.mod and mod_key in mdata_subset.mod:
            for col in mdata.mod[mod_key].obs.columns:
                if col not in mdata_subset.mod[mod_key].obs.columns:
                    mdata_subset.mod[mod_key].obs[col] = mdata.mod[mod_key].obs[col].values
    print(f"  Copied metadata to subset: {list(mdata_subset.obs.columns)}")

    # -------------------------------------------------------------------------
    # STEP 3: PREPARE DATA FOR MOFAPY2
    # -------------------------------------------------------------------------
    print("\nPreparing data matrices for MOFA+ (dense) ...")

    # Convert RNA to dense
    rna_sub = mdata_subset.mod['rna']
    if hasattr(rna_sub.X, "toarray"):
        rna_dense = rna_sub.X.toarray()
    else:
        rna_dense = rna_sub.X

    # Convert ATAC to dense
    atac_sub = mdata_subset.mod['atac']
    if hasattr(atac_sub.X, "toarray"):
        atac_dense = atac_sub.X.toarray()
    else:
        atac_dense = atac_sub.X

    print(f"  RNA dense shape:  {rna_dense.shape}")
    print(f"  ATAC dense shape: {atac_dense.shape}")

    # MOFA+ expects a list of views, each a list of groups
    # Here we treat the dataset as a single group with two views: RNA, ATAC
    data_mat = [[rna_dense], [atac_dense]]
    likelihoods = ["gaussian", "gaussian"]

    # -------------------------------------------------------------------------
    # STEP 4: TRAIN MOFA+ MODEL WITH MOFAPY2
    # -------------------------------------------------------------------------
    print(f"\nTraining MOFA+ model with {n_factors} factors...")
    print("This may take from minutes to hours depending on dataset size...")

    from mofapy2.run.entry_point import entry_point

    ent = entry_point()

    # Data options
    # We assume inputs are already reasonably scaled/logged from upstream
    ent.set_data_options(scale_views=False)

    # Load data
    ent.set_data_matrix(
        data_mat,
        likelihoods=likelihoods
    )

    # Model options
    ent.set_model_options(
        factors=n_factors,
        spikeslab_weights=True,
        ard_weights=True
    )

    # Training options
    ent.set_train_options(
        convergence_mode=convergence_mode,
        dropR2=0.001,
        gpu_mode=False,
        seed=seed,
        verbose=True
    )

    # Build and run
    model_file = Path(output_dir) / "model.hdf5"
    ent.build()
    ent.run()
    ent.save(outfile=str(model_file))

    print(f"\n MOFA+ training complete!")
    print(f"Model saved: {model_file}")

    # -------------------------------------------------------------------------
    # STEP 5: EXTRACT FACTORS WITH MOFAX AND ATTACH TO MUDATA
    # -------------------------------------------------------------------------
    print("\nExtracting factors and attaching to MuData...")

    import mofax as mfx

    model = mfx.mofa_model(str(model_file))
    factors_array = model.get_factors(df=False)   # shape: (n_cells, n_factors)
    model.close()

    if factors_array.shape[0] != mdata_subset.n_obs:
        raise RuntimeError(
            f"Number of cells in MOFA factors ({factors_array.shape[0]}) "
            f"does not match MuData ({mdata_subset.n_obs})"
        )

    mdata_subset.obsm["X_mofa"] = factors_array

    # -------------------------------------------------------------------------
    # STEP 6: SAVE INTEGRATED MUDATA AND FACTOR MATRIX
    # -------------------------------------------------------------------------
    output_path = Path(output_dir) / "mofa_integrated.h5mu"
    sanitize_mudata(mdata_subset, output_path)
    mdata_subset.write_h5mu(output_path)
    print(f"Saved integrated MuData: {output_path}")
    n_retained = factors_array.shape[1]

    factors_df = pd.DataFrame(
        factors_array,
        index=mdata_subset.obs_names,
        columns=[f"Factor{i+1}" for i in range(n_retained)]
    )
    factors_csv = Path(output_dir) / "mofa_factors.csv"
    factors_df.to_csv(factors_csv)
    print(f"Saved factor matrix: {factors_csv}")
    # -------------------------------------------------------------------------
    # STEP 7: SAVE BASIC STATS
    # -------------------------------------------------------------------------
    stats = {
        "mode": "high_memory",
        "n_cells": int(mdata_subset.n_obs),
        "n_factors": int(n_retained),
        "n_rna_features": int(len(rna_features)),
        "n_atac_features": int(len(atac_features)),
        "convergence_mode": convergence_mode,
        "seed": int(seed)
    }

    stats_json = Path(output_dir) / "mofa_stats.json"
    with open(stats_json, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"Saved statistics: {stats_json}")

    return output_path

def main():
    parser = argparse.ArgumentParser(description="MOFA+ High-Memory Integration")
    parser.add_argument("--mudata", required=True, help="Path to input MuData")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--n_factors", type=int, default=10, help="Number of factors")
    parser.add_argument(
        "--convergence_mode",
        default="fast",
        choices=["fast", "medium", "slow"],
    )
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print("Running MOFA+ in HIGH-MEMORY mode")
    result = run_mofa_high_memory(
        mudata_path=args.mudata,
        output_dir=args.output_dir,
        n_factors=args.n_factors,
        convergence_mode=args.convergence_mode,
        seed=args.seed,
    )
    print(f" Integration complete: {result}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
