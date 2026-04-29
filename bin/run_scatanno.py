#!/usr/bin/env python3
"""
run_scatanno.py — scATAnno ATAC cell type annotation for the FORGE pipeline.

Annotates ATAC query cells by integrating with a reference peak atlas via
spectral embedding + Harmony + KNN assignment. Requires a SnapATAC2 h5ad
with fragment_paired in obsm for peak re-counting against the reference.

Monkey-patches are applied for version compatibility within the container
(anndata >=0.10.9, harmonypy >=0.2.0, newer scipy/pandas).
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import snapatac2 as snap

# ======================================================================
# Monkey-patches for container compatibility (snapatac_extended.sif)
# ======================================================================

# 1. AxisArrays: anndata >=0.10.9 requires 'store' kwarg
from anndata._core.aligned_mapping import AxisArrays as _OrigAxisArrays
_orig_init = _OrigAxisArrays.__init__

def _patched_init(self, parent, *, axis, store=None):
    if store is None:
        store = {}
    _orig_init(self, parent, axis=axis, store=store)

_OrigAxisArrays.__init__ = _patched_init

# 2. spectral(): pandas Series can't index sparse matrices in newer scipy
import scATAnno.SnapATAC2_spectral as _snap_spectral
_orig_spectral = _snap_spectral.spectral

def _patched_spectral(data, n_comps=None, features="selected", **kwargs):
    if isinstance(features, str):
        if features in data.var:
            features = np.asarray(data.var[features]).astype(bool)
        else:
            raise NameError("Please call `select_features` first")
    return _orig_spectral(data, n_comps=n_comps, features=features, **kwargs)

_snap_spectral.spectral = _patched_spectral
import scATAnno.scATAnno_integration as _int_mod
_int_mod.spectral = _patched_spectral

# 3. harmony(): harmonypy >=0.2.0 returns (cells, dims) not (dims, cells)
import scATAnno.SnapATAC2_tools as _snap_tools

def _patched_harmony(adata, batch, use_dims=None, use_rep=None, inplace=True, **kwargs):
    import harmonypy
    if use_rep is None:
        use_rep = "X_spectral"
    mat = adata.obsm[use_rep]
    if isinstance(use_dims, int):
        use_dims = range(use_dims)
    mat = mat if use_dims is None else mat[:, use_dims]

    harmony_out = harmonypy.run_harmony(mat, adata.obs, batch, **kwargs)
    corrected = harmony_out.Z_corr
    if corrected.shape[0] == adata.n_obs:
        result = corrected
    else:
        result = corrected.T
    if inplace:
        adata.obsm[use_rep + "_harmony"] = np.asarray(result)
    else:
        return np.asarray(result)

_snap_tools.harmony = _patched_harmony
_int_mod.harmony = _patched_harmony

# 4. raw_assignment(): numpy/pandas scalar vs Series incompatibility
import scATAnno.scATAnno_assignment as _assign_mod

def _patched_raw_assignment(K, neighbors_labels):
    prediction = []
    uncertainty = []
    for i in range(K.shape[0]):
        c_label = neighbors_labels[i, :]
        c_D = K[i, :]
        c_df = pd.DataFrame({'label': c_label, 'dist': c_D})
        p = c_df.groupby('label').sum('dist') / np.sum(c_df['dist'])
        u = 1 - p
        pred_y = u.index[np.argmin(u.values)]
        uncertainty.append(float(u.min().min()))
        prediction.append(pred_y)
    return prediction, uncertainty

_assign_mod.raw_assignment = _patched_raw_assignment

# scATAnno imports (after patches)
from scATAnno.scATAnno_preprocess import load_reference_data
from scATAnno.scATAnno_integration import scATAnno_integrate, scATAnno_harmony
from scATAnno.scATAnno_assignment import scATAnno_annotate


def parse_args():
    p = argparse.ArgumentParser(description="scATAnno ATAC cell type annotation")
    p.add_argument("--query", required=True,
                   help="SnapATAC2 h5ad with fragment_paired in obsm")
    p.add_argument("--reference", required=True,
                   help="scATAnno reference atlas h5ad")
    p.add_argument("--output", default="scatanno_annotations.h5ad",
                   help="Output annotated h5ad")
    p.add_argument("--csv-output", default="scatanno_annotations.csv",
                   help="Output CSV with annotations")
    p.add_argument("--atlas", default="general",
                   help="Atlas name for scATAnno_annotate (e.g., PBMC, brain)")
    p.add_argument("--distance-threshold", type=int, default=95,
                   help="Distance percentile threshold for unknown assignment")
    p.add_argument("--uncertainty-threshold", type=float, default=0.5,
                   help="Uncertainty threshold for unknown assignment")
    p.add_argument("--n-dims", type=int, default=30,
                   help="Number of spectral dimensions")
    p.add_argument("--knn-neighbors", type=int, default=30,
                   help="Number of KNN neighbors for annotation")
    return p.parse_args()


def main():
    args = parse_args()

    # ================================================================
    # STEP 1: Load reference atlas
    # ================================================================
    print("=" * 60)
    print("STEP 1: Loading reference atlas")
    print("=" * 60)
    reference = load_reference_data(args.reference)
    ref_peaks = list(reference.var_names)
    print(f"  Reference: {reference.shape[0]} cells x {reference.shape[1]} peaks")
    print(f"  Cell types ({reference.obs['celltypes'].nunique()}):")
    print(reference.obs['celltypes'].value_counts().to_string(header=False))

    # ================================================================
    # STEP 2: Re-count query fragments against reference peaks
    # ================================================================
    print("\n" + "=" * 60)
    print("STEP 2: Re-counting query fragments against reference peaks")
    print("=" * 60)

    # Handle both AnnDataSet (.h5ads) and individual h5ad.
    # mode='r' keeps the source .h5ads read-only — avoids corrupting the
    # upstream publishDir copy and permits concurrent tasks to share it.
    if args.query.endswith('.h5ads'):
        query_snap = snap.read_dataset(args.query, mode='r')
    else:
        query_snap = snap.read(args.query, backed='r')

    print(f"  Query: {query_snap.shape[0]} cells x {query_snap.shape[1]} features")
    print(f"  obsm keys: {list(query_snap.obsm.keys())}")

    query_peak = snap.pp.make_peak_matrix(
        query_snap,
        use_rep=ref_peaks,
        counting_strategy='insertion',
        inplace=False
    )
    print(f"  Re-counted query: {query_peak.shape[0]} cells x {query_peak.shape[1]} peaks")

    if not args.query.endswith('.h5ads'):
        query_snap.close()

    assert query_peak.shape[1] == reference.shape[1], \
        f"Peak mismatch: query {query_peak.shape[1]} vs ref {reference.shape[1]}"
    print("  Peak sets aligned: OK")

    # ================================================================
    # STEP 3: Prepare query AnnData for scATAnno
    # ================================================================
    print("\n" + "=" * 60)
    print("STEP 3: Preparing query for scATAnno")
    print("=" * 60)

    if hasattr(query_peak, 'to_memory'):
        query = query_peak.to_memory()
    else:
        query = query_peak.copy()

    # Store original barcodes before scATAnno modifies them
    original_barcodes = query.obs_names.to_numpy().copy()

    query.obs['celltypes'] = 'query'
    query.obs['tissue'] = args.atlas
    query.obs['dataset'] = 'query'
    print(f"  Query prepared: {query.shape}")

    # ================================================================
    # STEP 4: Integrate reference + query (spectral embedding)
    # ================================================================
    print("\n" + "=" * 60)
    print("STEP 4: Integrating reference and query (spectral embedding)")
    print("=" * 60)
    integrated = scATAnno_integrate(
        reference_data=reference,
        query_data=query,
        variable_prefix='query',
        dim=args.n_dims,
        sample_size=10000,
        distance_metric='jaccard'
    )
    print(f"  Integrated shape: {integrated.shape}")
    print(f"  Dataset breakdown: {integrated.obs['dataset'].value_counts().to_dict()}")

    # ================================================================
    # STEP 5: Harmony batch correction
    # ================================================================
    print("\n" + "=" * 60)
    print("STEP 5: Running Harmony batch correction")
    print("=" * 60)
    integrated = scATAnno_harmony(integrated, batch_col='dataset', use_dims=args.n_dims)
    print(f"  Harmony complete. obsm keys: {list(integrated.obsm.keys())}")

    print("  Computing UMAP on integrated embedding...")
    sc.pp.neighbors(integrated, use_rep='X_spectral_harmony')
    sc.tl.umap(integrated)
    print("  UMAP done.")

    # ================================================================
    # STEP 6: Split and annotate
    # ================================================================
    print("\n" + "=" * 60)
    print("STEP 6: Running scATAnno annotation (KNN -> distance -> cluster)")
    print("=" * 60)

    ref_mask = integrated.obs['dataset'] == 'Atlas'
    query_mask = integrated.obs['dataset'] == 'query'
    reference_int = integrated[ref_mask].copy()
    query_int = integrated[query_mask].copy()
    print(f"  Reference: {reference_int.shape[0]} cells")
    print(f"  Query: {query_int.shape[0]} cells")

    query_annotated = scATAnno_annotate(
        reference=reference_int,
        query=query_int,
        reference_label_col='celltypes',
        atlas=args.atlas,
        distance_threshold=args.distance_threshold,
        uncertainty_threshold=args.uncertainty_threshold,
        cluster_col=None,
        UMAP=True,
        leiden_resolution=3,
        use_rep='X_spectral_harmony',
        knn_neighbors=args.knn_neighbors
    )

    # ================================================================
    # STEP 7: Results summary
    # ================================================================
    print("\n" + "=" * 60)
    print("ANNOTATION RESULTS")
    print("=" * 60)

    # Primary prediction column
    pred_col = 'pred_y'
    if pred_col in query_annotated.obs.columns:
        counts = query_annotated.obs[pred_col].value_counts()
        print(f"\n  {pred_col} ({len(counts)} types):")
        print(counts.to_string(header=False))

        n_unknown = (query_annotated.obs[pred_col].str.lower() == 'unknown').sum()
        total = query_annotated.shape[0]
        print(f"\n  Unknown fraction ({pred_col}): {n_unknown}/{total} ({100*n_unknown/total:.1f}%)")

    if 'cluster_annotation' in query_annotated.obs.columns:
        ca_counts = query_annotated.obs['cluster_annotation'].value_counts()
        n_unk_ca = (query_annotated.obs['cluster_annotation'].str.lower() == 'unknown').sum()
        print(f"  Unknown fraction (cluster_annotation): {n_unk_ca}/{total} ({100*n_unk_ca/total:.1f}%)")

    if 'Uncertainty_Score' in query_annotated.obs.columns:
        u = query_annotated.obs['Uncertainty_Score'].astype(float)
        print(f"\n  Uncertainty: mean={u.mean():.3f}, median={u.median():.3f}, "
              f"Q95={u.quantile(0.95):.3f}")

    # ================================================================
    # STEP 8: Restore barcodes and save
    # ================================================================
    print("\n" + "=" * 60)
    print("SAVING OUTPUTS")
    print("=" * 60)

    # scATAnno appends '_query' suffix to barcodes — strip it to restore original
    query_annotated.obs_names = pd.Index([
        bc.removesuffix('_query') for bc in query_annotated.obs_names
    ])
    print(f"  Restored barcodes (first 3): {list(query_annotated.obs_names[:3])}")

    # Add a unified cell_type_prediction column for downstream compatibility
    query_annotated.obs['cell_type_prediction'] = query_annotated.obs[pred_col]

    # Save h5ad
    query_annotated.write_h5ad(args.output)
    print(f"  h5ad: {args.output}")

    # Save CSV with annotation columns
    ann_cols = [c for c in query_annotated.obs.columns
                if any(k in c.lower() for k in [
                    'pred', 'uncertainty', 'leiden', 'annotation',
                    'celltype', 'cell_type', 'dataset'
                ])]
    query_annotated.obs[ann_cols].to_csv(args.csv_output)
    print(f"  CSV:  {args.csv_output}")

    print("\nDone!")


if __name__ == "__main__":
    main()
