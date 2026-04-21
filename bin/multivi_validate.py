#!/usr/bin/env python
"""
MultiVI Validation — Biological validation of gap-filled cells.

Validation hierarchy (matching Ashuach et al.):
  1. PRIMARY:    UMAP projection of imputed cells onto existing RNA & ATAC UMAPs
  2. SECONDARY:  50-NN concordance in joint latent space
  3. TERTIARY:   Marker gene dot plots (observed vs imputed)
  4. QUATERNARY: Cell type classification consistency
"""

import argparse
import json
import os
import warnings

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from sklearn.neighbors import NearestNeighbors

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# Module-level single-slot list used by helpers to know the runtime cell_type_key
# without threading it through every function signature. main() sets this once.
_CT_REQUESTED_KEY = []


def parse_args():
    p = argparse.ArgumentParser(description="Validate MultiVI gap-filled cells")
    p.add_argument("--gap_filled_mudata", required=True, help="Gap-filled MuData (.h5mu)")
    p.add_argument("--rna_h5ad", required=True, help="Original RNA h5ad with UMAP (post-integration)")
    p.add_argument("--atac_h5ad", required=True, help="Original ATAC h5ad with UMAP (annotated)")
    p.add_argument("--output_dir", default="validation", help="Output directory")
    p.add_argument("--cell_type_key", default="cell_type", help="Canonical cell type key (RNA side)")
    p.add_argument("--atac_cell_type_key", default="cell_type", help="Cell type key for ATAC")
    p.add_argument("--marker_genes", default=None, help="Comma-separated marker genes (or auto-detect for PBMC)")
    p.add_argument("--k_neighbors", type=int, default=50, help="k for NN concordance")
    return p.parse_args()


# ---------------------------------------------------------------------------
# 1. UMAP Projection
# ---------------------------------------------------------------------------

def project_onto_rna_umap(gap_mdata, rna_ref, cell_type_key, output_dir):
    """
    Project imputed-RNA cells onto the existing RNA UMAP.

    Uses scanpy.tl.ingest to map new cells into the reference embedding.
    """
    print("\n[1a] Projecting imputed-RNA cells onto RNA UMAP...")

    modality = np.array(gap_mdata.obs["modality"])
    atac_only_mask = modality == "accessibility"
    paired_mask = modality == "paired"

    if atac_only_mask.sum() == 0:
        print("  No ATAC-only cells to project. Skipping RNA UMAP projection.")
        return

    # Check if reference has UMAP
    has_umap = "X_umap" in rna_ref.obsm
    if not has_umap:
        print("  Reference RNA h5ad has no X_umap. Computing UMAP...")
        sc.pp.neighbors(rna_ref, n_neighbors=15)
        sc.tl.umap(rna_ref)

    # Build query AnnData from imputed RNA of ATAC-only cells
    imputed_rna_X = gap_mdata.mod["rna"].layers["multivi_imputed"][atac_only_mask]
    if issparse(imputed_rna_X):
        imputed_rna_X = imputed_rna_X.toarray()

    # Align features
    ref_genes = rna_ref.var_names
    query_genes = gap_mdata.mod["rna"].var_names
    common_genes = ref_genes.intersection(query_genes)

    if len(common_genes) < 100:
        print(f"  WARNING: Only {len(common_genes)} common genes. Skipping projection.")
        return

    print(f"  {len(common_genes)} common genes between reference and imputed data")

    # Subset reference and query to common genes
    ref_sub = rna_ref[:, common_genes].copy()
    query_idx = [list(query_genes).index(g) for g in common_genes]
    query_X = imputed_rna_X[:, query_idx]

    query_ad = ad.AnnData(
        X=query_X,
        obs=gap_mdata.obs[atac_only_mask].copy(),
    )
    query_ad.var_names = common_genes

    # Normalize query like reference
    sc.pp.normalize_total(query_ad, target_sum=1e4)
    sc.pp.log1p(query_ad)

    # Ensure reference is preprocessed
    if "neighbors" not in ref_sub.uns:
        sc.pp.highly_variable_genes(ref_sub, n_top_genes=2000, subset=True) if ref_sub.n_vars > 2000 else None
        sc.pp.neighbors(ref_sub)

    # Use ingest to project
    try:
        sc.tl.ingest(query_ad, ref_sub, obs="leiden" if "leiden" in ref_sub.obs.columns else None)

        # Plot: reference + projected
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Panel 1: reference colored by cell type
        ct_col = cell_type_key if cell_type_key in ref_sub.obs.columns else None
        sc.pl.umap(ref_sub, color=ct_col, ax=axes[0], show=False,
                    title="Reference RNA (observed)", frameon=False,
                    legend_loc="right margin", legend_fontsize=6)

        # Panel 2: projected imputed cells
        sc.pl.umap(query_ad, color="leiden" if "leiden" in query_ad.obs.columns else None,
                    ax=axes[1], show=False,
                    title="Imputed RNA (ATAC-only cells, projected)", frameon=False,
                    legend_loc="right margin", legend_fontsize=6)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "validation_rna_umap_projection.pdf"), dpi=150)
        plt.close()
        print("  Saved: validation_rna_umap_projection.pdf")

    except Exception as e:
        print(f"  WARNING: ingest projection failed: {e}")
        print("  Falling back to combined UMAP...")

        # Fallback: compute a fresh UMAP on combined data
        combined = ad.concat([ref_sub, query_ad], join="inner", label="source",
                             keys=["observed", "imputed"])
        sc.pp.highly_variable_genes(combined, n_top_genes=2000, subset=True) if combined.n_vars > 2000 else None
        sc.pp.neighbors(combined, n_neighbors=15)
        sc.tl.umap(combined)

        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(combined, color="source", ax=ax, show=False,
                    title="RNA UMAP: Observed vs Imputed", frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "validation_rna_umap_projection.pdf"), dpi=150)
        plt.close()


def project_onto_atac_umap(gap_mdata, atac_ref, atac_ct_key, output_dir):
    """
    Project imputed-ATAC cells onto the existing ATAC UMAP.
    """
    print("\n[1b] Projecting imputed-ATAC cells onto ATAC UMAP...")

    modality = np.array(gap_mdata.obs["modality"])
    rna_only_mask = modality == "expression"

    if rna_only_mask.sum() == 0:
        print("  No RNA-only cells to project. Skipping ATAC UMAP projection.")
        return

    has_umap = "X_umap" in atac_ref.obsm
    if not has_umap:
        print("  Reference ATAC h5ad has no X_umap. Computing UMAP...")
        if "X_lsi" in atac_ref.obsm:
            sc.pp.neighbors(atac_ref, use_rep="X_lsi")
        else:
            sc.pp.neighbors(atac_ref)
        sc.tl.umap(atac_ref)

    # Build query from imputed ATAC
    imputed_atac_X = gap_mdata.mod["atac"].layers["multivi_imputed"][rna_only_mask]
    if issparse(imputed_atac_X):
        imputed_atac_X = imputed_atac_X.toarray()

    ref_peaks = atac_ref.var_names
    query_peaks = gap_mdata.mod["atac"].var_names
    common_peaks = ref_peaks.intersection(query_peaks)

    if len(common_peaks) < 100:
        print(f"  WARNING: Only {len(common_peaks)} common peaks. Skipping projection.")
        return

    print(f"  {len(common_peaks)} common peaks")

    # Combined UMAP approach (more robust than ingest for ATAC)
    ref_sub = atac_ref[:, common_peaks].copy()
    query_idx = [list(query_peaks).index(p) for p in common_peaks]
    query_X = imputed_atac_X[:, query_idx]

    query_ad = ad.AnnData(X=query_X, obs=gap_mdata.obs[rna_only_mask].copy())
    query_ad.var_names = common_peaks

    combined = ad.concat([ref_sub, query_ad], join="inner", label="source",
                         keys=["observed", "imputed"])

    # TF-IDF + SVD (LSI) for ATAC
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.decomposition import TruncatedSVD

    X = combined.X
    if issparse(X):
        X = X.toarray()
    tfidf = TfidfTransformer().fit_transform(X)
    svd = TruncatedSVD(n_components=min(50, X.shape[1] - 1))
    X_lsi = svd.fit_transform(tfidf)
    combined.obsm["X_lsi"] = X_lsi[:, 1:]  # drop first component

    sc.pp.neighbors(combined, use_rep="X_lsi", n_neighbors=15)
    sc.tl.umap(combined)

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(combined, color="source", ax=ax, show=False,
                title="ATAC UMAP: Observed vs Imputed", frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "validation_atac_umap_projection.pdf"), dpi=150)
    plt.close()
    print("  Saved: validation_atac_umap_projection.pdf")


def plot_joint_umap(gap_mdata, output_dir):
    """Joint MultiVI UMAP colored by imputed/observed status."""
    print("\n[1c] Joint MultiVI UMAP...")

    if "X_MultiVI" not in gap_mdata.obsm:
        print("  No X_MultiVI found, skipping joint UMAP.")
        return

    sc.pp.neighbors(gap_mdata, use_rep="X_MultiVI", n_neighbors=30)
    sc.tl.umap(gap_mdata, min_dist=0.3)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    sc.pl.umap(gap_mdata, color="modality", ax=axes[0], show=False,
                title="Joint UMAP: Modality Source", frameon=False,
                palette={"paired": "#888888", "expression": "steelblue", "accessibility": "coral"})

    # Color by cell type if available
    ct_col = None
    for candidate in ["celltypist_prediction", "rna:celltypist_prediction", "cell_type"]:
        if candidate in gap_mdata.obs.columns:
            ct_col = candidate
            break
    if ct_col:
        sc.pl.umap(gap_mdata, color=ct_col, ax=axes[1], show=False,
                    title="Joint UMAP: Cell Type", frameon=False,
                    legend_loc="right margin", legend_fontsize=6)
    else:
        axes[1].text(0.5, 0.5, "No cell type annotations", ha="center", va="center",
                     transform=axes[1].transAxes)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "validation_joint_umap_imputed.pdf"), dpi=150)
    plt.close()
    print("  Saved: validation_joint_umap_imputed.pdf")


# ---------------------------------------------------------------------------
# 2. 50-NN Concordance
# ---------------------------------------------------------------------------

def compute_nn_concordance(gap_mdata, k=50):
    """
    For each imputed cell, find k nearest neighbors among observed (paired) cells
    in the joint latent space. Report fraction sharing the same cell type.
    """
    print(f"\n[2] Computing {k}-NN concordance in joint latent space...")

    if "X_MultiVI" not in gap_mdata.obsm:
        print("  No X_MultiVI found, skipping NN concordance.")
        return None

    Z = gap_mdata.obsm["X_MultiVI"]
    modality = np.array(gap_mdata.obs["modality"])

    paired_mask = modality == "paired"
    imputed_mask = ~paired_mask

    if imputed_mask.sum() == 0:
        print("  No imputed cells. Skipping.")
        return None

    # Find cell type column — honor the runtime --cell_type_key (canonical: 'cell_type'),
    # then its rna-prefixed variant, then legacy fallbacks for older cached MuData files.
    requested_key = _CT_REQUESTED_KEY[0] if _CT_REQUESTED_KEY else "cell_type"
    ct_col = None
    for candidate in [requested_key, f"rna:{requested_key}", "cell_type", "rna:cell_type",
                      "celltypist_prediction", "rna:celltypist_prediction",
                      "scanvi_prediction", "rna:scanvi_prediction"]:
        if candidate in gap_mdata.obs.columns:
            ct_col = candidate
            break

    if ct_col is None:
        print("  No cell type column found. Skipping NN concordance.")
        return None

    cell_types = gap_mdata.obs[ct_col].values

    # Build NN index on paired cells only
    Z_paired = Z[paired_mask]
    ct_paired = cell_types[paired_mask]
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(Z_paired)

    # Query with imputed cells
    Z_imputed = Z[imputed_mask]
    ct_imputed = cell_types[imputed_mask]
    mod_imputed = modality[imputed_mask]

    _, indices = nn.kneighbors(Z_imputed)

    records = []
    for i in range(len(Z_imputed)):
        neighbor_cts = ct_paired[indices[i]]
        query_ct = ct_imputed[i]

        if pd.isna(query_ct):
            concordance = np.nan
        else:
            concordance = np.mean(neighbor_cts == query_ct)

        records.append({
            "barcode": gap_mdata.obs_names[imputed_mask][i],
            "modality": mod_imputed[i],
            "cell_type": query_ct,
            "nn_concordance": concordance,
        })

    df = pd.DataFrame(records)

    # Summary
    for mod in ["expression", "accessibility"]:
        mod_df = df[df["modality"] == mod]
        if len(mod_df) > 0:
            print(f"  {mod}: median concordance = {mod_df['nn_concordance'].median():.3f}")

    df.to_csv(os.path.join(gap_mdata.uns.get("output_dir", "validation"), "validation_nn_concordance.csv"), index=False)
    return df


# ---------------------------------------------------------------------------
# 3. Marker Gene Dot Plots
# ---------------------------------------------------------------------------

PBMC_MARKERS = {
    "B cells": ["CD79A", "MS4A1", "CD19"],
    "T cells": ["CD3E", "CD3D", "IL7R"],
    "CD8 T": ["CD8A", "CD8B"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "Monocytes": ["CD14", "LYZ", "S100A9"],
    "DC": ["FCER1A", "CST3"],
    "Platelets": ["PPBP", "PF4"],
}


def _sanitize_groupby_column(adata, col):
    """Prep a dotplot groupby column: strip literal 'nan' string artifact, drop
    NaN cells, and rebuild Categorical so only observed categories are enumerated.

    Returns (new_adata_subset, n_unique_groups). new_adata_subset may be an
    empty view when every cell is unlabeled.
    """
    s = adata.obs[col].astype(object)
    # Treat a historical sanitizer artifact ('nan' as string) as true NaN.
    s = s.where(s != "nan", np.nan)
    keep = s.notna().values
    if not keep.all():
        adata = adata[keep].copy()
        s = s[keep]
    adata.obs[col] = pd.Categorical(s.astype(str))
    return adata, int(adata.obs[col].nunique())


def _nn_transfer_labels(gap_mdata, source_mask, target_mask, source_labels, k=10):
    """Majority-vote NN transfer in X_MultiVI space from source → target cells.

    Returns np.ndarray of transferred labels (str) aligned with target_mask cells;
    NaN where the target cell has no valid source neighbors.
    """
    Z = gap_mdata.obsm.get("X_MultiVI")
    if Z is None or source_mask.sum() == 0 or target_mask.sum() == 0:
        return np.array([np.nan] * int(target_mask.sum()), dtype=object)
    k = min(k, int(source_mask.sum()))
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(Z[source_mask])
    _, ind = nn.kneighbors(Z[target_mask])
    source_labels = np.asarray(source_labels)
    out = []
    for i in range(int(target_mask.sum())):
        neigh = source_labels[ind[i]]
        valid = [x for x in neigh if not (x is None or (isinstance(x, float) and np.isnan(x)))]
        if valid:
            out.append(pd.Series(valid).mode().iloc[0])
        else:
            out.append(np.nan)
    return np.array(out, dtype=object)


def plot_marker_dotplots(gap_mdata, marker_genes, cell_type_key, output_dir):
    """Side-by-side dot plots: observed (paired) vs imputed (ATAC-only) RNA.

    Observed panel groupby: cell_type_key (RNA annotation).
    Imputed panel groupby (hybrid): primary atac:cell_type (scATAnno,
    independent of MultiVI); fallback NN-transfer in X_MultiVI space from
    paired cells' RNA labels. Both label columns are persisted on imputed_ad
    as 'groupby_scatanno' and 'groupby_nn_transfer' for inspection.
    Panels are skipped (not crashed) when <2 groups remain after sanitation.
    """
    print("\n[3] Marker gene dot plots...")

    modality = np.array(gap_mdata.obs["modality"])
    if "multivi_imputed" not in gap_mdata.mod["rna"].layers:
        print("  No multivi_imputed layer found. Skipping.")
        return

    paired_mask = modality == "paired"
    atac_only_mask = modality == "accessibility"
    if atac_only_mask.sum() == 0:
        print("  No imputed RNA cells. Skipping dot plots.")
        return

    rna_X = gap_mdata.mod["rna"].X[paired_mask]
    if issparse(rna_X):
        rna_X = rna_X.toarray()
    observed_ad = ad.AnnData(
        X=rna_X,
        obs=gap_mdata.obs[paired_mask].copy(),
        var=gap_mdata.mod["rna"].var.copy(),
    )
    observed_ad.obs_names = gap_mdata.obs_names[paired_mask]

    imputed_X = gap_mdata.mod["rna"].layers["multivi_imputed"][atac_only_mask]
    if issparse(imputed_X):
        imputed_X = imputed_X.toarray()
    imputed_ad = ad.AnnData(
        X=imputed_X,
        obs=gap_mdata.obs[atac_only_mask].copy(),
        var=gap_mdata.mod["rna"].var.copy(),
    )
    imputed_ad.obs_names = gap_mdata.obs_names[atac_only_mask]

    available_markers = [g for g in marker_genes if g in observed_ad.var_names]
    if len(available_markers) == 0:
        print("  No marker genes found in var_names. Skipping.")
        return
    print(f"  Using {len(available_markers)} marker genes")

    sc.pp.normalize_total(observed_ad, target_sum=1e4); sc.pp.log1p(observed_ad)
    sc.pp.normalize_total(imputed_ad, target_sum=1e4); sc.pp.log1p(imputed_ad)

    # ---- Observed panel groupby: RNA cell_type ----
    observed_col = None
    for candidate in (cell_type_key, f"rna:{cell_type_key}", "celltypist_prediction"):
        if candidate in observed_ad.obs.columns:
            observed_col = candidate
            break
    if observed_col is None:
        print("  No RNA cell type column for observed panel. Skipping dot plots.")
        return
    observed_fixed, n_obs_groups = _sanitize_groupby_column(observed_ad, observed_col)

    # ---- Imputed panel groupby: hybrid (scATAnno primary, NN-transfer fallback) ----
    # Persist both candidate labels on imputed_ad for inspection.
    atac_candidates = ("atac:cell_type", "atac_cell_type", "atac:celltypist_prediction")
    imputed_scatanno_raw = None
    for candidate in atac_candidates:
        if candidate in imputed_ad.obs.columns:
            s = imputed_ad.obs[candidate].astype(object).where(
                imputed_ad.obs[candidate].astype(object) != "nan", np.nan
            )
            if s.notna().sum() > 0:
                imputed_scatanno_raw = (candidate, s.values)
                break
    if imputed_scatanno_raw is not None:
        imputed_ad.obs["groupby_scatanno"] = imputed_scatanno_raw[1]

    # NN-transfer labels (always computed when X_MultiVI is available, for comparison).
    paired_labels_raw = pd.Series(
        gap_mdata.obs[observed_col].astype(object).values
        if observed_col in gap_mdata.obs.columns
        else gap_mdata.mod["rna"].obs[observed_col].astype(object).values,
        index=gap_mdata.obs_names,
    )
    paired_labels_raw = paired_labels_raw.where(paired_labels_raw != "nan", np.nan)
    paired_labels = paired_labels_raw.values[paired_mask]
    nn_labels = _nn_transfer_labels(gap_mdata, paired_mask, atac_only_mask, paired_labels, k=10)
    imputed_ad.obs["groupby_nn_transfer"] = nn_labels

    # Primary = scATAnno; fallback = NN-transfer when scATAnno gives <2 groups.
    imputed_col = None
    imputed_label_source = None
    if "groupby_scatanno" in imputed_ad.obs.columns:
        s = imputed_ad.obs["groupby_scatanno"].dropna()
        if s.nunique() >= 2:
            imputed_col = "groupby_scatanno"
            imputed_label_source = imputed_scatanno_raw[0]
    if imputed_col is None:
        s = pd.Series(nn_labels).dropna()
        if s.nunique() >= 2:
            imputed_col = "groupby_nn_transfer"
            imputed_label_source = f"NN-transfer from {observed_col}"
    imputed_fixed, n_imp_groups = (imputed_ad, 0)
    if imputed_col is not None:
        imputed_fixed, n_imp_groups = _sanitize_groupby_column(imputed_ad, imputed_col)
    print(f"  observed groupby='{observed_col}' n={observed_fixed.n_obs} groups={n_obs_groups}")
    print(f"  imputed  groupby='{imputed_col}' (source='{imputed_label_source}') "
          f"n={imputed_fixed.n_obs} groups={n_imp_groups}")

    plot_observed = observed_fixed.n_obs > 0 and n_obs_groups >= 2
    plot_imputed = imputed_fixed.n_obs > 0 and n_imp_groups >= 2
    if not plot_observed and not plot_imputed:
        print("  Neither panel has >=2 groups after sanitation. Skipping dotplot.")
        return

    n_panels = int(plot_observed) + int(plot_imputed)
    fig, axes = plt.subplots(n_panels, 1, figsize=(max(12, len(available_markers) * 0.8), 5 * n_panels))
    if n_panels == 1:
        axes = [axes]
    i = 0
    if plot_observed:
        sc.pl.dotplot(observed_fixed, var_names=available_markers, groupby=observed_col,
                      ax=axes[i], show=False, title=f"Observed RNA (groupby={observed_col})")
        i += 1
    if plot_imputed:
        sc.pl.dotplot(imputed_fixed, var_names=available_markers, groupby=imputed_col,
                      ax=axes[i], show=False,
                      title=f"Imputed RNA (ATAC-only; groupby={imputed_col}, source={imputed_label_source})")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "validation_marker_dotplots.pdf"), dpi=150)
    plt.close()
    print("  Saved: validation_marker_dotplots.pdf")


# ---------------------------------------------------------------------------
# 4. Cell Type Classification Consistency
# ---------------------------------------------------------------------------

def check_celltype_consistency(gap_mdata, cell_type_key, atac_ct_key, output_dir):
    """
    For imputed cells, compare cell type labels from different sources.

    - ATAC-only cells: compare scATAnno ATAC label vs label inferred from
      imputed RNA (if available from CellTypist or from NN transfer).
    - RNA-only cells: compare RNA cell type vs ATAC NN label.
    """
    print("\n[4] Cell type classification consistency...")

    modality = np.array(gap_mdata.obs["modality"])

    # For ATAC-only cells, they should have ATAC-derived cell types
    atac_only_mask = modality == "accessibility"
    if atac_only_mask.sum() == 0:
        print("  No ATAC-only cells for consistency check.")
        return None

    # Check what annotations are available
    obs = gap_mdata.obs

    # Try to find ATAC cell type for ATAC-only cells
    atac_ct = None
    for candidate in [atac_ct_key, f"atac:{atac_ct_key}", "atac:cell_type", "scatanno_prediction"]:
        if candidate in obs.columns:
            atac_ct = obs[candidate].values[atac_only_mask]
            print(f"  Using ATAC cell type from: {candidate}")
            break

    # Transfer cell type from NN in joint latent space
    if "X_MultiVI" in gap_mdata.obsm:
        Z = gap_mdata.obsm["X_MultiVI"]
        paired_mask = modality == "paired"

        # Find RNA cell type for paired cells
        rna_ct_col = None
        for candidate in [cell_type_key, f"rna:{cell_type_key}", "celltypist_prediction"]:
            if candidate in obs.columns:
                rna_ct_col = candidate
                break

        if rna_ct_col and paired_mask.sum() > 0:
            nn = NearestNeighbors(n_neighbors=10, metric="euclidean").fit(Z[paired_mask])
            _, indices = nn.kneighbors(Z[atac_only_mask])

            paired_cts = obs[rna_ct_col].values[paired_mask]
            # Majority vote from 10-NN
            nn_predicted = []
            for i in range(atac_only_mask.sum()):
                neighbor_cts = paired_cts[indices[i]]
                valid = [ct for ct in neighbor_cts if not pd.isna(ct)]
                if valid:
                    nn_predicted.append(pd.Series(valid).mode().iloc[0])
                else:
                    nn_predicted.append(np.nan)
            nn_predicted = np.array(nn_predicted)

            # Compare ATAC cell type vs NN-predicted RNA cell type
            if atac_ct is not None:
                valid = ~pd.isna(atac_ct) & ~pd.isna(nn_predicted)
                if valid.sum() > 0:
                    agreement = np.mean(atac_ct[valid] == nn_predicted[valid])
                    print(f"  ATAC vs NN-RNA agreement: {agreement:.3f} ({valid.sum()} cells)")

                    result_df = pd.DataFrame({
                        "barcode": gap_mdata.obs_names[atac_only_mask],
                        "atac_cell_type": atac_ct,
                        "nn_rna_cell_type": nn_predicted,
                        "agree": atac_ct == nn_predicted,
                    })
                    result_df.to_csv(
                        os.path.join(output_dir, "validation_celltype_agreement.csv"),
                        index=False,
                    )
                    return result_df

    print("  Insufficient annotations for consistency check.")
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Stash runtime cell_type_key so helpers (e.g., compute_nn_concordance) can
    # look it up without threading it through every signature.
    _CT_REQUESTED_KEY.clear()
    _CT_REQUESTED_KEY.append(args.cell_type_key)

    print("Loading gap-filled MuData...")
    gap_mdata = mu.read_h5mu(args.gap_filled_mudata)
    gap_mdata.uns["output_dir"] = args.output_dir

    print("Loading reference RNA h5ad...")
    rna_ref = sc.read_h5ad(args.rna_h5ad)

    print("Loading reference ATAC h5ad...")
    atac_ref = sc.read_h5ad(args.atac_h5ad)

    modality = np.array(gap_mdata.obs["modality"])
    n_paired = (modality == "paired").sum()
    n_rna_only = (modality == "expression").sum()
    n_atac_only = (modality == "accessibility").sum()
    print(f"  {n_paired} paired, {n_rna_only} RNA-only, {n_atac_only} ATAC-only")

    # 1. UMAP projections
    project_onto_rna_umap(gap_mdata, rna_ref, args.cell_type_key, args.output_dir)
    project_onto_atac_umap(gap_mdata, atac_ref, args.atac_cell_type_key, args.output_dir)
    plot_joint_umap(gap_mdata, args.output_dir)

    # 2. NN concordance
    nn_df = compute_nn_concordance(gap_mdata, k=args.k_neighbors)
    if nn_df is not None:
        nn_df.to_csv(os.path.join(args.output_dir, "validation_nn_concordance.csv"), index=False)

    # 3. Marker dot plots
    if args.marker_genes:
        markers = [g.strip() for g in args.marker_genes.split(",")]
    else:
        # Default PBMC markers
        markers = []
        for genes in PBMC_MARKERS.values():
            markers.extend(genes)
    plot_marker_dotplots(gap_mdata, markers, args.cell_type_key, args.output_dir)

    # 4. Cell type consistency
    ct_df = check_celltype_consistency(
        gap_mdata, args.cell_type_key, args.atac_cell_type_key, args.output_dir
    )

    # Summary
    summary = {
        "n_paired": int(n_paired),
        "n_rna_only": int(n_rna_only),
        "n_atac_only": int(n_atac_only),
        "nn_concordance_median": float(nn_df["nn_concordance"].median()) if nn_df is not None else None,
        "celltype_agreement": float(ct_df["agree"].mean()) if ct_df is not None else None,
    }
    with open(os.path.join(args.output_dir, "validation_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nValidation complete. Results in {args.output_dir}/")


if __name__ == "__main__":
    main()
