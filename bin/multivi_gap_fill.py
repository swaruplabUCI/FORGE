#!/usr/bin/env python
"""
MultiVI Gap-Fill — Recover cells lost to QC by imputing their missing modality.

Procedure:
  1. Load the pre-intersection RNA and ATAC h5ad files to identify gap cells
     (cells present in only one modality after QC).
  2. Build an expanded MuData that includes gap cells with the missing modality
     zeroed out and a modality indicator set.
  3. Train MultiVI on the expanded dataset (gap cells included during training).
  4. Impute the missing modality using get_normalized_expression/accessibility.
  5. Score confidence via posterior sampling (n_samples=25).
  6. Output the gap-filled MuData with imputed layers and flags.

NOTE: This module does NOT re-run any downstream pipeline modules.
The gap-filled MuData is produced for future analysis.
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
import scvi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse, csr_matrix
from scvi.model import MULTIVI
from scvi.train._callbacks import LoudEarlyStopping

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def _sanitize_mudata_obs_for_h5py(mdata):
    """Coerce object-dtype obs columns to concrete types before write_h5mu.

    pd.concat of bool + NaN (from build_expanded_mudata gap-cell rows) yields
    object dtype; h5py then tries vlen-string serialization and fails on
    non-string scalars (e.g. 'rna:predicted_doublets' from rna_qc scrublet).
    Coercion preserves NaN:
      - All-bool non-null → bool with NaN→False.
      - Otherwise → pd.Categorical with NaN as code -1 (h5py serializes
        categorical as int codes + categories array; no vlen-string failure
        and downstream readers see pd.NA / NaN for missing, not the literal
        string 'nan' which would pollute groupby operations.
    """
    def _fix(obs_df):
        for col in list(obs_df.columns):
            s = obs_df[col]
            if s.dtype == object:
                nn = s.dropna()
                if len(nn) > 0 and nn.map(lambda x: isinstance(x, (bool, np.bool_))).all():
                    obs_df[col] = s.fillna(False).astype(bool)
                else:
                    # Stringify non-null values only; NaN stays NaN so that
                    # pd.Categorical assigns code -1 (serialized/read as NaN).
                    str_vals = [None if pd.isna(x) else str(x) for x in s]
                    obs_df[col] = pd.Categorical(str_vals)

    _fix(mdata.obs)
    if hasattr(mdata, "mod"):
        for _m in mdata.mod:
            _fix(mdata.mod[_m].obs)


def parse_args():
    p = argparse.ArgumentParser(description="MultiVI gap-filling for QC-filtered cells")
    p.add_argument("--mudata", required=True, help="Paired-only integrated MuData (.h5mu)")
    p.add_argument("--rna_h5ad", required=True, help="Post-QC RNA h5ad (all cells, pre-intersection)")
    p.add_argument("--atac_h5ad", required=True, help="Post-QC ATAC h5ad (all cells, pre-intersection)")
    p.add_argument("--output_dir", default="gap_fill", help="Output directory")
    p.add_argument("--n_epochs", type=int, default=200, help="Training epochs")
    p.add_argument("--batch_key", default="sample_id", help="Batch key")
    p.add_argument("--n_latent", type=int, default=20, help="Latent dims")
    p.add_argument("--modality_weights", default="equal")
    p.add_argument("--modality_penalty", default="Jeffreys")
    p.add_argument("--n_posterior_samples", type=int, default=25, help="Posterior samples for confidence")
    p.add_argument("--min_confidence", type=float, default=0.3, help="Min confidence threshold for flagging")
    p.add_argument("--cell_type_key", default="celltypist_prediction", help="Cell type key")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Barcode parsing utilities
# ---------------------------------------------------------------------------

import re

_SUFFIX_RE = re.compile(r"-\d+$")


def _strip_gem_well(bc):
    """Strip trailing -N gem well suffix (10x convention) for matching."""
    return _SUFFIX_RE.sub("", bc)


def _parse_colon_barcodes(obs_names):
    """Parse SAMPLE:BARCODE format → dict{(sample, norm_barcode): obs_name}."""
    lookup = {}
    for name in obs_names:
        if ":" in name:
            sample, raw_bc = name.split(":", 1)
            lookup[(sample, _strip_gem_well(raw_bc))] = name
        else:
            lookup[("", _strip_gem_well(name))] = name
    return lookup


def _parse_rna_concat_barcodes(rna_adata):
    """Parse RNA concat barcodes (BARCODE-SAMPLE from index_unique='-').

    Uses the 'sample' or 'sample_id' obs column to determine each cell's
    sample, then strips the '-SAMPLE' suffix from the obs_name to recover
    the raw barcode.

    Returns dict{(sample, norm_barcode): obs_name}.
    """
    # Determine sample column
    if "sample_id" in rna_adata.obs.columns:
        sample_col = "sample_id"
    elif "sample" in rna_adata.obs.columns:
        sample_col = "sample"
    else:
        raise ValueError("RNA h5ad has neither 'sample_id' nor 'sample' in obs")

    lookup = {}
    for obs_name, sample in zip(rna_adata.obs_names, rna_adata.obs[sample_col]):
        sample = str(sample)
        # Strip the '-SAMPLE' suffix added by anndata index_unique='-'
        if obs_name.endswith("-" + sample):
            raw_bc = obs_name[: -(len(sample) + 1)]
        else:
            raw_bc = obs_name
        lookup[(sample, _strip_gem_well(raw_bc))] = obs_name
    return lookup


# ---------------------------------------------------------------------------
# Gap cell identification
# ---------------------------------------------------------------------------

def identify_gap_cells(mdata, rna_adata, atac_adata):
    """
    Identify cells present in only one modality after QC.

    Handles multi-sample datasets by matching on (sample, normalized_barcode)
    tuples rather than raw string comparison.

    Returns:
      rna_only_obs_names: list of obs_names from rna_adata (for direct subsetting)
      atac_only_obs_names: list of obs_names from atac_adata (for direct subsetting)
      rna_to_mudata: dict mapping rna obs_name → MuData-format barcode
      atac_to_mudata: dict mapping atac obs_name → MuData-format barcode
    """
    # Build (sample, norm_barcode) lookups for each source
    mudata_lookup = _parse_colon_barcodes(mdata.obs_names)
    rna_lookup = _parse_rna_concat_barcodes(rna_adata)
    atac_lookup = _parse_colon_barcodes(atac_adata.obs_names)

    paired_keys = set(mudata_lookup.keys())
    rna_keys = set(rna_lookup.keys())
    atac_keys = set(atac_lookup.keys())

    # Gap cells: in one modality but not the other
    rna_only_keys = rna_keys - atac_keys
    atac_only_keys = atac_keys - rna_keys

    # Map gap cells back to their original obs_names and MuData-format barcodes
    rna_only_obs_names = []
    rna_to_mudata = {}
    for key in sorted(rna_only_keys):
        rna_obs = rna_lookup[key]
        sample, norm_bc = key
        mudata_bc = f"{sample}:{norm_bc}" if sample else norm_bc
        rna_only_obs_names.append(rna_obs)
        rna_to_mudata[rna_obs] = mudata_bc

    atac_only_obs_names = []
    atac_to_mudata = {}
    for key in sorted(atac_only_keys):
        atac_obs = atac_lookup[key]
        sample, norm_bc = key
        mudata_bc = f"{sample}:{norm_bc}" if sample else norm_bc
        atac_only_obs_names.append(atac_obs)
        atac_to_mudata[atac_obs] = mudata_bc

    print(f"  Paired cells (in MuData): {len(paired_keys)}")
    print(f"  Total RNA cells (post-QC): {len(rna_keys)}")
    print(f"  Total ATAC cells (post-QC): {len(atac_keys)}")
    print(f"  RNA-only gap cells: {len(rna_only_obs_names)}")
    print(f"  ATAC-only gap cells: {len(atac_only_obs_names)}")

    n_samples = len({k[0] for k in rna_keys | atac_keys})
    if n_samples > 1:
        print(f"  Samples detected: {n_samples}")

    return rna_only_obs_names, atac_only_obs_names, rna_to_mudata, atac_to_mudata


def _populate_gap_obs(gap_obs, paired_obs, src_obs, mod_prefix):
    """Fill gap_obs with columns from paired_obs, resolving src_obs via mod-prefix strip.

    paired_obs has MuData-joint columns like 'rna:celltypist_prediction'.
    src_obs (raw h5ad .obs) has bare columns like 'celltypist_prediction'.
    For each paired column, strip the modality prefix and look up the bare name;
    fall back to exact-match if no prefix was present. Missing → NaN.
    """
    prefix = f"{mod_prefix}:"
    for col in paired_obs.columns:
        bare = col[len(prefix):] if col.startswith(prefix) else col
        if bare in src_obs.columns:
            gap_obs[col] = src_obs[bare].values
        elif col in src_obs.columns:
            gap_obs[col] = src_obs[col].values
        else:
            gap_obs[col] = np.nan


def build_expanded_mudata(mdata, rna_adata, atac_adata,
                          rna_only_obs_names, atac_only_obs_names,
                          rna_to_mudata, atac_to_mudata,
                          cell_type_key="celltypist_prediction"):
    """
    Build a MuData that includes paired cells + gap cells.

    Gap cells are looked up by their original obs_names in the source h5ads,
    then assigned MuData-format barcodes via the mapping dicts.
    """
    from scipy.sparse import csr_matrix, vstack

    # Start with the existing paired data
    rna_X = mdata.mod["rna"].X.copy()
    atac_X = mdata.mod["atac"].X.copy()
    rna_var = mdata.mod["rna"].var.copy()
    atac_var = mdata.mod["atac"].var.copy()
    paired_obs = mdata.obs.copy()

    if issparse(rna_X):
        rna_X_dense = False
    else:
        rna_X_dense = True

    n_genes = rna_X.shape[1]
    n_peaks = atac_X.shape[1]

    # Get feature names for alignment
    rna_genes = mdata.mod["rna"].var_names.values
    atac_peaks = mdata.mod["atac"].var_names.values

    all_barcodes = list(mdata.obs_names)
    all_modality = ["paired"] * len(all_barcodes)
    all_rna_rows = [rna_X]
    all_atac_rows = [atac_X]

    obs_records = [paired_obs]

    # Add RNA-only gap cells (need ATAC imputed)
    if len(rna_only_obs_names) > 0:
        print(f"  Adding {len(rna_only_obs_names)} RNA-only cells...")
        available = [bc for bc in rna_only_obs_names if bc in rna_adata.obs_names]
        if len(available) < len(rna_only_obs_names):
            print(f"    WARNING: {len(rna_only_obs_names) - len(available)} RNA-only barcodes not found in RNA h5ad")

        if len(available) > 0:
            rna_sub = rna_adata[available, :]
            mudata_bcs = [rna_to_mudata[bc] for bc in available]

            # Align genes to MuData's gene set
            common_genes = [g for g in rna_genes if g in rna_sub.var_names]
            if len(common_genes) < n_genes:
                rna_aligned = np.zeros((len(available), n_genes), dtype=np.float32)
                sub_gene_idx = {g: i for i, g in enumerate(rna_sub.var_names)}
                for j, gene in enumerate(rna_genes):
                    if gene in sub_gene_idx:
                        col = rna_sub.X[:, sub_gene_idx[gene]]
                        rna_aligned[:, j] = col.toarray().flatten() if issparse(col) else col.flatten()
            else:
                rna_aligned = rna_sub[:, rna_genes].X
                if issparse(rna_aligned):
                    rna_aligned = rna_aligned.toarray()

            atac_zeros = np.zeros((len(available), n_peaks), dtype=np.float32)

            if not rna_X_dense:
                all_rna_rows.append(csr_matrix(rna_aligned))
                all_atac_rows.append(csr_matrix(atac_zeros))
            else:
                all_rna_rows.append(rna_aligned)
                all_atac_rows.append(atac_zeros)

            all_barcodes.extend(mudata_bcs)
            all_modality.extend(["expression"] * len(available))

            # Propagate obs columns from source h5ad
            gap_obs = pd.DataFrame(index=mudata_bcs)
            src_obs = rna_sub.obs.copy()
            src_obs.index = mudata_bcs
            _populate_gap_obs(gap_obs, paired_obs, src_obs, mod_prefix="rna")
            obs_records.append(gap_obs)

    # Add ATAC-only gap cells (need RNA imputed)
    if len(atac_only_obs_names) > 0:
        print(f"  Adding {len(atac_only_obs_names)} ATAC-only cells...")
        available = [bc for bc in atac_only_obs_names if bc in atac_adata.obs_names]
        if len(available) < len(atac_only_obs_names):
            print(f"    WARNING: {len(atac_only_obs_names) - len(available)} ATAC-only barcodes not found in ATAC h5ad")

        if len(available) > 0:
            atac_sub = atac_adata[available, :]
            mudata_bcs = [atac_to_mudata[bc] for bc in available]

            common_peaks = [p for p in atac_peaks if p in atac_sub.var_names]
            if len(common_peaks) < n_peaks:
                atac_aligned = np.zeros((len(available), n_peaks), dtype=np.float32)
                sub_peak_idx = {p: i for i, p in enumerate(atac_sub.var_names)}
                for j, peak in enumerate(atac_peaks):
                    if peak in sub_peak_idx:
                        col = atac_sub.X[:, sub_peak_idx[peak]]
                        atac_aligned[:, j] = col.toarray().flatten() if issparse(col) else col.flatten()
            else:
                atac_aligned = atac_sub[:, atac_peaks].X
                if issparse(atac_aligned):
                    atac_aligned = atac_aligned.toarray()

            rna_zeros = np.zeros((len(available), n_genes), dtype=np.float32)

            if not rna_X_dense:
                all_rna_rows.append(csr_matrix(rna_zeros))
                all_atac_rows.append(csr_matrix(atac_aligned))
            else:
                all_rna_rows.append(rna_zeros)
                all_atac_rows.append(atac_aligned)

            all_barcodes.extend(mudata_bcs)
            all_modality.extend(["accessibility"] * len(available))

            gap_obs = pd.DataFrame(index=mudata_bcs)
            src_obs = atac_sub.obs.copy()
            src_obs.index = mudata_bcs
            _populate_gap_obs(gap_obs, paired_obs, src_obs, mod_prefix="atac")
            obs_records.append(gap_obs)

    # Stack matrices
    if not rna_X_dense:
        from scipy.sparse import vstack as sparse_vstack
        final_rna = sparse_vstack(all_rna_rows)
        final_atac = sparse_vstack(all_atac_rows)
    else:
        final_rna = np.vstack(all_rna_rows)
        final_atac = np.vstack(all_atac_rows)

    # Build new MuData
    all_obs = pd.concat(obs_records)
    all_obs.index = all_barcodes
    all_obs["modality"] = all_modality

    rna_ad = ad.AnnData(X=final_rna, obs=all_obs.copy(), var=rna_var)
    rna_ad.obs_names = all_barcodes

    atac_ad = ad.AnnData(X=final_atac, obs=all_obs.copy(), var=atac_var)
    atac_ad.obs_names = all_barcodes

    expanded = mu.MuData({"rna": rna_ad, "atac": atac_ad})
    expanded.obs["modality"] = all_modality

    # Guarantee an unprefixed cell-type column on expanded.obs so downstream
    # validators (multivi_validate.py) can find it directly regardless of
    # MuData's join-level prefixing behavior. Source of truth for cell_type_key,
    # in priority order:
    #   1. all_obs[f"rna:{cell_type_key}"] — populated for paired cells from
    #      paired_obs (mdata.obs) and for gap cells by _populate_gap_obs.
    #   2. all_obs[cell_type_key]          — unprefixed variant (some paths use this).
    #   3. mdata.mod["rna"].obs[cell_type_key] by barcode match — for paired cells
    #      when the joint obs lost the column.
    rna_ct_col = f"rna:{cell_type_key}"
    ct_values = None
    if rna_ct_col in all_obs.columns:
        ct_values = all_obs[rna_ct_col].values
    elif cell_type_key in all_obs.columns:
        ct_values = all_obs[cell_type_key].values
    elif (hasattr(mdata, "mod") and "rna" in mdata.mod
            and cell_type_key in mdata.mod["rna"].obs.columns):
        # Reindex mod-level cell types to our all_barcodes (paired subset matches;
        # gap cells end up NaN and can be filled downstream if needed).
        src = mdata.mod["rna"].obs[cell_type_key]
        ct_values = src.reindex(all_barcodes).values
    if ct_values is not None:
        expanded.obs[cell_type_key] = ct_values
        expanded.mod["rna"].obs[cell_type_key] = ct_values

    # Propagate cell_type_source (provenance) if present in paired_obs. Gap cells
    # inherit the same provenance since the annotation tool is a pipeline-wide choice.
    src_col = f"rna:cell_type_source"
    if src_col in all_obs.columns:
        expanded.obs["cell_type_source"] = all_obs[src_col].values
    elif "cell_type_source" in all_obs.columns:
        expanded.obs["cell_type_source"] = all_obs["cell_type_source"].values

    # True vs MultiVI-imputed distinction: paired cells have observed cell types
    # on both modalities; gap cells have observed cell type on one side only and
    # get a MultiVI-imputed representation on the other.
    modality_arr = np.array(all_modality)
    expanded.obs["cell_type_imputed"] = modality_arr != "paired"

    print(f"  Expanded MuData: {expanded.n_obs} cells "
          f"({sum(m == 'paired' for m in all_modality)} paired, "
          f"{sum(m == 'expression' for m in all_modality)} RNA-only, "
          f"{sum(m == 'accessibility' for m in all_modality)} ATAC-only)")

    return expanded


# ---------------------------------------------------------------------------
# Confidence scoring
# ---------------------------------------------------------------------------

def _welford_mean_var(sampler, n_samples):
    """
    Streaming mean/variance (Welford) across n_samples draws from `sampler()`.

    Returns (mean, var) with the same shape as one draw and population variance
    (ddof=0, matching np.var default). Memory: 3 arrays of the draw shape —
    stacking 25 draws of a (n_cells × n_peaks) ATAC matrix would otherwise
    blow past the per-task RAM budget on BD-scale datasets.
    """
    n = 0
    mean = None
    M2 = None
    for _ in range(n_samples):
        x = sampler().astype(np.float32, copy=False)
        n += 1
        if mean is None:
            mean = np.zeros_like(x)
            M2 = np.zeros_like(x)
        delta = x - mean
        mean += delta / n
        M2 += delta * (x - mean)
    var = M2 / n
    return mean, var


def compute_confidence(mvi, mdata, n_samples=25):
    """
    Compute per-cell, per-feature confidence via posterior sampling.

    Draws n_samples from the posterior and computes variance. Uses Welford
    streaming so peak memory is ~3 × single-sample size rather than
    n_samples × single-sample. Confidence = 1 - CV² (clamped to [0, 1]).
    """
    print(f"  Drawing {n_samples} posterior samples for confidence scoring (Welford streaming)...")

    rna_mean, rna_var = _welford_mean_var(
        lambda: mvi.get_normalized_expression(return_numpy=True), n_samples
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        rna_cv2 = np.where(rna_mean > 0, rna_var / (rna_mean**2), 0)
    rna_confidence = np.clip(1 - rna_cv2, 0, 1)

    atac_mean, atac_var = _welford_mean_var(
        lambda: mvi.get_normalized_accessibility(return_numpy=True), n_samples
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        atac_cv2 = np.where(atac_mean > 0, atac_var / (atac_mean**2), 0)
    atac_confidence = np.clip(1 - atac_cv2, 0, 1)

    return rna_confidence, atac_confidence, rna_mean, atac_mean


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_confidence_distributions(mdata, rna_conf, atac_conf, output_dir):
    """Histograms of per-cell mean confidence, split by modality."""
    modality = np.array(mdata.obs["modality"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # RNA confidence for ATAC-only cells (these are the imputed ones)
    atac_only = modality == "accessibility"
    if atac_only.sum() > 0:
        mean_conf = np.mean(rna_conf[atac_only], axis=1)
        axes[0].hist(mean_conf, bins=50, color="steelblue", alpha=0.7, edgecolor="white")
        axes[0].axvline(x=np.median(mean_conf), color="red", ls="--", label=f"Median: {np.median(mean_conf):.3f}")
        axes[0].set_title(f"RNA Imputation Confidence\n(ATAC-only cells, n={atac_only.sum()})")
        axes[0].set_xlabel("Mean Confidence per Cell")
        axes[0].set_ylabel("Count")
        axes[0].legend()

    # ATAC confidence for RNA-only cells
    rna_only = modality == "expression"
    if rna_only.sum() > 0:
        mean_conf = np.mean(atac_conf[rna_only], axis=1)
        axes[1].hist(mean_conf, bins=50, color="coral", alpha=0.7, edgecolor="white")
        axes[1].axvline(x=np.median(mean_conf), color="red", ls="--", label=f"Median: {np.median(mean_conf):.3f}")
        axes[1].set_title(f"ATAC Imputation Confidence\n(RNA-only cells, n={rna_only.sum()})")
        axes[1].set_xlabel("Mean Confidence per Cell")
        axes[1].set_ylabel("Count")
        axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gap_fill_confidence_distributions.pdf"), dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Load data
    print("Loading paired MuData...")
    mdata = mu.read_h5mu(args.mudata)

    print("Loading pre-intersection RNA h5ad...")
    rna_adata = sc.read_h5ad(args.rna_h5ad)

    print("Loading pre-intersection ATAC h5ad...")
    atac_adata = sc.read_h5ad(args.atac_h5ad)

    # Step 2: Identify gap cells
    print("\nIdentifying gap cells...")
    rna_only_obs, atac_only_obs, rna_to_mudata, atac_to_mudata = identify_gap_cells(
        mdata, rna_adata, atac_adata
    )

    if len(rna_only_obs) == 0 and len(atac_only_obs) == 0:
        print("\nNo gap cells found — all cells are paired. Nothing to impute.")
        mdata.obs["rna_is_imputed"] = False
        mdata.obs["atac_is_imputed"] = False
        mdata.obs["cell_type_imputed"] = False
        # Guarantee an unprefixed cell-type column for downstream validators.
        rna_ct_col = f"rna:{args.cell_type_key}"
        if args.cell_type_key not in mdata.obs.columns:
            if rna_ct_col in mdata.obs.columns:
                mdata.obs[args.cell_type_key] = mdata.obs[rna_ct_col].values
            elif args.cell_type_key in mdata.mod["rna"].obs.columns:
                mdata.obs[args.cell_type_key] = mdata.mod["rna"].obs[args.cell_type_key].reindex(mdata.obs_names).values
        from h5ad_compat import sanitize_mudata
        sanitize_mudata(mdata, os.path.join(args.output_dir, "multivi_gap_filled.h5mu"))
        _sanitize_mudata_obs_for_h5py(mdata)
        mdata.write_h5mu(os.path.join(args.output_dir, "multivi_gap_filled.h5mu"))

        summary = {"n_paired": mdata.n_obs, "n_rna_only": 0, "n_atac_only": 0, "gap_cells_found": False}
        with open(os.path.join(args.output_dir, "gap_fill_stats.json"), "w") as f:
            json.dump(summary, f, indent=2)
        return

    # Step 3: Build expanded MuData
    print("\nBuilding expanded MuData with gap cells...")
    expanded = build_expanded_mudata(
        mdata, rna_adata, atac_adata,
        rna_only_obs, atac_only_obs,
        rna_to_mudata, atac_to_mudata,
        cell_type_key=args.cell_type_key,
    )

    # Step 4: Train MultiVI on expanded data
    print(f"\nSetting up and training MultiVI ({args.n_epochs} epochs)...")
    batch_col = args.batch_key
    if batch_col not in expanded.mod["rna"].obs.columns:
        expanded.mod["rna"].obs[batch_col] = "batch0"

    # Safety net: if any gap cells still have NaN batch_col (source h5ad
    # lacked the column), fill with the first valid value.
    # Also cast to str — scvi 1.4 _make_column_categorical fails when the
    # column is already Categorical with Index-typed categories.
    for mod_key in ["rna", "atac"]:
        col = expanded.mod[mod_key].obs[batch_col]
        n_na = col.isna().sum()
        if n_na > 0:
            fill_val = col.dropna().iloc[0]
            expanded.mod[mod_key].obs[batch_col] = col.fillna(fill_val).astype(str)
            print(f"  WARNING: Filled {n_na} NaN {batch_col} values in {mod_key} with '{fill_val}'")
        else:
            expanded.mod[mod_key].obs[batch_col] = col.astype(str)

    MULTIVI.setup_mudata(
        expanded,
        rna_layer=None,
        atac_layer=None,
        protein_layer=None,
        batch_key=batch_col,
        modalities={
            "rna_layer": "rna",
            "atac_layer": "atac",
            "protein_layer": None,
            "batch_key": "rna",
        },
    )

    mvi = MULTIVI(
        expanded,
        modality_weights=args.modality_weights,
        modality_penalty=args.modality_penalty,
        n_latent=args.n_latent,
        region_factors=True,
        gene_likelihood="zinb",
        dispersion="gene",
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    )

    # scvi 1.4.2 hardcodes early_stopping_patience=50 in MULTIVI.train's TrainRunner
    # call, so passing it as a kwarg collides at the call site. Disable scvi's default
    # and supply our own LoudEarlyStopping with patience=10.
    mvi.train(
        max_epochs=args.n_epochs,
        adversarial_mixing=True,
        early_stopping=False,
        check_val_every_n_epoch=1,
        callbacks=[LoudEarlyStopping(
            monitor="reconstruction_loss_validation",
            patience=10,
            mode="min",
            warmup_epochs=0,
        )],
    )

    # Step 5: Get latent representation and imputations
    print("\nComputing latent representation and imputations...")
    expanded.obsm["X_MultiVI"] = mvi.get_latent_representation(modality="joint")

    # Imputed values for ALL cells
    expr_imputed = mvi.get_normalized_expression(return_numpy=True)
    acc_imputed = mvi.get_normalized_accessibility(return_numpy=True)

    # Store imputed layers
    expanded.mod["rna"].layers["multivi_imputed"] = expr_imputed
    expanded.mod["atac"].layers["multivi_imputed"] = acc_imputed

    # Set imputation flags
    modality = np.array(expanded.obs["modality"])
    expanded.obs["rna_is_imputed"] = modality == "accessibility"   # RNA was imputed for ATAC-only cells
    expanded.obs["atac_is_imputed"] = modality == "expression"     # ATAC was imputed for RNA-only cells

    # Step 6: Confidence scoring
    print("\nComputing confidence scores via posterior sampling...")
    rna_conf, atac_conf, _, _ = compute_confidence(mvi, expanded, args.n_posterior_samples)

    expanded.mod["rna"].layers["imputation_confidence"] = rna_conf
    expanded.mod["atac"].layers["imputation_confidence"] = atac_conf

    # Step 7: Save
    print("\nSaving gap-filled MuData...")
    out_path = os.path.join(args.output_dir, "multivi_gap_filled.h5mu")
    try:
        from h5ad_compat import sanitize_mudata
        sanitize_mudata(expanded, out_path)
    except ImportError:
        pass
    _sanitize_mudata_obs_for_h5py(expanded)
    expanded.write_h5mu(out_path)

    # Save model
    model_dir = os.path.join(args.output_dir, "multivi_gap_fill_model")
    mvi.save(model_dir, overwrite=True)

    # Step 8: Summary statistics
    n_rna_only = int((modality == "expression").sum())
    n_atac_only = int((modality == "accessibility").sum())
    n_paired = int((modality == "paired").sum())

    rna_conf_imputed = rna_conf[modality == "accessibility"] if n_atac_only > 0 else np.array([])
    atac_conf_imputed = atac_conf[modality == "expression"] if n_rna_only > 0 else np.array([])

    summary = {
        "n_paired": n_paired,
        "n_rna_only": n_rna_only,
        "n_atac_only": n_atac_only,
        "n_total": int(expanded.n_obs),
        "gap_cells_found": True,
        "rna_imputation_confidence": {
            "mean": float(np.mean(rna_conf_imputed)) if len(rna_conf_imputed) > 0 else None,
            "median": float(np.median(rna_conf_imputed)) if len(rna_conf_imputed) > 0 else None,
        },
        "atac_imputation_confidence": {
            "mean": float(np.mean(atac_conf_imputed)) if len(atac_conf_imputed) > 0 else None,
            "median": float(np.median(atac_conf_imputed)) if len(atac_conf_imputed) > 0 else None,
        },
        "min_confidence_threshold": args.min_confidence,
        "n_posterior_samples": args.n_posterior_samples,
    }
    with open(os.path.join(args.output_dir, "gap_fill_stats.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # Step 9: Plots
    print("\nGenerating plots...")
    plot_confidence_distributions(expanded, rna_conf, atac_conf, args.output_dir)

    print(f"\nGap-filling complete. Results in {args.output_dir}/")
    print(f"  Total cells: {expanded.n_obs} ({n_paired} paired + {n_rna_only} RNA-only + {n_atac_only} ATAC-only)")


if __name__ == "__main__":
    main()
