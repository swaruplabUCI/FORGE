#!/usr/bin/env python
"""
MultiVI Masking Sweep — Cell-level artificial unpairing benchmark.

Following Ashuach et al. (Nature Methods 2023): randomly unpair a fraction of
paired multiome cells, retrain MultiVI, then measure how well the model imputes
the held-out modality.  Sweeps over unpairing fractions {25%, 50%, 75%} with
multiple random seeds.

Metrics (matching the paper):
  RNA  — Spearman r (raw), Spearman r (smoothed via 50-NN)
  ATAC — AUPRC (raw), Spearman r (smoothed via 50-NN)
  Integration — LISI, ARI, silhouette
"""

import argparse
import json
import os
import warnings

import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, silhouette_score
from sklearn.neighbors import NearestNeighbors
from scvi.model import MULTIVI

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def parse_args():
    p = argparse.ArgumentParser(description="MultiVI masking sweep benchmark")
    p.add_argument("--mudata", required=True, help="Pre-MultiVI integrated MuData (.h5mu)")
    p.add_argument("--output_dir", default="masking_sweep", help="Output directory")
    p.add_argument("--fractions", default="0.25,0.50,0.75", help="Comma-separated unpairing fractions")
    p.add_argument("--seeds", default="42,123,456", help="Comma-separated random seeds")
    p.add_argument("--n_epochs", type=int, default=200, help="Training epochs per run")
    p.add_argument("--batch_key", default="sample_id", help="Batch key in obs")
    p.add_argument("--cell_type_key", default="celltypist_prediction", help="Cell type key for stratification")
    p.add_argument("--n_latent", type=int, default=20, help="Latent space dimensionality")
    p.add_argument("--modality_weights", default="equal")
    p.add_argument("--modality_penalty", default="Jeffreys")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Artificial unpairing
# ---------------------------------------------------------------------------

def artificially_unpair(mdata, fraction, seed):
    """
    Take a fully paired MuData and artificially unpair a fraction of cells.

    For each unpairing cell:
      - The cell appears twice: once as expression-only, once as accessibility-only.
      - The missing modality is zeroed out.

    Returns: (modified_mdata, held_out_rna, held_out_atac, unpaired_mask)
      held_out_rna:  dict {barcode: dense array} for RNA ground truth of ATAC-only cells
      held_out_atac: dict {barcode: dense array} for ATAC ground truth of RNA-only cells
    """
    import copy
    rng = np.random.RandomState(seed)

    n_cells = mdata.n_obs
    n_unpair = int(n_cells * fraction)
    unpair_idx = rng.choice(n_cells, n_unpair, replace=False)
    unpair_mask = np.zeros(n_cells, dtype=bool)
    unpair_mask[unpair_idx] = True

    barcodes = mdata.obs_names.values
    rna_X = mdata.mod["rna"].X.copy()
    atac_X = mdata.mod["atac"].X.copy()

    if issparse(rna_X):
        rna_X = rna_X.toarray()
    if issparse(atac_X):
        atac_X = atac_X.toarray()

    # Store ground truth for unpairing cells
    held_out_rna = {}
    held_out_atac = {}
    for idx in unpair_idx:
        bc = barcodes[idx]
        held_out_rna[bc] = rna_X[idx].copy()
        held_out_atac[bc] = atac_X[idx].copy()

    # Build the new dataset: paired cells + expression-only + accessibility-only
    # Paired cells (not unpairing)
    paired_idx = np.where(~unpair_mask)[0]

    # For unpairing cells, create two copies each
    # Expression-only copies: zero out ATAC
    rna_only_rna = rna_X[unpair_idx]
    rna_only_atac = np.zeros_like(atac_X[unpair_idx])

    # Accessibility-only copies: zero out RNA
    atac_only_rna = np.zeros_like(rna_X[unpair_idx])
    atac_only_atac = atac_X[unpair_idx]

    # Concatenate
    new_rna_X = np.vstack([rna_X[paired_idx], rna_only_rna, atac_only_rna])
    new_atac_X = np.vstack([atac_X[paired_idx], rna_only_atac, atac_only_atac])

    # Build new barcodes
    paired_bcs = barcodes[paired_idx]
    rna_only_bcs = np.array([f"{bc}_rna_only" for bc in barcodes[unpair_idx]])
    atac_only_bcs = np.array([f"{bc}_atac_only" for bc in barcodes[unpair_idx]])
    new_barcodes = np.concatenate([paired_bcs, rna_only_bcs, atac_only_bcs])

    # Build modality indicator
    modality = (
        ["paired"] * len(paired_idx)
        + ["expression"] * n_unpair
        + ["accessibility"] * n_unpair
    )

    # Propagate obs columns
    paired_obs = mdata.obs.iloc[paired_idx].copy()
    rna_only_obs = mdata.obs.iloc[unpair_idx].copy()
    atac_only_obs = mdata.obs.iloc[unpair_idx].copy()

    paired_obs.index = paired_bcs
    rna_only_obs.index = rna_only_bcs
    atac_only_obs.index = atac_only_bcs

    new_obs = pd.concat([paired_obs, rna_only_obs, atac_only_obs])
    new_obs["modality"] = modality

    # Build new AnnData objects for each modality
    import anndata as ad

    rna_adata = ad.AnnData(
        X=new_rna_X,
        obs=new_obs.copy(),
        var=mdata.mod["rna"].var.copy(),
    )
    rna_adata.obs_names = new_barcodes.copy()

    atac_adata = ad.AnnData(
        X=new_atac_X,
        obs=new_obs.copy(),
        var=mdata.mod["atac"].var.copy(),
    )
    atac_adata.obs_names = new_barcodes.copy()

    # Build MuData
    new_mdata = mu.MuData({"rna": rna_adata, "atac": atac_adata})
    new_mdata.obs["modality"] = modality

    # Propagate batch key to rna modality obs if needed
    if "sample_id" in mdata.mod["rna"].obs.columns:
        new_mdata.mod["rna"].obs["sample_id"] = new_obs["sample_id"].values if "sample_id" in new_obs.columns else "batch0"

    return new_mdata, held_out_rna, held_out_atac, unpair_mask


# ---------------------------------------------------------------------------
# Smoothed ground truth (50-NN average)
# ---------------------------------------------------------------------------

def compute_smoothed_ground_truth(mdata, latent_key="X_MultiVI", k=50):
    """
    For each cell, average the observed values of its k nearest neighbors
    in the joint latent space.  Returns smoothed RNA and ATAC matrices.
    """
    Z = mdata.obsm[latent_key]
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(Z)
    _, indices = nn.kneighbors(Z)
    # Exclude self (index 0)
    neighbor_idx = indices[:, 1:]

    rna_X = mdata.mod["rna"].X
    atac_X = mdata.mod["atac"].X
    if issparse(rna_X):
        rna_X = rna_X.toarray()
    if issparse(atac_X):
        atac_X = atac_X.toarray()

    smoothed_rna = np.zeros_like(rna_X)
    smoothed_atac = np.zeros_like(atac_X)

    for i in range(len(Z)):
        smoothed_rna[i] = rna_X[neighbor_idx[i]].mean(axis=0)
        smoothed_atac[i] = atac_X[neighbor_idx[i]].mean(axis=0)

    return smoothed_rna, smoothed_atac


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_imputation_metrics(
    imputed, ground_truth, modality, cell_barcodes, cell_types=None
):
    """
    Compute per-cell imputation metrics.

    For RNA: Spearman r per cell (across genes)
    For ATAC: AUPRC per cell (binary ground truth vs predicted probabilities)

    Returns a DataFrame with one row per cell.
    """
    records = []
    for i, bc in enumerate(cell_barcodes):
        gt = ground_truth[i]
        pred = imputed[i]

        rec = {"barcode": bc, "modality": modality}
        if cell_types is not None:
            rec["cell_type"] = cell_types[i]

        if modality == "rna":
            # Spearman correlation across genes
            if gt.std() > 0 and pred.std() > 0:
                rec["spearman_r"] = spearmanr(gt, pred).correlation
            else:
                rec["spearman_r"] = np.nan
        elif modality == "atac":
            # AUPRC: binarize ground truth at > 0
            gt_binary = (gt > 0).astype(int)
            if gt_binary.sum() > 0 and gt_binary.sum() < len(gt_binary):
                rec["auprc"] = average_precision_score(gt_binary, pred)
            else:
                rec["auprc"] = np.nan

        records.append(rec)

    return pd.DataFrame(records)


def compute_smoothed_metrics(imputed, smoothed_gt, modality, cell_barcodes):
    """Spearman r against smoothed ground truth (50-NN average)."""
    records = []
    for i, bc in enumerate(cell_barcodes):
        gt = smoothed_gt[i]
        pred = imputed[i]
        rec = {"barcode": bc, "modality": modality}
        if gt.std() > 0 and pred.std() > 0:
            rec["spearman_r_smoothed"] = spearmanr(gt, pred).correlation
        else:
            rec["spearman_r_smoothed"] = np.nan
        records.append(rec)
    return pd.DataFrame(records)


def compute_lisi(latent, labels, perplexity=30):
    """
    Compute the Local Inverse Simpson's Index for modality mixing.
    Simplified implementation using kNN.
    """
    from sklearn.neighbors import NearestNeighbors

    k = int(3 * perplexity)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(latent)
    _, indices = nn.kneighbors(latent)

    unique_labels = np.unique(labels)
    label_to_int = {l: i for i, l in enumerate(unique_labels)}
    int_labels = np.array([label_to_int[l] for l in labels])

    lisi_scores = []
    for i in range(len(latent)):
        neighbor_labels = int_labels[indices[i]]
        freqs = np.bincount(neighbor_labels, minlength=len(unique_labels)) / k
        simpson = np.sum(freqs**2)
        lisi_scores.append(1.0 / simpson if simpson > 0 else 1.0)

    return np.array(lisi_scores)


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_degradation_curves(all_results, output_dir):
    """Plot metric vs unpairing fraction — the key masking sweep figure."""
    df = all_results.copy()

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: RNA Spearman r (raw)
    rna_df = df[df["modality"] == "rna"].dropna(subset=["spearman_r"])
    if len(rna_df) > 0:
        sns.lineplot(
            data=rna_df, x="fraction", y="spearman_r",
            ci=95, marker="o", ax=axes[0], color="steelblue",
        )
        axes[0].set_title("RNA Imputation (raw)")
        axes[0].set_xlabel("Unpairing Fraction")
        axes[0].set_ylabel("Spearman r")
        axes[0].set_ylim(0, 1)
        axes[0].axhline(y=0.57, color="gray", ls="--", alpha=0.5, label="Paper (0.57)")
        axes[0].legend()

    # Panel 2: ATAC AUPRC (raw)
    atac_df = df[df["modality"] == "atac"].dropna(subset=["auprc"])
    if len(atac_df) > 0:
        sns.lineplot(
            data=atac_df, x="fraction", y="auprc",
            ci=95, marker="o", ax=axes[1], color="coral",
        )
        axes[1].set_title("ATAC Imputation (raw)")
        axes[1].set_xlabel("Unpairing Fraction")
        axes[1].set_ylabel("AUPRC")
        axes[1].set_ylim(0, 1)
        axes[1].axhline(y=0.41, color="gray", ls="--", alpha=0.5, label="Paper (0.41)")
        axes[1].legend()

    # Panel 3: Smoothed Spearman r (both modalities)
    smoothed_cols = [c for c in df.columns if "smoothed" in c]
    if smoothed_cols:
        for mod, color, label in [("rna", "steelblue", "RNA"), ("atac", "coral", "ATAC")]:
            mod_df = df[df["modality"] == mod].dropna(subset=["spearman_r_smoothed"])
            if len(mod_df) > 0:
                sns.lineplot(
                    data=mod_df, x="fraction", y="spearman_r_smoothed",
                    ci=95, marker="o", ax=axes[2], color=color, label=label,
                )
        axes[2].set_title("Imputation (smoothed, 50-NN)")
        axes[2].set_xlabel("Unpairing Fraction")
        axes[2].set_ylabel("Spearman r (smoothed)")
        axes[2].set_ylim(0, 1)
        axes[2].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "masking_degradation_curves.pdf"), dpi=150)
    plt.close()


def plot_celltype_breakdown(all_results, output_dir):
    """Per-cell-type metric breakdown at each fraction."""
    df = all_results.copy()

    for mod, metric in [("rna", "spearman_r"), ("atac", "auprc")]:
        mod_df = df[(df["modality"] == mod) & df["cell_type"].notna()].dropna(subset=[metric])
        if len(mod_df) == 0:
            continue

        fig, ax = plt.subplots(figsize=(12, 6))
        sns.boxplot(
            data=mod_df, x="cell_type", y=metric, hue="fraction",
            ax=ax, palette="viridis",
        )
        ax.set_title(f"{mod.upper()} Imputation by Cell Type")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(
            os.path.join(output_dir, f"masking_celltype_{mod}.pdf"), dpi=150
        )
        plt.close()


def plot_integration_quality(integration_metrics, output_dir):
    """Plot LISI and silhouette across fractions."""
    df = pd.DataFrame(integration_metrics)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sns.barplot(data=df, x="fraction", y="lisi_mean", hue="seed", ax=axes[0])
    axes[0].set_title("Modality LISI (higher = better mixing)")
    axes[0].set_ylabel("Mean LISI")

    sns.barplot(data=df, x="fraction", y="silhouette", hue="seed", ax=axes[1])
    axes[1].set_title("Cell Type Silhouette (higher = better separation)")
    axes[1].set_ylabel("Silhouette Score")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "masking_integration_quality.pdf"), dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    fractions = [float(f) for f in args.fractions.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]

    print(f"Loading MuData from {args.mudata}...")
    mdata_orig = mu.read_h5mu(args.mudata)
    print(f"  Cells: {mdata_orig.n_obs}")
    print(f"  RNA features: {mdata_orig.mod['rna'].n_vars}")
    print(f"  ATAC features: {mdata_orig.mod['atac'].n_vars}")

    all_results = []
    integration_metrics = []

    for frac in fractions:
        for seed in seeds:
            print(f"\n{'='*60}")
            print(f"  Fraction: {frac}, Seed: {seed}")
            print(f"{'='*60}")

            # Step 1: Artificially unpair
            print("  Artificially unpairing cells...")
            new_mdata, held_out_rna, held_out_atac, _ = artificially_unpair(
                mdata_orig, frac, seed
            )
            n_paired = (np.array(new_mdata.obs["modality"]) == "paired").sum()
            n_rna_only = (np.array(new_mdata.obs["modality"]) == "expression").sum()
            n_atac_only = (np.array(new_mdata.obs["modality"]) == "accessibility").sum()
            print(f"  Paired: {n_paired}, RNA-only: {n_rna_only}, ATAC-only: {n_atac_only}")

            # Step 2: Register and train MultiVI
            print("  Setting up MultiVI...")
            batch_col = args.batch_key
            if batch_col not in new_mdata.mod["rna"].obs.columns:
                new_mdata.mod["rna"].obs[batch_col] = "batch0"

            MULTIVI.setup_mudata(
                new_mdata,
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
                new_mdata,
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

            print(f"  Training MultiVI for {args.n_epochs} epochs...")
            mvi.train(max_epochs=args.n_epochs, adversarial_mixing=True)

            # Step 3: Get latent representation and imputed values
            print("  Computing latent representation and imputations...")
            Z = mvi.get_latent_representation(modality="joint")
            new_mdata.obsm["X_MultiVI"] = Z

            expr_imputed = mvi.get_normalized_expression(return_numpy=True)
            acc_imputed = mvi.get_normalized_accessibility(return_numpy=True)

            # Step 4: Evaluate RNA imputation (ATAC-only cells -> imputed RNA)
            print("  Evaluating RNA imputation on ATAC-only cells...")
            atac_only_mask = np.array(new_mdata.obs["modality"]) == "accessibility"
            atac_only_bcs = new_mdata.obs_names[atac_only_mask]
            # Map back to original barcodes
            orig_bcs = [bc.replace("_atac_only", "") for bc in atac_only_bcs]

            rna_gt = np.array([held_out_rna[bc] for bc in orig_bcs])
            rna_imp = expr_imputed[atac_only_mask]

            # Cell types for stratification
            ct_key = args.cell_type_key
            cell_types = None
            for candidate in [ct_key, f"rna:{ct_key}", "celltypist_prediction", "rna:celltypist_prediction"]:
                if candidate in new_mdata.obs.columns:
                    cell_types = new_mdata.obs[candidate].values[atac_only_mask]
                    break

            rna_metrics = compute_imputation_metrics(
                rna_imp, rna_gt, "rna", atac_only_bcs.values, cell_types
            )
            rna_metrics["fraction"] = frac
            rna_metrics["seed"] = seed

            # Step 5: Evaluate ATAC imputation (RNA-only cells -> imputed ATAC)
            print("  Evaluating ATAC imputation on RNA-only cells...")
            rna_only_mask = np.array(new_mdata.obs["modality"]) == "expression"
            rna_only_bcs = new_mdata.obs_names[rna_only_mask]
            orig_bcs_atac = [bc.replace("_rna_only", "") for bc in rna_only_bcs]

            atac_gt = np.array([held_out_atac[bc] for bc in orig_bcs_atac])
            atac_imp = acc_imputed[rna_only_mask]

            cell_types_atac = None
            for candidate in [ct_key, f"rna:{ct_key}", "celltypist_prediction", "rna:celltypist_prediction"]:
                if candidate in new_mdata.obs.columns:
                    cell_types_atac = new_mdata.obs[candidate].values[rna_only_mask]
                    break

            atac_metrics = compute_imputation_metrics(
                atac_imp, atac_gt, "atac", rna_only_bcs.values, cell_types_atac
            )
            atac_metrics["fraction"] = frac
            atac_metrics["seed"] = seed

            # Step 6: Smoothed ground truth metrics
            print("  Computing smoothed ground truth (50-NN)...")
            smoothed_rna, smoothed_atac = compute_smoothed_ground_truth(
                new_mdata, latent_key="X_MultiVI", k=50
            )

            rna_smoothed_metrics = compute_smoothed_metrics(
                rna_imp, smoothed_rna[atac_only_mask], "rna", atac_only_bcs.values
            )
            atac_smoothed_metrics = compute_smoothed_metrics(
                atac_imp, smoothed_atac[rna_only_mask], "atac", rna_only_bcs.values
            )

            # Merge smoothed metrics
            rna_metrics = rna_metrics.merge(
                rna_smoothed_metrics[["barcode", "spearman_r_smoothed"]],
                on="barcode", how="left",
            )
            atac_metrics = atac_metrics.merge(
                atac_smoothed_metrics[["barcode", "spearman_r_smoothed"]],
                on="barcode", how="left",
            )

            all_results.append(rna_metrics)
            all_results.append(atac_metrics)

            # Step 7: Integration quality metrics
            print("  Computing integration quality metrics...")
            modality_labels = np.array(new_mdata.obs["modality"])
            lisi_scores = compute_lisi(Z, modality_labels)

            # Cell type silhouette (if available)
            sil = np.nan
            if cell_types is not None:
                ct_all = None
                for candidate in [ct_key, f"rna:{ct_key}", "celltypist_prediction", "rna:celltypist_prediction"]:
                    if candidate in new_mdata.obs.columns:
                        ct_all = new_mdata.obs[candidate].values
                        break
                if ct_all is not None:
                    valid = ~pd.isna(ct_all)
                    if valid.sum() > 100:
                        sil = silhouette_score(Z[valid], ct_all[valid])

            integration_metrics.append({
                "fraction": frac,
                "seed": seed,
                "lisi_mean": float(np.mean(lisi_scores)),
                "lisi_median": float(np.median(lisi_scores)),
                "silhouette": float(sil),
            })

            print(f"  RNA Spearman r (median): {rna_metrics['spearman_r'].median():.3f}")
            print(f"  ATAC AUPRC (median): {atac_metrics['auprc'].median():.3f}")
            if "spearman_r_smoothed" in rna_metrics.columns:
                print(f"  RNA Spearman r smoothed (median): {rna_metrics['spearman_r_smoothed'].median():.3f}")

    # Combine all results
    all_results_df = pd.concat(all_results, ignore_index=True)
    all_results_df.to_csv(
        os.path.join(args.output_dir, "masking_sweep_results.csv"), index=False
    )

    # Summary statistics
    summary = {"fractions": fractions, "seeds": seeds, "n_cells": int(mdata_orig.n_obs)}
    for frac in fractions:
        frac_df = all_results_df[all_results_df["fraction"] == frac]
        rna_df = frac_df[frac_df["modality"] == "rna"]
        atac_df = frac_df[frac_df["modality"] == "atac"]
        summary[f"frac_{frac}"] = {
            "rna_spearman_r_median": float(rna_df["spearman_r"].median()) if len(rna_df) > 0 else None,
            "rna_spearman_r_smoothed_median": float(rna_df["spearman_r_smoothed"].median()) if "spearman_r_smoothed" in rna_df.columns and len(rna_df) > 0 else None,
            "atac_auprc_median": float(atac_df["auprc"].median()) if len(atac_df) > 0 else None,
            "atac_spearman_r_smoothed_median": float(atac_df["spearman_r_smoothed"].median()) if "spearman_r_smoothed" in atac_df.columns and len(atac_df) > 0 else None,
        }
    summary["integration"] = integration_metrics

    with open(os.path.join(args.output_dir, "masking_sweep_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    # Plots
    print("\nGenerating plots...")
    plot_degradation_curves(all_results_df, args.output_dir)
    if "cell_type" in all_results_df.columns:
        plot_celltype_breakdown(all_results_df, args.output_dir)
    plot_integration_quality(integration_metrics, args.output_dir)

    print(f"\nMasking sweep complete. Results in {args.output_dir}/")


if __name__ == "__main__":
    main()
