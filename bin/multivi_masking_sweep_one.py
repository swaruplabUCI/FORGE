#!/usr/bin/env python
"""
MultiVI Masking Sweep — single (fraction, seed) fit.

Runs one iteration of the artificial-unpairing benchmark, then writes three
small files to --output_dir:
  fit_frac{F}_seed{S}_rna.csv
  fit_frac{F}_seed{S}_atac.csv
  fit_frac{F}_seed{S}_integration.json

A separate aggregator (`multivi_masking_sweep_aggregate.py`) merges these into
the final `masking_sweep_results.csv` + `masking_sweep_summary.json` + plots.

Process isolation is the point: one fit per python invocation so the scvi-tools
MultiVI host-RAM leak cannot accumulate across fits.
"""

import argparse
import json
import os
import sys
import warnings

import muon as mu
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multivi_masking_sweep import (
    artificially_unpair,
    compute_imputation_metrics,
    compute_lisi,
    compute_smoothed_ground_truth,
    compute_smoothed_metrics,
)
from scvi.model import MULTIVI
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def parse_args():
    p = argparse.ArgumentParser(description="MultiVI masking sweep — single fit")
    p.add_argument("--mudata", required=True, help="Pre-MultiVI integrated MuData (.h5mu)")
    p.add_argument("--output_dir", default=".", help="Output directory")
    p.add_argument("--fraction", type=float, required=True, help="Unpairing fraction for this fit")
    p.add_argument("--seed", type=int, required=True, help="Random seed for this fit")
    p.add_argument("--n_epochs", type=int, default=200, help="Training epochs")
    p.add_argument("--batch_key", default="sample_id")
    p.add_argument("--cell_type_key", default="celltypist_prediction")
    p.add_argument("--n_latent", type=int, default=20)
    p.add_argument("--modality_weights", default="equal")
    p.add_argument("--modality_penalty", default="Jeffreys")
    return p.parse_args()


def fit_tag(frac, seed):
    return f"frac{frac}_seed{seed}"


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    frac = args.fraction
    seed = args.seed
    tag = fit_tag(frac, seed)

    print(f"Loading MuData from {args.mudata}...")
    mdata_orig = mu.read_h5mu(args.mudata)
    print(f"  Cells: {mdata_orig.n_obs}")
    print(f"  RNA features: {mdata_orig.mod['rna'].n_vars}")
    print(f"  ATAC features: {mdata_orig.mod['atac'].n_vars}")

    print(f"\n{'='*60}\n  Fraction: {frac}, Seed: {seed}\n{'='*60}")

    print("  Artificially unpairing cells...")
    new_mdata, held_out_rna, held_out_atac, _ = artificially_unpair(mdata_orig, frac, seed)
    n_paired = (np.array(new_mdata.obs["modality"]) == "paired").sum()
    n_rna_only = (np.array(new_mdata.obs["modality"]) == "expression").sum()
    n_atac_only = (np.array(new_mdata.obs["modality"]) == "accessibility").sum()
    print(f"  Paired: {n_paired}, RNA-only: {n_rna_only}, ATAC-only: {n_atac_only}")

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

    print("  Computing latent representation and imputations...")
    Z = mvi.get_latent_representation(modality="joint")
    new_mdata.obsm["X_MultiVI"] = Z

    expr_imputed = mvi.get_normalized_expression(return_numpy=True)
    acc_imputed = mvi.get_normalized_accessibility(return_numpy=True)

    print("  Evaluating RNA imputation on ATAC-only cells...")
    atac_only_mask = np.array(new_mdata.obs["modality"]) == "accessibility"
    atac_only_bcs = new_mdata.obs_names[atac_only_mask]
    orig_bcs = [bc.replace("_atac_only", "") for bc in atac_only_bcs]

    rna_gt = np.array([held_out_rna[bc] for bc in orig_bcs])
    rna_imp = expr_imputed[atac_only_mask]

    ct_key = args.cell_type_key
    cell_types = None
    for candidate in [ct_key, f"rna:{ct_key}", "celltypist_prediction", "rna:celltypist_prediction"]:
        if candidate in new_mdata.obs.columns:
            cell_types = new_mdata.obs[candidate].values[atac_only_mask]
            break

    rna_metrics = compute_imputation_metrics(rna_imp, rna_gt, "rna", atac_only_bcs.values, cell_types)
    rna_metrics["fraction"] = frac
    rna_metrics["seed"] = seed

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

    atac_metrics = compute_imputation_metrics(atac_imp, atac_gt, "atac", rna_only_bcs.values, cell_types_atac)
    atac_metrics["fraction"] = frac
    atac_metrics["seed"] = seed

    print("  Computing smoothed ground truth (50-NN)...")
    smoothed_rna, smoothed_atac = compute_smoothed_ground_truth(new_mdata, latent_key="X_MultiVI", k=50)

    rna_smoothed_metrics = compute_smoothed_metrics(rna_imp, smoothed_rna[atac_only_mask], "rna", atac_only_bcs.values)
    atac_smoothed_metrics = compute_smoothed_metrics(atac_imp, smoothed_atac[rna_only_mask], "atac", rna_only_bcs.values)

    rna_metrics = rna_metrics.merge(
        rna_smoothed_metrics[["barcode", "spearman_r_smoothed"]], on="barcode", how="left"
    )
    atac_metrics = atac_metrics.merge(
        atac_smoothed_metrics[["barcode", "spearman_r_smoothed"]], on="barcode", how="left"
    )

    print("  Computing integration quality metrics...")
    modality_labels = np.array(new_mdata.obs["modality"])
    lisi_scores = compute_lisi(Z, modality_labels)

    sil = np.nan
    ct_all = None
    for candidate in [ct_key, f"rna:{ct_key}", "celltypist_prediction", "rna:celltypist_prediction"]:
        if candidate in new_mdata.obs.columns:
            ct_all = new_mdata.obs[candidate].values
            break
    if ct_all is not None:
        valid = ~pd.isna(ct_all)
        if valid.sum() > 100:
            sil = silhouette_score(Z[valid], ct_all[valid])

    integration = {
        "fraction": float(frac),
        "seed": int(seed),
        "lisi_mean": float(np.mean(lisi_scores)),
        "lisi_median": float(np.median(lisi_scores)),
        "silhouette": float(sil),
    }

    print(f"  RNA Spearman r (median): {rna_metrics['spearman_r'].median():.3f}")
    print(f"  ATAC AUPRC (median): {atac_metrics['auprc'].median():.3f}")
    if "spearman_r_smoothed" in rna_metrics.columns:
        print(f"  RNA Spearman r smoothed (median): {rna_metrics['spearman_r_smoothed'].median():.3f}")

    rna_path = os.path.join(args.output_dir, f"fit_{tag}_rna.csv")
    atac_path = os.path.join(args.output_dir, f"fit_{tag}_atac.csv")
    int_path = os.path.join(args.output_dir, f"fit_{tag}_integration.json")

    rna_metrics.to_csv(rna_path, index=False)
    atac_metrics.to_csv(atac_path, index=False)
    with open(int_path, "w") as f:
        json.dump(integration, f, indent=2)

    print(f"\nFit {tag} complete. Outputs in {args.output_dir}/")


if __name__ == "__main__":
    main()
