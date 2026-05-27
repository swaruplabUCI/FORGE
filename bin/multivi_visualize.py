#!/usr/bin/env python

import argparse
import os
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scvi.model import MULTIVI
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from sklearn.feature_extraction.text import TfidfTransformer


def parse_args():
    p = argparse.ArgumentParser(description="Visualize MultiVI integration results")
    p.add_argument("--mudata", required=True, help="MultiVI-integrated MuData file")
    p.add_argument("--model_dir", required=True, help="Saved MultiVI model directory")
    p.add_argument("--output_dir", default="multivi_plots", help="Output directory for plots")
    p.add_argument("--cell_type_key", default="cell_type", help="Cell type annotation key")
    p.add_argument("--batch_key", default="sample_id", help="Batch key")
    p.add_argument("--run_imputation", type=str, default="true", help="Run imputation analysis")
    p.add_argument("--run_differential", type=str, default="false", help="Run differential analysis")
    return p.parse_args()


def compute_umap_and_cluster(mdata):
    """Compute UMAP and Leiden clustering on MultiVI latent space"""
    print("Computing UMAP on MultiVI latent representation...")

    # Use MultiVI latent representation for neighbors
    if "X_MultiVI" not in mdata.obsm_keys():
        raise ValueError("X_MultiVI not found in mdata.obsm. Cannot compute neighbors/UMAP.")

    sc.pp.neighbors(mdata, use_rep="X_MultiVI", n_neighbors=30)
    sc.tl.umap(mdata, min_dist=0.3)
    sc.tl.leiden(mdata, resolution=0.5, key_added="leiden_multivi")

    return mdata


def plot_modality_integration(mdata, output_dir, cell_type_key="cell_type"):
    """Create UMAP plots colored by modality, batch, clusters, and cell type"""
    print("Creating modality integration plots...")

    # Plot 1: Color by modality (if exists)
    if "modality" in mdata.obs.columns:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        sc.pl.umap(
            mdata,
            color="modality",
            ax=ax,
            show=False,
            title="MultiVI: Modality Integration",
            legend_loc="right margin",
        )
        plt.savefig(os.path.join(output_dir, "umap_modality.pdf"), bbox_inches="tight")
        plt.close()

    # Plot 2: Color by batch
    if "batch" in mdata.obs.columns or "sample_id" in mdata.obs.columns:
        batch_key = "batch" if "batch" in mdata.obs.columns else "sample_id"
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        sc.pl.umap(
            mdata,
            color=batch_key,
            ax=ax,
            show=False,
            title="MultiVI: Batch Integration",
            legend_loc="right margin",
        )
        plt.savefig(os.path.join(output_dir, "umap_batch.pdf"), bbox_inches="tight")
        plt.close()

    # Plot 3: Color by Leiden clusters
    if "leiden_multivi" in mdata.obs.columns:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        sc.pl.umap(
            mdata,
            color="leiden_multivi",
            ax=ax,
            show=False,
            title="MultiVI: Leiden Clustering",
            legend_loc="right margin",
        )
        plt.savefig(os.path.join(output_dir, "umap_leiden.pdf"), bbox_inches="tight")
        plt.close()

    # Plot 4: Color by cell type - ENHANCED DETECTION
    # Try multiple possible cell type column names
    cell_type_found = False
    for possible_key in [
        cell_type_key,
        "celltypist_prediction",
        "scanvi_prediction",
        "cell_type",
        "celltype",
        "rna:celltypist_prediction",
        "rna:scanvi_prediction",
        "atac:scanvi_prediction",
    ]:
        if possible_key in mdata.obs.columns:
            _vals = mdata.obs[possible_key].dropna().unique()
            if len(_vals) == 1 and _vals[0] == 'Unknown':
                continue  # skip all-Unknown columns
            print(f"Found cell type annotations in column: {possible_key}")

            # Inject explicit color palette — scanpy's built-in godsnot_102 palette
            # tops out at 102 colors; CellTypist brain labels routinely produce 100-200
            # categories, causing a grayscale fallback.  Mirror the ATAC injection logic.
            import matplotlib.colors as _mc_rna
            _rna_cats = sorted(mdata.obs[possible_key].dropna().astype(str).unique())
            _n_rna = len(_rna_cats)
            _cmap_rna = plt.get_cmap(
                "gist_ncar" if _n_rna > 102 else "turbo" if _n_rna > 20 else "tab20"
            )
            _rna_hex = [_mc_rna.to_hex(_cmap_rna(i / max(_n_rna - 1, 1))) for i in range(_n_rna)]
            if ':' in possible_key:
                _mod_name_rna, _obs_k_rna = possible_key.split(':', 1)
                if hasattr(mdata, 'mod') and _mod_name_rna in mdata.mod:
                    mdata.mod[_mod_name_rna].uns[f'{_obs_k_rna}_colors'] = _rna_hex
            else:
                mdata.uns[f'{possible_key}_colors'] = _rna_hex
            print(f"  Injected {_n_rna}-color palette for '{possible_key}'")

            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            # For very high category counts suppress the legend (unreadable anyway);
            # the palette itself is still colour-correct in the scatter.
            _legend_loc_rna = "right margin" if _n_rna <= 60 else "none"
            sc.pl.umap(
                mdata,
                color=possible_key,
                ax=ax,
                show=False,
                title=f"MultiVI: RNA Cell Types ({_n_rna} types)",
                legend_loc=_legend_loc_rna,
                legend_fontsize=6,
                legend_fontoutline=1,
                frameon=False,
            )
            plt.savefig(os.path.join(output_dir, "umap_celltype.pdf"), bbox_inches="tight")
            plt.close()

            # Also save cell type key for later use
            mdata.uns["cell_type_key_used"] = possible_key
            cell_type_found = True
            break

    if not cell_type_found:
        print("Warning: No cell type annotations found for UMAP visualization")

    # Plot 4b: ATAC cell type (FIX-67)
    atac_ct_key = None
    for candidate in ['atac:cell_type', 'atac:scanvi_prediction']:
        if candidate in mdata.obs.columns:
            atac_ct_key = candidate
            break
    if atac_ct_key is None and hasattr(mdata, 'mod') and 'atac' in mdata.mod:
        if 'cell_type' in mdata.mod['atac'].obs.columns:
            mdata.obs['atac_cell_type'] = mdata.mod['atac'].obs['cell_type']
            atac_ct_key = 'atac_cell_type'
    if atac_ct_key:
        print(f"Plotting ATAC cell type UMAP using: {atac_ct_key}")
        # Inject full palette so high-N ATAC types (e.g. 60+ brain subclasses) get distinct colors
        import matplotlib.colors as _mc
        _atac_cats = sorted(mdata.obs[atac_ct_key].dropna().astype(str).unique())
        _n_atac = len(_atac_cats)
        _cmap_atac = plt.get_cmap("gist_ncar" if _n_atac > 102 else "turbo" if _n_atac > 20 else "tab20")
        _atac_hex = [_mc.to_hex(_cmap_atac(i / max(_n_atac - 1, 1))) for i in range(_n_atac)]
        if ':' in atac_ct_key:
            _mod_name, _obs_k = atac_ct_key.split(':', 1)
            if _mod_name in mdata.mod:
                mdata.mod[_mod_name].uns[f'{_obs_k}_colors'] = _atac_hex
        else:
            mdata.uns[f'{atac_ct_key}_colors'] = _atac_hex
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        sc.pl.umap(
            mdata, color=atac_ct_key, ax=ax, show=False,
            title="MultiVI: ATAC Cell Types",
            legend_loc="right margin", legend_fontsize=6,
            legend_fontoutline=1, frameon=False,
        )
        plt.savefig(os.path.join(output_dir, "umap_atac_celltype.pdf"), bbox_inches="tight")
        plt.close()
        print(f"Saved: umap_atac_celltype.pdf")


def evaluate_imputation(mdata, mvi_model, output_dir):
    """Evaluate cross-modal imputation quality"""
    print("Evaluating cross-modal imputation...")

    metrics = {}

    # Get normalized expression and accessibility from model
    expr_norm = mvi_model.get_normalized_expression(return_numpy=True)
    acc_norm = mvi_model.get_normalized_accessibility(return_numpy=True)

    # Identify RNA-only and ATAC-only cells if modality info exists
    if "modality" in mdata.obs.columns:
        rna_only_mask = mdata.obs["modality"] == "expression"
        atac_only_mask = mdata.obs["modality"] == "accessibility"
        paired_mask = mdata.obs["modality"] == "paired"

        # For paired cells, compare observed vs imputed
        if paired_mask.sum() > 0:
            print(f"Found {paired_mask.sum()} paired cells for validation...")

            # Sample genes and peaks for correlation analysis
            n_features = min(100, mdata.mod["rna"].n_vars)
            gene_indices = np.random.choice(mdata.mod["rna"].n_vars, n_features, replace=False)

            correlations = []
            for idx in gene_indices:
                obs = mdata.mod["rna"].X[paired_mask, idx].toarray().flatten()
                imp = expr_norm[paired_mask, idx]
                if obs.std() > 0 and imp.std() > 0:
                    corr = pearsonr(obs, imp)[0]
                    correlations.append(corr)

            metrics["mean_gene_correlation"] = np.mean(correlations)
            metrics["median_gene_correlation"] = np.median(correlations)

            # Create correlation scatter plot for top variable genes
            fig, axes = plt.subplots(2, 2, figsize=(10, 10))

            # Select top variable genes
            sc.pp.highly_variable_genes(mdata.mod["rna"], n_top_genes=4, subset=False)
            top_genes = mdata.mod["rna"].var_names[mdata.mod["rna"].var["highly_variable"]][:4]

            for i, gene in enumerate(top_genes):
                ax = axes[i // 2, i % 2]
                gene_idx = mdata.mod["rna"].var_names.get_loc(gene)
                obs = mdata.mod["rna"].X[paired_mask, gene_idx].toarray().flatten()
                imp = expr_norm[paired_mask, gene_idx]

                ax.scatter(obs, imp, alpha=0.5, s=1)
                ax.set_xlabel(f"Observed {gene}")
                ax.set_ylabel(f"Imputed {gene}")

                # Add correlation to plot
                if obs.std() > 0 and imp.std() > 0:
                    corr = pearsonr(obs, imp)[0]
                    ax.text(0.05, 0.95, f"r = {corr:.3f}", transform=ax.transAxes)

            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "imputation_correlation_genes.pdf"))
            plt.close()

    # Save metrics
    pd.DataFrame([metrics]).to_csv("multivi_imputation_metrics.csv", index=False)

    return metrics


def run_differential_analysis(mdata, mvi_model, cell_type_key, output_dir):
    """Run differential expression/accessibility analysis"""
    print("Running differential analysis between cell types...")

    if cell_type_key not in mdata.obs.columns:
        print(f"Cell type key '{cell_type_key}' not found. Skipping differential analysis.")
        return None, None

    # Get unique cell types
    cell_types = mdata.obs[cell_type_key].unique()

    if len(cell_types) < 2:
        print("Need at least 2 cell types for differential analysis.")
        return None, None

    # Example: Compare first two cell types
    ct1, ct2 = cell_types[0], cell_types[1]
    print(f"Comparing {ct1} vs {ct2}...")

    idx1 = mdata.obs[cell_type_key] == ct1
    idx2 = mdata.obs[cell_type_key] == ct2

    # Differential expression
    de_results = mvi_model.differential_expression(
        idx1=idx1,
        idx2=idx2,
    )
    de_results.to_csv("multivi_de_results.csv")

    # Differential accessibility
    da_results = mvi_model.differential_accessibility(
        idx1=idx1,
        idx2=idx2,
    )
    da_results.to_csv("multivi_da_results.csv")

    # Create volcano plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # DE volcano
    ax = axes[0]
    de_results["log10_bayes"] = np.abs(de_results["bayes_factor"])
    significant = de_results["is_de_fdr_0.05"]

    ax.scatter(
        de_results.loc[~significant, "lfc_mean"],
        de_results.loc[~significant, "log10_bayes"],
        alpha=0.5,
        s=1,
        c="gray",
    )
    ax.scatter(
        de_results.loc[significant, "lfc_mean"],
        de_results.loc[significant, "log10_bayes"],
        alpha=0.5,
        s=1,
        c="red",
    )
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel("|Bayes Factor|")
    ax.set_title(f"DE: {ct1} vs {ct2}")
    ax.axhline(y=3, color="k", linestyle="--", alpha=0.3)
    ax.axvline(x=0, color="k", linestyle="--", alpha=0.3)

    # DA volcano
    ax = axes[1]
    da_results["log10_bayes"] = np.abs(da_results["bayes_factor"])
    significant = da_results["is_da_fdr"]

    ax.scatter(
        da_results.loc[~significant, "effect_size"],
        da_results.loc[~significant, "log10_bayes"],
        alpha=0.5,
        s=1,
        c="gray",
    )
    ax.scatter(
        da_results.loc[significant, "effect_size"],
        da_results.loc[significant, "log10_bayes"],
        alpha=0.5,
        s=1,
        c="blue",
    )
    ax.set_xlabel("Effect Size")
    ax.set_ylabel("|Bayes Factor|")
    ax.set_title(f"DA: {ct1} vs {ct2}")
    ax.axhline(y=3, color="k", linestyle="--", alpha=0.3)
    ax.axvline(x=0, color="k", linestyle="--", alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "differential_volcano.pdf"))
    plt.close()

    return de_results, da_results


def create_summary_report(mdata, metrics, output_dir):
    """Create HTML summary report"""
    print("Creating summary report...")

    html_content = f"""
    <html>
    <head>
        <title>MultiVI Integration Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1 {{ color: #333; }}
            h2 {{ color: #666; }}
            .metric {{ background-color: #f0f0f0; padding: 10px; margin: 10px 0; }}
            .plot {{ margin: 20px 0; }}
            img {{ max-width: 100%; height: auto; }}
        </style>
    </head>
    <body>
        <h1>MultiVI Integration Report</h1>

        <h2>Dataset Overview</h2>
        <div class="metric">
            <p><b>Total cells:</b> {mdata.n_obs}</p>
            <p><b>RNA features:</b> {mdata.mod['rna'].n_vars if 'rna' in mdata.mod else 'N/A'}</p>
            <p><b>ATAC features:</b> {mdata.mod['atac'].n_vars if 'atac' in mdata.mod else 'N/A'}</p>
            <p><b>Number of batches:</b> {mdata.obs['sample_id'].nunique() if 'sample_id' in mdata.obs else 'N/A'}</p>
            <p><b>Number of clusters:</b> {mdata.obs['leiden_multivi'].nunique() if 'leiden_multivi' in mdata.obs else 'N/A'}</p>
        </div>

        <h2>Integration Quality</h2>
        <div class="plot">
            <h3>UMAP - Modality Integration</h3>
            <img src="multivi_plots/umap_modality.pdf" alt="Modality UMAP">
        </div>

        <div class="plot">
            <h3>UMAP - Batch Integration</h3>
            <img src="multivi_plots/umap_batch.pdf" alt="Batch UMAP">
        </div>

        <div class="plot">
            <h3>UMAP - Cell Types</h3>
            <img src="multivi_plots/umap_celltype.pdf" alt="Cell Type UMAP">
        </div>
    """

    if metrics:
        html_content += f"""
        <h2>Imputation Metrics</h2>
        <div class="metric">
            <p><b>Mean gene correlation:</b> {metrics.get('mean_gene_correlation', 'N/A'):.3f}</p>
            <p><b>Median gene correlation:</b> {metrics.get('median_gene_correlation', 'N/A'):.3f}</p>
        </div>

        <div class="plot">
            <h3>Imputation Correlation</h3>
            <img src="multivi_plots/imputation_correlation_genes.pdf" alt="Imputation Correlation">
        </div>
        """

    html_content += """
    </body>
    </html>
    """

    with open("multivi_visualization_report.html", "w") as f:
        f.write(html_content)


def main():
    args = parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load integrated MuData and model
    print("Loading MultiVI-integrated MuData...")
    mdata = mu.read_h5mu(args.mudata)

    print("Loading MultiVI model...")
    mvi = MULTIVI.load(args.model_dir, adata=mdata)

    # Compute UMAP and clustering
    mdata = compute_umap_and_cluster(mdata)

    # Detect cell type key
    cell_type_key = args.cell_type_key
    for possible_key in [
        "cell_type",
        "scanvi_prediction",
        "rna:scanvi_prediction",
        "atac:scanvi_prediction",
    ]:
        if possible_key in mdata.obs.columns:
            print(f"Detected cell type key: {possible_key}")
            cell_type_key = possible_key
            break

    # Create integration plots with detected cell type key
    plot_modality_integration(mdata, args.output_dir, cell_type_key)

    # Evaluate imputation quality
    metrics = {}
    if args.run_imputation.lower() == "true":
        metrics = evaluate_imputation(mdata, mvi, args.output_dir)
    else:
        # Write empty metrics CSV so Nextflow output declaration is satisfied
        print("Imputation skipped — writing empty metrics CSV.")
        pd.DataFrame(columns=["mean_gene_correlation", "median_gene_correlation"]).to_csv(
            "multivi_imputation_metrics.csv", index=False
        )

    # Run differential analysis with detected cell type key
    if args.run_differential.lower() == "true":
        de_results, da_results = run_differential_analysis(
            mdata, mvi, cell_type_key, args.output_dir  # Use detected key
        )

    # Create summary report
    create_summary_report(mdata, metrics, args.output_dir)

    print("MultiVI visualization complete!")


if __name__ == "__main__":
    main()
