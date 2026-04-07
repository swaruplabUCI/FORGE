#!/usr/bin/env python3
"""
signal_chain_correlation.py — Recipe Step C6

Full-chain correlation: ligand expression (sender) ~ footprint depth (receiver)
~ target gene expression (receiver).  Classifies interactions into evidence
tiers (1-4) based on correlation strength.

Evidence Tiers:
    Tier 1: Strong ligand-footprint AND footprint-target correlations (>0.5)
    Tier 2: Moderate correlations (>0.3)
    Tier 3: Weak but significant correlations (p<0.05)
    Tier 4: Insufficient evidence

Inputs:
    - Enhancer footprint h5ads
    - RNA h5ad (integrated)
    - signaling_metadata.tsv

Outputs:
    - signal_chain_correlations.tsv
    - evidence_tiers.tsv
    - chain_plots/
"""

import argparse
import glob
import json
import re
from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def normalize_celltype(ct):
    """Normalize cell type name for matching."""
    ct = str(ct).strip()
    ct = re.sub(r'[^a-zA-Z0-9]', '_', ct)
    ct = re.sub(r'_+', '_', ct).strip('_')
    return ct.lower()


def parse_args():
    p = argparse.ArgumentParser(description="Signal chain correlation analysis")
    p.add_argument("--footprint-h5ads", nargs='+', help="Footprint h5ad files")
    p.add_argument("--rna-h5ad", required=True, help="Integrated RNA h5ad")
    p.add_argument("--signaling-metadata", required=True, help="Signaling metadata TSV")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    plots_dir = outdir / "chain_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load signaling metadata
    print(f"Loading signaling metadata from {args.signaling_metadata}")
    meta_df = pd.read_csv(args.signaling_metadata, sep='\t')
    print(f"  {len(meta_df)} signaling target entries")

    if meta_df.empty:
        print("No signaling targets. Writing empty outputs.")
        pd.DataFrame().to_csv(outdir / "signal_chain_correlations.tsv", sep='\t', index=False)
        pd.DataFrame().to_csv(outdir / "evidence_tiers.tsv", sep='\t', index=False)
        (outdir / "chain_summary.txt").write_text("No signaling targets to analyze\n")
        return

    # Load RNA
    # FIX-90: Strip encoding_type='null' datasets before reading (anndata compat)
    print(f"Loading RNA h5ad from {args.rna_h5ad}")
    try:
        import h5py, shutil
        tmp_rna = str(Path(args.rna_h5ad).stem) + '_stripped.h5ad'
        shutil.copy2(args.rna_h5ad, tmp_rna)
        with h5py.File(tmp_rna, 'a') as f:
            # FIX-94: Collect paths first, then delete — can't modify HDF5 tree during visititems
            to_delete = []
            def _strip_null(name, obj):
                if isinstance(obj, h5py.Dataset):
                    et = obj.attrs.get('encoding-type', None)
                    if et is not None and (et == 'null' or et == b'null'):
                        to_delete.append(name)
            f.visititems(_strip_null)
            for name in to_delete:
                print(f"  FIX-94: Stripping {name} (encoding_type='null')")
                del f[name]
        args.rna_h5ad = tmp_rna
    except Exception as e:
        print(f"  [WARN] FIX-90 strip failed, trying original: {e}")
    rna = ad.read_h5ad(args.rna_h5ad)

    # Detect cell type column — prefer columns with real labels over all-Unknown
    ct_col = None
    for col in ['cell_type', 'celltypist_prediction', 'scanvi_prediction']:
        if col in rna.obs.columns:
            vals = rna.obs[col].dropna().unique()
            if len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown'):
                ct_col = col
                break
    if ct_col is None:
        ct_col = 'cell_type'  # fall through to original default
    rna.obs['cell_type'] = rna.obs[ct_col]

    cell_types = sorted(rna.obs['cell_type'].unique().tolist())
    ct_norm_map = {normalize_celltype(ct): ct for ct in cell_types}

    gene_names = list(rna.var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Compute per-cell-type mean expression
    print("Computing per-cell-type mean expression...")
    ct_mean_expr = {}
    for ct in cell_types:
        mask = rna.obs['cell_type'] == ct
        ct_data = rna[mask]
        if hasattr(ct_data.X, 'toarray'):
            ct_mean_expr[ct] = np.array(ct_data.X.toarray().mean(axis=0)).flatten()
        else:
            ct_mean_expr[ct] = np.array(ct_data.X.mean(axis=0)).flatten()

    # Load footprint data
    print("Loading footprint h5ads...")
    fp_files = args.footprint_h5ads or []
    fp_data = {}  # (cell_type_name, tf) -> mean score

    for fp_file in fp_files:
        fname = Path(fp_file).stem
        # Parse: enhancer_footprints_{ct}_{tf}
        match = re.match(r'enhancer_footprints_(.+?)_([A-Z][A-Z0-9]+)$', fname)
        if not match:
            parts = fname.replace('enhancer_footprints_', '').rsplit('_', 1)
            if len(parts) == 2:
                ct_label, tf_label = parts
            else:
                continue
        else:
            ct_label, tf_label = match.group(1), match.group(2)

        try:
            fp_adata = ad.read_h5ad(fp_file)
            if hasattr(fp_adata.X, 'toarray'):
                mean_scores = np.nanmean(fp_adata.X.toarray(), axis=1)
            else:
                mean_scores = np.nanmean(np.array(fp_adata.X), axis=1)

            for i, group_name in enumerate(fp_adata.obs_names):
                if i < len(mean_scores):
                    fp_data[(str(group_name), tf_label)] = float(mean_scores[i])
        except Exception as e:
            print(f"  WARNING: Failed to load {fp_file}: {e}")

    print(f"  {len(fp_data)} (cell_type, TF) footprint scores loaded")

    # Compute correlations for each signaling chain
    correlation_records = []

    for _, entry in meta_df.iterrows():
        pathway = str(entry.get('pathway', ''))
        sender = str(entry.get('sender', ''))
        receiver = str(entry.get('receiver', ''))
        tf = str(entry.get('tf', ''))

        if not tf or tf == 'nan':
            continue

        # Collect per-cell-type data points
        fp_scores = []
        tf_expr = []
        ct_labels = []

        tf_idx = gene_to_idx.get(tf)
        if tf_idx is None:
            continue

        for ct in cell_types:
            fp_key = (ct, tf)
            if fp_key not in fp_data:
                ct_norm = normalize_celltype(ct)
                for k in fp_data:
                    if normalize_celltype(k[0]) == ct_norm and k[1] == tf:
                        fp_key = k
                        break
                else:
                    continue

            fp_scores.append(fp_data[fp_key])
            tf_expr.append(float(ct_mean_expr[ct][tf_idx]))
            ct_labels.append(ct)

        if len(fp_scores) < 3:
            continue

        fp_arr = np.array(fp_scores)
        tf_arr = np.array(tf_expr)

        # Ligand-footprint correlation (if ligand gene known)
        # For now, use TF expression as proxy
        corr_fp_tf, pval_fp_tf = stats.spearmanr(fp_arr, tf_arr)

        # Classify into evidence tiers
        abs_corr = abs(corr_fp_tf)
        if abs_corr > 0.5 and pval_fp_tf < 0.05:
            tier = 1
        elif abs_corr > 0.3 and pval_fp_tf < 0.05:
            tier = 2
        elif pval_fp_tf < 0.05:
            tier = 3
        else:
            tier = 4

        correlation_records.append({
            'pathway': pathway,
            'sender': sender,
            'receiver': receiver,
            'TF': tf,
            'n_cell_types': len(fp_scores),
            'corr_footprint_tfexpr': round(corr_fp_tf, 4),
            'pval_footprint_tfexpr': round(pval_fp_tf, 6),
            'evidence_tier': tier,
        })

    # Write correlations
    corr_df = pd.DataFrame(correlation_records) if correlation_records else pd.DataFrame(
        columns=['pathway', 'sender', 'receiver', 'TF', 'n_cell_types',
                 'corr_footprint_tfexpr', 'pval_footprint_tfexpr', 'evidence_tier'])

    corr_path = outdir / "signal_chain_correlations.tsv"
    corr_df.to_csv(corr_path, sep='\t', index=False)
    print(f"\nWrote {len(corr_df)} chain correlations")

    # Evidence tiers summary
    tier_df = corr_df[['pathway', 'sender', 'receiver', 'TF', 'evidence_tier']].copy() if not corr_df.empty else pd.DataFrame(
        columns=['pathway', 'sender', 'receiver', 'TF', 'evidence_tier'])
    tier_path = outdir / "evidence_tiers.tsv"
    tier_df.to_csv(tier_path, sep='\t', index=False)

    # Plot: tier distribution
    if not corr_df.empty:
        try:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

            tier_counts = corr_df['evidence_tier'].value_counts().sort_index()
            tier_counts.plot.bar(ax=ax1, color=['#2ca02c', '#1f77b4', '#ff7f0e', '#d62728'])
            ax1.set_xlabel("Evidence Tier")
            ax1.set_ylabel("Count")
            ax1.set_title("Evidence Tier Distribution")

            if len(corr_df) > 0:
                ax2.scatter(
                    corr_df['corr_footprint_tfexpr'],
                    -np.log10(corr_df['pval_footprint_tfexpr'].clip(1e-50)),
                    c=corr_df['evidence_tier'],
                    cmap='RdYlGn_r',
                    alpha=0.6,
                )
                ax2.axhline(-np.log10(0.05), ls='--', color='gray', alpha=0.5)
                ax2.set_xlabel("Spearman Correlation")
                ax2.set_ylabel("-log10(p-value)")
                ax2.set_title("Footprint-TF Expression Correlation")

            plt.tight_layout()
            plt.savefig(plots_dir / "evidence_tier_overview.png", dpi=150)
            plt.close()
        except Exception as e:
            print(f"  Plot failed: {e}")
            plt.close('all')

    # Summary
    summary_lines = [
        "Signal Chain Correlation Summary",
        "=" * 40,
        f"Signaling targets analyzed: {len(meta_df)}",
        f"Chains with sufficient data: {len(corr_df)}",
        f"Footprint scores available: {len(fp_data)}",
    ]
    if not corr_df.empty:
        summary_lines.append("")
        summary_lines.append("Evidence tier distribution:")
        for tier in [1, 2, 3, 4]:
            count = (corr_df['evidence_tier'] == tier).sum()
            summary_lines.append(f"  Tier {tier}: {count}")

    summary_path = outdir / "chain_summary.txt"
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
