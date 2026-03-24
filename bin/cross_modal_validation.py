#!/usr/bin/env python3
"""
cross_modal_validation.py — Recipe Step B6

Correlate footprint depth with TF expression and target gene expression
across cell states.  Per-cell-type mean footprint score x mean TF expression
x mean target gene expression -> Spearman correlations -> classify as
activator/repressor/indirect per recipe criteria.

Inputs:
    - Enhancer footprint h5ads (from ENHANCER_FOOTPRINTING)
    - RNA h5ad (integrated, with cell_type in obs)
    - Target genes JSON ({TF: [genes]})

Outputs:
    - cross_modal_validation.tsv
    - validation_plots/
    - validation_summary.txt
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


def parse_args():
    p = argparse.ArgumentParser(description="Cross-modal validation of enhancer footprints")
    p.add_argument("--footprint-dir", required=True, help="Directory containing footprint h5ads")
    p.add_argument("--footprint-h5ads", nargs='+', help="Footprint h5ad files")
    p.add_argument("--rna-h5ad", required=True, help="Integrated RNA h5ad")
    p.add_argument("--target-genes", required=True, help="Target genes JSON")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def normalize_celltype(ct):
    """Normalize cell type name for matching across modalities."""
    ct = str(ct).strip()
    ct = re.sub(r'[^a-zA-Z0-9]', '_', ct)
    ct = re.sub(r'_+', '_', ct).strip('_')
    return ct.lower()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    plots_dir = outdir / "validation_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Load target genes
    print(f"Loading target genes from {args.target_genes}")
    with open(args.target_genes) as f:
        target_genes = json.load(f)
    print(f"  {len(target_genes)} TFs with target gene lists")

    # Load RNA data
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
    if 'cell_type' not in rna.obs.columns:
        # Try alternative column names
        for col in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_annotation']:
            if col in rna.obs.columns:
                rna.obs['cell_type'] = rna.obs[col]
                break
    if 'cell_type' not in rna.obs.columns:
        raise ValueError("RNA h5ad missing cell_type column")

    cell_types = sorted(rna.obs['cell_type'].unique().tolist())
    print(f"  {len(cell_types)} cell types, {rna.n_obs} cells, {rna.n_vars} genes")

    # Build cell-type normalized lookup
    ct_norm_map = {normalize_celltype(ct): ct for ct in cell_types}

    # Compute per-cell-type mean expression for all genes
    print("Computing per-cell-type mean expression...")
    ct_mean_expr = {}
    for ct in cell_types:
        mask = rna.obs['cell_type'] == ct
        ct_data = rna[mask]
        if hasattr(ct_data.X, 'toarray'):
            ct_mean_expr[ct] = np.array(ct_data.X.toarray().mean(axis=0)).flatten()
        else:
            ct_mean_expr[ct] = np.array(ct_data.X.mean(axis=0)).flatten()

    gene_names = list(rna.var_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Load footprint h5ads and extract per-cell-type mean footprint scores
    print("Loading footprint h5ads...")
    fp_files = args.footprint_h5ads if args.footprint_h5ads else glob.glob(
        str(Path(args.footprint_dir) / "enhancer_footprints_*.h5ad"))

    fp_data = {}  # (ct, tf) -> mean_footprint_score
    for fp_file in fp_files:
        fname = Path(fp_file).stem  # e.g., enhancer_footprints_CD4_T_STAT1
        # Parse cell_type and tf from filename
        match = re.match(r'enhancer_footprints_(.+?)_([A-Z][A-Z0-9]+)$', fname)
        if not match:
            # Try alternative: split on last underscore that precedes a TF-like name
            parts = fname.replace('enhancer_footprints_', '').rsplit('_', 1)
            if len(parts) == 2:
                ct_label, tf_label = parts
            else:
                print(f"  Skipping unparseable file: {fname}")
                continue
        else:
            ct_label, tf_label = match.group(1), match.group(2)

        try:
            fp_adata = ad.read_h5ad(fp_file)
            # Average footprint score across regions and scales
            if hasattr(fp_adata.X, 'toarray'):
                mean_scores = np.nanmean(fp_adata.X.toarray(), axis=1)
            else:
                mean_scores = np.nanmean(np.array(fp_adata.X), axis=1)

            # Map to cell types via obs_names
            for i, group_name in enumerate(fp_adata.obs_names):
                if i < len(mean_scores):
                    fp_data[(str(group_name), tf_label)] = float(mean_scores[i])
        except Exception as e:
            print(f"  WARNING: Failed to load {fp_file}: {e}")

    print(f"  Loaded {len(fp_data)} (cell_type, TF) footprint scores")

    # Cross-modal correlation
    validation_records = []

    for tf, genes in target_genes.items():
        if tf not in gene_to_idx:
            continue

        tf_idx = gene_to_idx[tf]

        for gene in genes:
            if gene not in gene_to_idx:
                continue
            gene_idx = gene_to_idx[gene]

            # Collect per-cell-type values
            ct_fp_scores = []
            ct_tf_expr = []
            ct_gene_expr = []
            ct_labels = []

            for ct in cell_types:
                fp_key = (ct, tf)
                if fp_key not in fp_data:
                    # Try normalized match
                    ct_norm = normalize_celltype(ct)
                    for k in fp_data:
                        if normalize_celltype(k[0]) == ct_norm and k[1] == tf:
                            fp_key = k
                            break
                    else:
                        continue

                ct_fp_scores.append(fp_data[fp_key])
                ct_tf_expr.append(float(ct_mean_expr[ct][tf_idx]))
                ct_gene_expr.append(float(ct_mean_expr[ct][gene_idx]))
                ct_labels.append(ct)

            if len(ct_fp_scores) < 3:
                continue

            fp_arr = np.array(ct_fp_scores)
            tf_arr = np.array(ct_tf_expr)
            gene_arr = np.array(ct_gene_expr)

            # Spearman correlations
            corr_fp_tf, pval_fp_tf = stats.spearmanr(fp_arr, tf_arr)
            corr_fp_gene, pval_fp_gene = stats.spearmanr(fp_arr, gene_arr)
            corr_tf_gene, pval_tf_gene = stats.spearmanr(tf_arr, gene_arr)

            # Classify per recipe criteria
            if corr_fp_tf > 0.3 and corr_fp_gene > 0.3 and corr_tf_gene > 0.3:
                classification = 'activator'
            elif corr_fp_tf > 0.3 and corr_fp_gene < -0.3:
                classification = 'repressor'
            elif abs(corr_fp_tf) < 0.3 and abs(corr_fp_gene) > 0.3:
                classification = 'indirect'
            else:
                classification = 'unclassified'

            validation_records.append({
                'TF': tf,
                'gene': gene,
                'n_cell_types': len(ct_fp_scores),
                'corr_footprint_tfexpr': round(corr_fp_tf, 4),
                'pval_footprint_tfexpr': round(pval_fp_tf, 6),
                'corr_footprint_geneexpr': round(corr_fp_gene, 4),
                'pval_footprint_geneexpr': round(pval_fp_gene, 6),
                'corr_tfexpr_geneexpr': round(corr_tf_gene, 4),
                'pval_tfexpr_geneexpr': round(pval_tf_gene, 6),
                'classification': classification,
            })

    # Write validation table
    val_df = pd.DataFrame(validation_records) if validation_records else pd.DataFrame(
        columns=['TF', 'gene', 'n_cell_types', 'corr_footprint_tfexpr',
                 'pval_footprint_tfexpr', 'corr_footprint_geneexpr',
                 'pval_footprint_geneexpr', 'corr_tfexpr_geneexpr',
                 'pval_tfexpr_geneexpr', 'classification'])
    val_path = outdir / "cross_modal_validation.tsv"
    val_df.to_csv(val_path, sep='\t', index=False)
    print(f"\nWrote {len(val_df)} validation entries to {val_path}")

    # Classification summary
    if not val_df.empty:
        class_counts = val_df['classification'].value_counts()
        print("\nClassification summary:")
        for cls, count in class_counts.items():
            print(f"  {cls}: {count}")

        # Plot: classification pie chart
        try:
            fig, ax = plt.subplots(figsize=(8, 6))
            class_counts.plot.pie(ax=ax, autopct='%1.1f%%')
            ax.set_title("TF-Gene Classification Distribution")
            ax.set_ylabel("")
            plt.tight_layout()
            plt.savefig(plots_dir / "classification_distribution.png", dpi=150)
            plt.close()
        except Exception as e:
            print(f"  Plot failed: {e}")
            plt.close('all')

        # Plot: correlation scatter for top TFs
        unique_tfs = val_df['TF'].unique()[:10]
        for tf in unique_tfs:
            tf_rows = val_df[val_df['TF'] == tf]
            if len(tf_rows) < 3:
                continue
            try:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
                ax1.scatter(tf_rows['corr_footprint_tfexpr'],
                           tf_rows['corr_footprint_geneexpr'],
                           c=['green' if c == 'activator' else 'red' if c == 'repressor'
                              else 'gray' for c in tf_rows['classification']])
                ax1.axhline(0, color='k', ls='--', alpha=0.3)
                ax1.axvline(0, color='k', ls='--', alpha=0.3)
                ax1.set_xlabel("Corr(footprint, TF expr)")
                ax1.set_ylabel("Corr(footprint, gene expr)")
                ax1.set_title(f"{tf}: Cross-modal correlations")

                ax2.barh(tf_rows['gene'][:20], tf_rows['corr_footprint_geneexpr'][:20],
                        color=['green' if c == 'activator' else 'red' if c == 'repressor'
                               else 'gray' for c in tf_rows['classification'][:20]])
                ax2.set_xlabel("Corr(footprint, gene expr)")
                ax2.set_title(f"{tf}: Target gene correlations")
                plt.tight_layout()
                plt.savefig(plots_dir / f"{tf}_validation.png", dpi=150)
                plt.close()
            except Exception as e:
                print(f"  Plot failed for {tf}: {e}")
                plt.close('all')

    # Summary
    summary_lines = [
        "Cross-Modal Validation Summary",
        "=" * 40,
        f"Total TF-gene pairs tested: {len(val_df)}",
        f"Unique TFs: {val_df['TF'].nunique() if not val_df.empty else 0}",
        f"Unique target genes: {val_df['gene'].nunique() if not val_df.empty else 0}",
        f"Cell types in RNA: {len(cell_types)}",
        f"Footprint scores loaded: {len(fp_data)}",
    ]
    if not val_df.empty:
        summary_lines.append("")
        summary_lines.append("Classifications:")
        for cls, count in val_df['classification'].value_counts().items():
            summary_lines.append(f"  {cls}: {count}")

    summary_path = outdir / "validation_summary.txt"
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
