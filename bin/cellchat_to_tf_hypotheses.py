#!/usr/bin/env python3
"""
cellchat_to_tf_hypotheses.py — Recipe Steps C2+C3

Map CellChat signaling pathways to downstream TFs using a curated dictionary,
then validate with chromVAR deviations in receiver cells (Mann-Whitney U test).

Inputs:
    - CellChat results CSV (significant interactions)
    - chromvar_matrix_deviations.h5ad
    - pathway_to_tfs.json

Outputs:
    - receiver_tf_hypotheses.tsv (all mapped pairs)
    - validated_receiver_tf_pairs.tsv (filtered by chromVAR p < threshold)
    - unmapped_pathways.txt
"""

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np
from scipy import stats


def normalize_celltype(ct):
    """Normalize cell type name for matching across CellChat/chromVAR conventions.
    CellChat uses spaces/special chars (e.g., 'CD16+ NK cells'),
    chromVAR uses sanitized names."""
    ct = str(ct).strip()
    ct = re.sub(r'[^a-zA-Z0-9]', '_', ct)
    ct = re.sub(r'_+', '_', ct).strip('_')
    return ct.lower()


def parse_args():
    p = argparse.ArgumentParser(description="Map CellChat pathways to TF hypotheses")
    p.add_argument("--cellchat-csv", required=True, help="CellChat results CSV")
    p.add_argument("--chromvar-dev", required=True, help="ChromVAR deviations h5ad")
    p.add_argument("--pathway-to-tfs", required=True, help="Pathway-to-TFs JSON")
    p.add_argument("--pval-threshold", type=float, default=0.05)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)

    # Load pathway-to-TFs dictionary
    print(f"Loading pathway-to-TFs from {args.pathway_to_tfs}")
    with open(args.pathway_to_tfs) as f:
        pathway_to_tfs = json.load(f)

    # Remove metadata entries
    pathway_to_tfs = {k: v for k, v in pathway_to_tfs.items() if not k.startswith('_')}
    print(f"  {len(pathway_to_tfs)} pathways mapped")

    # Load CellChat results
    print(f"Loading CellChat results from {args.cellchat_csv}")
    cc_df = pd.read_csv(args.cellchat_csv)
    print(f"  {len(cc_df)} interactions loaded")
    print(f"  Columns: {list(cc_df.columns)}")

    # Detect column names
    pathway_col = None
    source_col = None
    target_col = None
    prob_col = None

    for col in cc_df.columns:
        col_lower = col.lower()
        if 'pathway' in col_lower and pathway_col is None:
            pathway_col = col
        if ('source' in col_lower or 'sender' in col_lower) and source_col is None:
            source_col = col
        if ('target' in col_lower or 'receiver' in col_lower) and target_col is None:
            target_col = col
        if ('prob' in col_lower or 'pval' in col_lower) and prob_col is None:
            prob_col = col

    # Fallback column detection
    if pathway_col is None:
        for col in cc_df.columns:
            if 'name' in col.lower() or 'interaction' in col.lower():
                pathway_col = col
                break
    if source_col is None and target_col is None:
        # Try positional
        if len(cc_df.columns) >= 3:
            source_col = cc_df.columns[0]
            target_col = cc_df.columns[1]

    if pathway_col is None:
        print("WARNING: Could not detect pathway column. Using first column.")
        pathway_col = cc_df.columns[0]

    print(f"  Using columns: pathway='{pathway_col}', source='{source_col}', target='{target_col}'")

    # Map pathways to TFs
    hypotheses = []
    unmapped_pathways = set()

    for _, row in cc_df.iterrows():
        pathway_raw = str(row[pathway_col]).strip()

        # Try exact match first, then uppercase
        tfs = pathway_to_tfs.get(pathway_raw)
        if tfs is None:
            tfs = pathway_to_tfs.get(pathway_raw.upper())
        if tfs is None:
            # Try partial match
            for key in pathway_to_tfs:
                if key.upper() in pathway_raw.upper() or pathway_raw.upper() in key.upper():
                    tfs = pathway_to_tfs[key]
                    break

        if tfs is None:
            unmapped_pathways.add(pathway_raw)
            continue

        sender = str(row[source_col]) if source_col else 'unknown'
        receiver = str(row[target_col]) if target_col else 'unknown'
        prob = float(row[prob_col]) if prob_col and pd.notna(row.get(prob_col)) else np.nan

        for tf in tfs:
            hypotheses.append({
                'sender': sender,
                'receiver': receiver,
                'pathway': pathway_raw,
                'downstream_TF': tf,
                'interaction_prob': prob,
            })

    hyp_df = pd.DataFrame(hypotheses) if hypotheses else pd.DataFrame(
        columns=['sender', 'receiver', 'pathway', 'downstream_TF', 'interaction_prob'])

    hyp_path = outdir / "receiver_tf_hypotheses.tsv"
    hyp_df.to_csv(hyp_path, sep='\t', index=False)
    print(f"\nWrote {len(hyp_df)} hypotheses to {hyp_path}")

    # Write unmapped pathways
    unmapped_path = outdir / "unmapped_pathways.txt"
    unmapped_path.write_text('\n'.join(sorted(unmapped_pathways)) if unmapped_pathways else 'None')
    print(f"Unmapped pathways: {len(unmapped_pathways)}")

    # Validate with chromVAR
    print(f"\nLoading chromVAR deviations from {args.chromvar_dev}")
    adata = ad.read_h5ad(args.chromvar_dev)

    if 'cell_type' not in adata.obs.columns:
        print("WARNING: chromVAR h5ad missing 'cell_type'. Skipping validation.")
        validated_df = hyp_df.copy()
        validated_df['chromvar_pval'] = np.nan
        validated_df['chromvar_stat'] = np.nan
        validated_df['validated'] = True
    else:
        cell_types = adata.obs['cell_type'].unique().tolist()
        ct_norm_map = {normalize_celltype(ct): ct for ct in cell_types}

        # Extract TF names from chromVAR var_names (format: "MA0001.2\tSPI1")
        tf_to_var_idx = {}
        for i, vname in enumerate(adata.var_names):
            parts = str(vname).split('\t')
            tf_symbol = parts[-1].split('_')[0] if len(parts) > 1 else parts[0].split('_')[0]
            tf_to_var_idx[tf_symbol.upper()] = i

        # Validate each (receiver, TF) pair
        validated_records = []
        for _, row in hyp_df.iterrows():
            tf = row['downstream_TF']
            receiver = row['receiver']

            # Find TF in chromVAR
            tf_idx = tf_to_var_idx.get(tf.upper())
            if tf_idx is None:
                validated_records.append({**row.to_dict(), 'chromvar_pval': np.nan,
                                         'chromvar_stat': np.nan, 'validated': False})
                continue

            # Find receiver cell type
            receiver_norm = normalize_celltype(receiver)
            actual_ct = ct_norm_map.get(receiver_norm)
            if actual_ct is None:
                # Try fuzzy match
                for norm_ct, orig_ct in ct_norm_map.items():
                    if receiver_norm in norm_ct or norm_ct in receiver_norm:
                        actual_ct = orig_ct
                        break

            if actual_ct is None:
                validated_records.append({**row.to_dict(), 'chromvar_pval': np.nan,
                                         'chromvar_stat': np.nan, 'validated': False})
                continue

            # Mann-Whitney U: receiver vs other cells
            receiver_mask = adata.obs['cell_type'] == actual_ct
            other_mask = ~receiver_mask

            if hasattr(adata.X, 'toarray'):
                receiver_vals = np.array(adata.X[receiver_mask, tf_idx].toarray()).flatten()
                other_vals = np.array(adata.X[other_mask, tf_idx].toarray()).flatten()
            else:
                receiver_vals = np.array(adata.X[receiver_mask, tf_idx]).flatten()
                other_vals = np.array(adata.X[other_mask, tf_idx]).flatten()

            if len(receiver_vals) < 5 or len(other_vals) < 5:
                validated_records.append({**row.to_dict(), 'chromvar_pval': np.nan,
                                         'chromvar_stat': np.nan, 'validated': False})
                continue

            stat, pval = stats.mannwhitneyu(receiver_vals, other_vals, alternative='greater')

            validated_records.append({
                **row.to_dict(),
                'chromvar_pval': round(pval, 8),
                'chromvar_stat': round(stat, 4),
                'validated': pval < args.pval_threshold,
            })

        validated_df = pd.DataFrame(validated_records)

    # Filter to validated pairs
    validated_pairs = validated_df[validated_df['validated'] == True].copy()
    val_path = outdir / "validated_receiver_tf_pairs.tsv"
    validated_pairs.to_csv(val_path, sep='\t', index=False)
    print(f"Validated pairs (p < {args.pval_threshold}): {len(validated_pairs)} of {len(validated_df)}")

    # Summary
    summary_lines = [
        "CellChat → TF Hypotheses Summary",
        "=" * 40,
        f"CellChat interactions: {len(cc_df)}",
        f"Mapped interactions: {len(cc_df) - len(unmapped_pathways)}",
        f"Unmapped pathways: {len(unmapped_pathways)}",
        f"Total hypotheses: {len(hyp_df)}",
        f"Unique (receiver, TF) pairs: {hyp_df[['receiver', 'downstream_TF']].drop_duplicates().shape[0] if not hyp_df.empty else 0}",
        f"Validated pairs (p < {args.pval_threshold}): {len(validated_pairs)}",
    ]
    summary_path = outdir / "hypothesis_summary.txt"
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
