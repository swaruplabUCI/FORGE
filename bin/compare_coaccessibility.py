#!/usr/bin/env python3
"""
compare_coaccessibility.py — Phase 4C

Compare co-accessibility scores between two conditions (Shi et al. methodology).
Mirrors the comparative analysis in PiD_AD_putativeEnhancer.Rmd:
  - Pearson correlation of shared peak-pair co-accessibility scores
  - Hexbin scatter plot
  - Identification of condition-specific enhancer-gene links

Inputs:
    - Cicero connections from control condition
    - Cicero connections from treatment condition

Outputs:
    - coaccessibility_correlation.tsv (per-peak-pair comparison)
    - condition_specific_enhancers.bed (differential enhancers)
    - coaccessibility_comparison.png (hexbin scatter)
    - coaccessibility_summary.txt
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


def parse_args():
    p = argparse.ArgumentParser(description="Compare co-accessibility between conditions")
    p.add_argument("--connections-ctrl", required=True, help="Control Cicero connections TSV")
    p.add_argument("--connections-trt", required=True, help="Treatment Cicero connections TSV")
    p.add_argument("--control-label", default="Control", help="Control condition label")
    p.add_argument("--treatment-label", default="Treatment", help="Treatment condition label")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def load_connections(path):
    """Load Cicero connections and standardize column names."""
    df = pd.read_csv(path, sep='\t')

    # Detect columns
    cols = df.columns.tolist()
    peak1_col = cols[0]
    peak2_col = cols[1]
    coaccess_col = cols[2]

    df = df.rename(columns={
        peak1_col: 'Peak1', peak2_col: 'Peak2', coaccess_col: 'coaccess'
    })
    df['coaccess'] = pd.to_numeric(df['coaccess'], errors='coerce')
    df = df.dropna(subset=['coaccess'])

    # Create canonical pair key (sorted so Peak1-Peak2 == Peak2-Peak1)
    df['pair_key'] = df.apply(
        lambda r: tuple(sorted([str(r['Peak1']), str(r['Peak2'])])), axis=1)
    df['pair_key'] = df['pair_key'].map(lambda x: f"{x[0]}||{x[1]}")

    return df


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading control connections: {args.connections_ctrl}")
    ctrl = load_connections(args.connections_ctrl)
    print(f"  {len(ctrl)} connections")

    print(f"Loading treatment connections: {args.connections_trt}")
    trt = load_connections(args.connections_trt)
    print(f"  {len(trt)} connections")

    # Merge on shared peak pairs
    ctrl_dedup = ctrl.groupby('pair_key')['coaccess'].mean().reset_index()
    ctrl_dedup.columns = ['pair_key', 'coaccess_ctrl']

    trt_dedup = trt.groupby('pair_key')['coaccess'].mean().reset_index()
    trt_dedup.columns = ['pair_key', 'coaccess_trt']

    merged = ctrl_dedup.merge(trt_dedup, on='pair_key', how='outer')
    merged['coaccess_ctrl'] = merged['coaccess_ctrl'].fillna(0)
    merged['coaccess_trt'] = merged['coaccess_trt'].fillna(0)

    n_shared = ((merged['coaccess_ctrl'] != 0) & (merged['coaccess_trt'] != 0)).sum()
    n_ctrl_only = ((merged['coaccess_ctrl'] != 0) & (merged['coaccess_trt'] == 0)).sum()
    n_trt_only = ((merged['coaccess_ctrl'] == 0) & (merged['coaccess_trt'] != 0)).sum()

    print(f"\n  Shared peak pairs: {n_shared}")
    print(f"  {args.control_label}-only pairs: {n_ctrl_only}")
    print(f"  {args.treatment_label}-only pairs: {n_trt_only}")

    # Compute delta
    merged['delta_coaccess'] = merged['coaccess_trt'] - merged['coaccess_ctrl']

    # Split pair_key back to Peak1, Peak2
    merged[['Peak1', 'Peak2']] = merged['pair_key'].str.split('||', expand=True, regex=False)

    # Pearson correlation on shared pairs with positive co-accessibility
    shared_mask = (merged['coaccess_ctrl'] > 0) & (merged['coaccess_trt'] > 0)
    shared = merged[shared_mask]

    if len(shared) > 2:
        r, p_val = stats.pearsonr(shared['coaccess_ctrl'], shared['coaccess_trt'])
        print(f"\n  Pearson r = {r:.4f}, p = {p_val:.2e} (n={len(shared)} shared pairs)")
    else:
        r, p_val = np.nan, np.nan
        print("  Not enough shared pairs for correlation")

    # Write correlation TSV
    out_cols = ['Peak1', 'Peak2', 'coaccess_ctrl', 'coaccess_trt', 'delta_coaccess']
    corr_path = outdir / 'coaccessibility_correlation.tsv'
    merged[out_cols].to_csv(corr_path, sep='\t', index=False)
    print(f"  Wrote {len(merged)} pair comparisons to {corr_path}")

    # Identify condition-specific enhancers (top differential pairs)
    # Use absolute delta > 0.1 and at least one positive coaccess
    diff_mask = (merged['delta_coaccess'].abs() > 0.1) & (
        (merged['coaccess_ctrl'] > 0) | (merged['coaccess_trt'] > 0))
    diff_pairs = merged[diff_mask].copy()
    diff_pairs['condition'] = np.where(
        diff_pairs['delta_coaccess'] > 0, args.treatment_label, args.control_label)

    # Extract unique enhancer peaks from differential pairs
    enh_records = []
    for _, row in diff_pairs.iterrows():
        for peak_str in [row['Peak1'], row['Peak2']]:
            parts = str(peak_str).replace('-', '_').split('_')
            if len(parts) >= 3:
                try:
                    enh_records.append({
                        'chr': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'condition': row['condition'],
                        'delta': row['delta_coaccess'],
                    })
                except ValueError:
                    continue

    if enh_records:
        enh_df = pd.DataFrame(enh_records).drop_duplicates(subset=['chr', 'start', 'end'])
        enh_df = enh_df.sort_values(['chr', 'start'])
    else:
        enh_df = pd.DataFrame(columns=['chr', 'start', 'end', 'condition', 'delta'])

    enh_path = outdir / 'condition_specific_enhancers.bed'
    enh_df.to_csv(enh_path, sep='\t', index=False, header=False)
    print(f"  Wrote {len(enh_df)} condition-specific enhancer peaks to {enh_path}")

    # Hexbin scatter plot (Shi et al. style)
    fig, ax = plt.subplots(figsize=(8, 8))

    if len(shared) > 10:
        hb = ax.hexbin(
            shared['coaccess_ctrl'], shared['coaccess_trt'],
            gridsize=50, cmap='YlOrRd', mincnt=1)
        plt.colorbar(hb, ax=ax, label='Count')
    else:
        ax.scatter(merged['coaccess_ctrl'], merged['coaccess_trt'],
                   alpha=0.3, s=5, color='steelblue')

    # Diagonal line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, '--', color='gray', alpha=0.5, linewidth=1)

    ax.set_xlabel(f'{args.control_label} co-accessibility', fontsize=12)
    ax.set_ylabel(f'{args.treatment_label} co-accessibility', fontsize=12)
    ax.set_title(
        f'Co-accessibility: {args.control_label} vs {args.treatment_label}\n'
        f'Pearson r = {r:.3f} (n = {len(shared):,} shared pairs)',
        fontsize=13)
    ax.set_aspect('equal')

    plt.tight_layout()
    plot_path = outdir / 'coaccessibility_comparison.png'
    plt.savefig(plot_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved scatter plot: {plot_path}")

    # Summary
    summary_lines = [
        "Co-Accessibility Comparison Summary",
        "=" * 45,
        f"Control: {args.control_label}",
        f"Treatment: {args.treatment_label}",
        f"",
        f"Control connections: {len(ctrl)}",
        f"Treatment connections: {len(trt)}",
        f"Total unique peak pairs: {len(merged)}",
        f"Shared pairs (both > 0): {n_shared}",
        f"{args.control_label}-only pairs: {n_ctrl_only}",
        f"{args.treatment_label}-only pairs: {n_trt_only}",
        f"",
        f"Pearson r: {r:.4f}",
        f"Pearson p-value: {p_val:.2e}",
        f"",
        f"Condition-specific enhancer peaks: {len(enh_df)}",
        f"  {args.control_label}-enriched: {(enh_df['condition'] == args.control_label).sum() if not enh_df.empty else 0}",
        f"  {args.treatment_label}-enriched: {(enh_df['condition'] == args.treatment_label).sum() if not enh_df.empty else 0}",
    ]

    summary_path = outdir / 'coaccessibility_summary.txt'
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
