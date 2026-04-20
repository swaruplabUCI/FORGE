#!/usr/bin/env python
"""
Aggregate per-fit MultiVI masking-sweep outputs into the canonical results.

Reads every `fit_frac*_seed*_rna.csv`, `fit_frac*_seed*_atac.csv`, and
`fit_frac*_seed*_integration.json` under --inputs (glob), then emits:
  <output_dir>/masking_sweep_results.csv
  <output_dir>/masking_sweep_summary.json
  <output_dir>/masking_degradation_curves.pdf
  <output_dir>/masking_celltype_{rna,atac}.pdf  (if cell_type present)
  <output_dir>/masking_integration_quality.pdf
"""

import argparse
import glob
import json
import os
import sys
import warnings

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multivi_masking_sweep import (
    plot_celltype_breakdown,
    plot_degradation_curves,
    plot_integration_quality,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate per-fit masking-sweep outputs")
    p.add_argument("--inputs", required=True, help="Glob for per-fit files, e.g. 'fits/*'")
    p.add_argument("--output_dir", default="masking_sweep", help="Output directory")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    rna_files = sorted(glob.glob(os.path.join(args.inputs, "fit_*_rna.csv")))
    atac_files = sorted(glob.glob(os.path.join(args.inputs, "fit_*_atac.csv")))
    int_files = sorted(glob.glob(os.path.join(args.inputs, "fit_*_integration.json")))

    if not rna_files or not atac_files or not int_files:
        raise SystemExit(
            f"No fit outputs found under {args.inputs}. "
            f"Expected fit_*_rna.csv / fit_*_atac.csv / fit_*_integration.json. "
            f"Got rna={len(rna_files)} atac={len(atac_files)} int={len(int_files)}"
        )

    print(f"Aggregating {len(rna_files)} RNA, {len(atac_files)} ATAC, {len(int_files)} integration files")

    all_metric_dfs = []
    for f in rna_files + atac_files:
        all_metric_dfs.append(pd.read_csv(f))
    all_results_df = pd.concat(all_metric_dfs, ignore_index=True)

    integration_metrics = []
    for f in int_files:
        with open(f) as fh:
            integration_metrics.append(json.load(fh))

    fractions = sorted({float(r["fraction"]) for r in integration_metrics})
    seeds = sorted({int(r["seed"]) for r in integration_metrics})

    results_csv = os.path.join(args.output_dir, "masking_sweep_results.csv")
    all_results_df.to_csv(results_csv, index=False)

    summary = {"fractions": fractions, "seeds": seeds}
    for frac in fractions:
        frac_df = all_results_df[all_results_df["fraction"] == frac]
        rna_df = frac_df[frac_df["modality"] == "rna"]
        atac_df = frac_df[frac_df["modality"] == "atac"]
        summary[f"frac_{frac}"] = {
            "rna_spearman_r_median": float(rna_df["spearman_r"].median()) if len(rna_df) else None,
            "rna_spearman_r_smoothed_median": (
                float(rna_df["spearman_r_smoothed"].median())
                if "spearman_r_smoothed" in rna_df.columns and len(rna_df)
                else None
            ),
            "atac_auprc_median": float(atac_df["auprc"].median()) if len(atac_df) else None,
            "atac_spearman_r_smoothed_median": (
                float(atac_df["spearman_r_smoothed"].median())
                if "spearman_r_smoothed" in atac_df.columns and len(atac_df)
                else None
            ),
        }
    summary["integration"] = integration_metrics

    summary_path = os.path.join(args.output_dir, "masking_sweep_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("\nGenerating plots...")
    plot_degradation_curves(all_results_df, args.output_dir)
    if "cell_type" in all_results_df.columns:
        plot_celltype_breakdown(all_results_df, args.output_dir)
    plot_integration_quality(integration_metrics, args.output_dir)

    print(f"\nAggregation complete. Outputs in {args.output_dir}/")


if __name__ == "__main__":
    main()
