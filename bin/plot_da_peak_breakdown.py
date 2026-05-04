#!/usr/bin/env python3
"""
plot_da_peak_breakdown.py — Shi Fig 2D equivalent.

Stacked bar of differential peak categories (Promoter / Exonic / Intronic / Distal)
per cell type, separately for gain and loss directions. Reads:
  - DA_peaks_<ct>__<trt>_vs_<ctrl>.csv (peak, log2FC, pval, padj)
  - peaks_annotated.tsv (peak_id, peakType) from annotate_peak_types.py

Outputs:
  da_peak_breakdown.{pdf,png}
  da_peak_breakdown.tsv  (long table: cell_type, direction, peakType, count, percent)
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PEAK_TYPES = ["Promoter", "Exonic", "Intronic", "Distal"]
PEAK_TYPE_COLORS = {
    "Promoter": "#d62728",
    "Exonic":   "#9467bd",
    "Intronic": "#2ca02c",
    "Distal":   "#1f77b4",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--diff-dir", required=True, help="Dir with DA_peaks_*.csv")
    p.add_argument("--peak-annotation", required=True, help="peaks_annotated.tsv")
    p.add_argument("--treatment", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--fdr-cutoff", type=float, default=0.05)
    p.add_argument("--lfc-cutoff", type=float, default=0.5)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    annot = pd.read_csv(args.peak_annotation, sep="\t")
    if "peakType" not in annot.columns:
        raise SystemExit(f"--peak-annotation missing 'peakType' column: {args.peak_annotation}")
    annot_map = dict(zip(annot["peak_id"].astype(str), annot["peakType"].astype(str)))

    suffix = f"__{args.treatment}_vs_{args.control}.csv"
    csvs = sorted(Path(args.diff_dir).glob(f"DA_peaks_*{suffix}"))
    if not csvs:
        raise SystemExit(f"[2D] no DA_peaks CSVs matching suffix '{suffix}' in {args.diff_dir}")

    rows = []
    for csv in csvs:
        ct = csv.name[len("DA_peaks_"):-len(suffix)]
        df = pd.read_csv(csv)
        # Some upstream variants name the columns differently; harmonize.
        rename = {}
        if "log2FC" in df.columns:
            rename["log2FC"] = "logfoldchange"
        if "padj" in df.columns:
            rename["padj"] = "pval_adj"
        df = df.rename(columns=rename)
        if "peak" not in df.columns:
            df = df.rename(columns={df.columns[0]: "peak"})

        df = df.dropna(subset=["pval_adj", "logfoldchange"])
        sig = df[(df["pval_adj"] < args.fdr_cutoff) &
                 (df["logfoldchange"].abs() >= args.lfc_cutoff)].copy()
        if sig.empty:
            print(f"[2D] {ct}: no significant peaks; skipping bar contribution", flush=True)
            continue
        sig["peakType"] = sig["peak"].astype(str).map(annot_map).fillna("Distal")
        sig["direction"] = np.where(sig["logfoldchange"] > 0,
                                    f"up in {args.treatment}",
                                    f"up in {args.control}")
        for direction, sub in sig.groupby("direction"):
            counts = sub["peakType"].value_counts()
            total = int(counts.sum())
            for pt in PEAK_TYPES:
                cnt = int(counts.get(pt, 0))
                rows.append({
                    "cell_type": ct,
                    "direction": direction,
                    "peakType": pt,
                    "count": cnt,
                    "percent": (cnt / total * 100) if total else 0.0,
                })
    if not rows:
        raise SystemExit("[2D] no significant peaks found across any cell type")

    long_df = pd.DataFrame(rows)
    long_df.to_csv(outdir / "da_peak_breakdown.tsv", sep="\t", index=False)

    cts = sorted(long_df["cell_type"].unique())
    directions = [f"up in {args.treatment}", f"up in {args.control}"]

    fig, axes = plt.subplots(1, 2, figsize=(max(10, len(cts) * 0.9), 5),
                             sharey=True)
    for ax, direction in zip(axes, directions):
        sub = long_df[long_df["direction"] == direction]
        wide = sub.pivot_table(index="cell_type", columns="peakType",
                               values="percent", aggfunc="sum",
                               fill_value=0).reindex(index=cts, columns=PEAK_TYPES,
                                                     fill_value=0)
        bottom = np.zeros(len(wide))
        for pt in PEAK_TYPES:
            ax.bar(wide.index, wide[pt], bottom=bottom,
                   color=PEAK_TYPE_COLORS[pt], label=pt, edgecolor="white",
                   linewidth=0.4)
            bottom = bottom + wide[pt].values
        ax.set_title(direction, fontsize=10)
        ax.set_ylabel("% of significant peaks")
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", rotation=45)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("right")
    axes[-1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
                    title="peakType", frameon=False)
    fig.suptitle(f"Differential peak biotype breakdown "
                 f"(FDR<{args.fdr_cutoff}, |log2FC|>{args.lfc_cutoff})", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"da_peak_breakdown.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[2D] wrote da_peak_breakdown.{{pdf,png,tsv}} for {len(cts)} cell types",
          flush=True)


if __name__ == "__main__":
    main()
