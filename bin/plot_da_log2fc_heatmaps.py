#!/usr/bin/env python3
"""
plot_da_log2fc_heatmaps.py — Shi Fig 2E equivalent.

Per-celltype heatmap of differential peak log2FC, with rows split into
Promoter and Distal blocks. Distal rows are restricted to peaks that
participate in a CCAN promoter-distal link (Cicero output). The right margin
is annotated with the most-frequent nearest gene per row block, replacing
Shi's GREAT GO sidebar.

Inputs:
  --diff-dir          DA_peaks_*.csv directory
  --peak-annotation   peaks_annotated.tsv (peakType + nearestGene from annotate_peak_types)
  --gene-links        ccan_enhancer_gene_links.tsv (defines Distal subset)

Outputs:
  <ct>_log2fc.{pdf,png}                — one per cell type
  log2fc_heatmaps_summary.tsv          — per-ct row counts, top genes
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--diff-dir", required=True)
    p.add_argument("--peak-annotation", required=True)
    p.add_argument("--gene-links", required=True,
                   help="ccan_enhancer_gene_links.tsv — defines distal-CCAN peaks")
    p.add_argument("--treatment", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--fdr-cutoff", type=float, default=0.05)
    p.add_argument("--lfc-cutoff", type=float, default=0.5)
    p.add_argument("--top-genes-per-block", type=int, default=8)
    p.add_argument("--max-rows-per-block", type=int, default=80)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    annot = pd.read_csv(args.peak_annotation, sep="\t")
    annot_idx = annot.set_index("peak_id")
    if "peakType" not in annot.columns or "nearestGene" not in annot.columns:
        raise SystemExit("--peak-annotation needs peakType and nearestGene columns")

    links = pd.read_csv(args.gene_links, sep="\t")
    distal_peaks = set((links["chr"].astype(str) + ":" +
                       links["start"].astype(str) + "-" +
                       links["end"].astype(str)).unique())

    suffix = f"__{args.treatment}_vs_{args.control}.csv"
    csvs = sorted(Path(args.diff_dir).glob(f"DA_peaks_*{suffix}"))
    if not csvs:
        raise SystemExit(f"[2E] no DA CSVs in {args.diff_dir}")

    summary_rows = []
    for csv in csvs:
        ct = csv.name[len("DA_peaks_"):-len(suffix)]
        df = pd.read_csv(csv)
        df = df.rename(columns={"log2FC": "logfoldchange", "padj": "pval_adj"})
        if "peak" not in df.columns:
            df = df.rename(columns={df.columns[0]: "peak"})
        df = df.dropna(subset=["pval_adj", "logfoldchange"])
        sig = df[(df["pval_adj"] < args.fdr_cutoff) &
                 (df["logfoldchange"].abs() >= args.lfc_cutoff)].copy()
        if sig.empty:
            continue

        sig = sig.merge(
            annot_idx[["peakType", "nearestGene"]],
            left_on="peak", right_index=True, how="left")
        sig["peakType"] = sig["peakType"].fillna("Distal")

        prom = sig[sig["peakType"] == "Promoter"].copy()
        dist = sig[sig["peak"].isin(distal_peaks)].copy()

        if prom.empty and dist.empty:
            continue

        # Sort each block by signed log2FC, cap row counts so the figure stays
        # legible even when n_sig is in the thousands.
        def _cap(b):
            b = b.sort_values("logfoldchange")
            n = len(b)
            if n > args.max_rows_per_block:
                k = args.max_rows_per_block // 2
                b = pd.concat([b.head(k), b.tail(k)])
            return b
        prom_disp = _cap(prom)
        dist_disp = _cap(dist)

        n_total = max(len(prom_disp) + len(dist_disp), 1)
        height_inches = max(4, n_total * 0.06)
        fig, ax = plt.subplots(figsize=(4, height_inches))
        rows = pd.concat([prom_disp.assign(_block="Promoter"),
                         dist_disp.assign(_block="Distal")])
        if rows.empty:
            plt.close(fig)
            continue

        vmax = max(abs(rows["logfoldchange"]).max(), 1e-3)
        mat = rows[["logfoldchange"]].values
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks([0])
        ax.set_xticklabels([f"{args.treatment}\nvs {args.control}"], rotation=0)
        ax.set_yticks([])

        # Block separator
        n_prom = len(prom_disp)
        if 0 < n_prom < n_total:
            ax.axhline(n_prom - 0.5, color="black", linewidth=0.7)
        ax.text(-0.6, n_prom / 2, "Promoter", rotation=90, va="center", ha="right",
                fontsize=8)
        ax.text(-0.6, n_prom + (n_total - n_prom) / 2, "Distal", rotation=90,
                va="center", ha="right", fontsize=8)

        cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
        cbar.set_label("log2 FC")

        prom_top = (prom["nearestGene"].value_counts().head(args.top_genes_per_block)
                    .index.tolist() if not prom.empty else [])
        dist_top = (dist["nearestGene"].value_counts().head(args.top_genes_per_block)
                    .index.tolist() if not dist.empty else [])
        gene_text = "Top promoter genes:\n  " + ", ".join(prom_top)
        gene_text += "\n\nTop distal-CCAN genes:\n  " + ", ".join(dist_top)
        ax.set_title(f"{ct}\nFDR<{args.fdr_cutoff}  |log2FC|>{args.lfc_cutoff}",
                     fontsize=9)
        fig.text(1.05, 0.5, gene_text, transform=ax.transAxes, fontsize=7,
                va="center")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{ct}_log2fc.{ext}", dpi=200, bbox_inches="tight")
        plt.close(fig)

        summary_rows.append({
            "cell_type": ct,
            "n_promoter": int(len(prom)),
            "n_distal": int(len(dist)),
            "promoter_top_genes": ",".join(prom_top),
            "distal_top_genes": ",".join(dist_top),
        })
        print(f"[2E] {ct} — promoter rows={len(prom)} distal rows={len(dist)}",
              flush=True)

    if not summary_rows:
        raise SystemExit("[2E] no significant peaks across any cell type")

    pd.DataFrame(summary_rows).to_csv(
        outdir / "log2fc_heatmaps_summary.tsv", sep="\t", index=False)
    print(f"[2E] wrote {len(summary_rows)} per-ct heatmaps", flush=True)


if __name__ == "__main__":
    main()
