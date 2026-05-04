#!/usr/bin/env python3
"""
plot_tf_differential_volcano.py — Shi Fig 4B / 5A / 5C / 5E equivalent.

For each tf_differential_<celltype>_<trt>_vs_<ctrl>.csv produced by
differential_tf_accessibility.py, render a volcano plot of TF motif
log-fold-change vs -log10(adjusted p), with the top motifs labelled.

Inputs:
  --diff-tf-dir  results/differential_tf/
  --treatment, --control  condition labels (used to build expected filenames)

Outputs (in --outdir):
  <ct>_volcano.{pdf,png}                — one per cell type
  combined_volcano.{pdf,png}            — 4-panel grid for the Shi lineages
  volcano_manifest.json                 — file listing for downstream consumers
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Map Shi lineages → our broad classes (best analogue per lineage)
SHI_LINEAGE_MAP = {
    "EX (Excitatory)":   "IT-ET_Glut",
    "ASC (Astrocyte)":   "Astro-Epen",
    "MG (Microglia)":    "Immune",
    "ODC (Oligodendrocyte)": "OPC-Oligo",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--diff-tf-dir", required=True)
    p.add_argument("--treatment", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--outdir", default=".")
    p.add_argument("--fdr-cutoff", type=float, default=0.05)
    p.add_argument("--lfc-cutoff", type=float, default=0.5)
    p.add_argument("--top-labels", type=int, default=10)
    return p.parse_args()


def _clean_motif(s):
    return str(s).replace("\t", " ").replace("\n", " ").strip()


def _label_text(motif):
    parts = motif.split()
    if len(parts) >= 2 and parts[0].startswith("MA"):
        return parts[1]
    return parts[-1] if parts else motif


def render_volcano(ax, df, fdr, lfc, top_n, treatment, control, title):
    if df.empty:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
        return

    df = df.copy()
    df["motif"] = df["motif"].map(_clean_motif)
    df["neglog10p"] = -np.log10(df["pval_adj"].clip(lower=1e-300))
    df["sig_up"] = (df["pval_adj"] < fdr) & (df["logfoldchange"] >= lfc)
    df["sig_dn"] = (df["pval_adj"] < fdr) & (df["logfoldchange"] <= -lfc)
    df["sig"] = df["sig_up"] | df["sig_dn"]

    ax.scatter(df.loc[~df["sig"], "logfoldchange"],
               df.loc[~df["sig"], "neglog10p"],
               s=8, color="#bbbbbb", alpha=0.5, linewidth=0)
    ax.scatter(df.loc[df["sig_up"], "logfoldchange"],
               df.loc[df["sig_up"], "neglog10p"],
               s=14, color="#d62728", alpha=0.8, linewidth=0,
               label=f"up in {treatment}")
    ax.scatter(df.loc[df["sig_dn"], "logfoldchange"],
               df.loc[df["sig_dn"], "neglog10p"],
               s=14, color="#1f77b4", alpha=0.8, linewidth=0,
               label=f"up in {control}")

    ax.axhline(-np.log10(fdr), color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(lfc, color="gray", linestyle="--", linewidth=0.5)
    ax.axvline(-lfc, color="gray", linestyle="--", linewidth=0.5)

    sig = df[df["sig"]].copy()
    if not sig.empty:
        sig["score_rank"] = sig["logfoldchange"].abs() + sig["neglog10p"] / 50.0
        for _, row in sig.nlargest(top_n, "score_rank").iterrows():
            ax.annotate(_label_text(row["motif"]),
                       xy=(row["logfoldchange"], row["neglog10p"]),
                       xytext=(3, 3), textcoords="offset points", fontsize=7)

    ax.set_xlabel(f"log2 FC ({treatment} / {control})")
    ax.set_ylabel("-log10 adjusted P")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper left", fontsize=7, framealpha=0.7)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    suffix = f"_{args.treatment}_vs_{args.control}.csv"
    csvs = sorted(Path(args.diff_tf_dir).glob(f"tf_differential_*{suffix}"))
    if not csvs:
        raise SystemExit(f"[volcano] no matching CSVs in {args.diff_tf_dir}")

    rendered = []
    df_by_ct = {}
    for csv in csvs:
        ct = csv.name[len("tf_differential_"):-len(suffix)]
        df = pd.read_csv(csv)
        df_by_ct[ct] = df

        fig, ax = plt.subplots(figsize=(6, 5))
        render_volcano(ax, df, args.fdr_cutoff, args.lfc_cutoff, args.top_labels,
                      args.treatment, args.control,
                      f"{ct} — {args.treatment} vs {args.control}")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            out = outdir / f"{ct}_volcano.{ext}"
            fig.savefig(out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        rendered.append(ct)
        print(f"[volcano] {ct} — n_motifs={len(df)} "
              f"n_sig={int((df['pval_adj']<args.fdr_cutoff).sum())}", flush=True)

    # Combined 4-panel for Shi lineages
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (label, target_ct) in zip(axes.flat, SHI_LINEAGE_MAP.items()):
        df = df_by_ct.get(target_ct, pd.DataFrame())
        title = f"{label} ↔ {target_ct}"
        render_volcano(ax, df, args.fdr_cutoff, args.lfc_cutoff, args.top_labels,
                      args.treatment, args.control, title)
    fig.suptitle(f"Differential TF activity — {args.treatment} vs {args.control}",
                 fontsize=12, y=1.0)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"combined_volcano.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "rendered_celltypes": rendered,
        "shi_lineage_map": SHI_LINEAGE_MAP,
        "fdr_cutoff": args.fdr_cutoff,
        "lfc_cutoff": args.lfc_cutoff,
        "treatment": args.treatment,
        "control": args.control,
    }
    (outdir / "volcano_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[volcano] wrote {len(rendered)} per-ct volcanoes + combined panel", flush=True)


if __name__ == "__main__":
    main()
