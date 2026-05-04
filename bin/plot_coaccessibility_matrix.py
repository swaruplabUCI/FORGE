#!/usr/bin/env python3
"""
plot_coaccessibility_matrix.py — Shi Fig 2C cleanup.

Shi 2C shows a cell-type × condition correlation matrix of Cicero
co-accessibility. Our stratified Cicero is run per condition (not per
cell type), so the most-faithful condition-level matrix is degenerate
(2×2). We instead produce a richer comparison that exposes the same
structural question (does the regulatory network differ between
conditions?) using the data we have:

  - Pearson r of shared peak-pair coaccess between control and treatment
  - Hexbin density scatter of paired coaccess values (the existing plot)
  - r and shared-pair counts split by peak-pair biotype using the
    peak annotation file (Promoter-Distal, Distal-Distal, etc.)
  - r split by genomic distance bins (0-50 kb, 50-200 kb, 200-500 kb)

Inputs:
  --connections-ctrl     Cicero connections TSV (control)
  --connections-trt      Cicero connections TSV (treatment)
  --peak-annotation      peaks_annotated.tsv (peakType column required)

Outputs (in --outdir):
  coacc_correlation_matrix.{pdf,png}     — multi-panel summary figure
  coacc_correlation_summary.tsv          — long table of r per stratum
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


DISTANCE_BINS = [(0, 50_000), (50_000, 200_000), (200_000, 500_000)]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--connections-ctrl", required=True)
    p.add_argument("--connections-trt", required=True)
    p.add_argument("--peak-annotation", required=True)
    p.add_argument("--control-label", default="control")
    p.add_argument("--treatment-label", default="treatment")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def load_connections(path):
    df = pd.read_csv(path, sep="\t")
    cols = list(df.columns)
    df = df.rename(columns={cols[0]: "Peak1", cols[1]: "Peak2", cols[2]: "coaccess"})
    df["coaccess"] = pd.to_numeric(df["coaccess"], errors="coerce")
    return df.dropna(subset=["coaccess"])


def _peak_midpoint(peak_id):
    try:
        chrom, span = peak_id.split(":")
        start, end = span.split("-")
        return chrom, (int(start) + int(end)) // 2
    except Exception:
        return None, None


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ctrl = load_connections(args.connections_ctrl)
    trt = load_connections(args.connections_trt)

    # canonical pair key
    for d in (ctrl, trt):
        a = np.minimum(d["Peak1"], d["Peak2"])
        b = np.maximum(d["Peak1"], d["Peak2"])
        d["pair_key"] = a + "||" + b

    ctrl_d = ctrl.groupby("pair_key", as_index=False)["coaccess"].mean()
    trt_d = trt.groupby("pair_key", as_index=False)["coaccess"].mean()
    merged = ctrl_d.merge(trt_d, on="pair_key", suffixes=("_ctrl", "_trt"), how="inner")
    print(f"[2C] shared peak pairs: {len(merged):,}", flush=True)

    annot = pd.read_csv(args.peak_annotation, sep="\t")
    pt_map = dict(zip(annot["peak_id"].astype(str), annot["peakType"].astype(str)))

    pairs = merged["pair_key"].str.split("||", expand=True, regex=False)
    pairs.columns = ["Peak1", "Peak2"]
    merged["pt1"] = pairs["Peak1"].map(pt_map).fillna("Distal")
    merged["pt2"] = pairs["Peak2"].map(pt_map).fillna("Distal")
    merged["pair_type"] = merged.apply(
        lambda r: "-".join(sorted([r["pt1"], r["pt2"]])), axis=1)

    # Distance bin (skip pairs from different chromosomes)
    chr1, mid1 = zip(*pairs["Peak1"].map(_peak_midpoint))
    chr2, mid2 = zip(*pairs["Peak2"].map(_peak_midpoint))
    chr1 = pd.Series(chr1); chr2 = pd.Series(chr2)
    mid1 = pd.Series(mid1); mid2 = pd.Series(mid2)
    same_chr = (chr1 == chr2)
    distance = (mid2 - mid1).abs().where(same_chr, np.nan)
    merged["distance"] = distance.values

    # Stratified r values
    rows = []
    rows.append(_record("all", "all", merged))
    for pt, sub in merged.groupby("pair_type"):
        if len(sub) >= 50:
            rows.append(_record("pair_type", pt, sub))
    for lo, hi in DISTANCE_BINS:
        sub = merged[(merged["distance"] >= lo) & (merged["distance"] < hi)]
        if len(sub) >= 50:
            rows.append(_record("distance_kb", f"{lo//1000}-{hi//1000}", sub))

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "coacc_correlation_summary.tsv", sep="\t", index=False)

    # Render multi-panel figure
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3)
    ax_hex = fig.add_subplot(gs[0, 0])
    ax_mat = fig.add_subplot(gs[0, 1])
    ax_pt  = fig.add_subplot(gs[1, 0])
    ax_dist = fig.add_subplot(gs[1, 1])
    ax_sum = fig.add_subplot(gs[:, 2])

    # 1. Hexbin
    if len(merged) > 10:
        hb = ax_hex.hexbin(merged["coaccess_ctrl"], merged["coaccess_trt"],
                          gridsize=50, cmap="YlOrRd", mincnt=1)
        plt.colorbar(hb, ax=ax_hex, label="count")
    r_all = summary.loc[summary["stratum"] == "all", "pearson_r"].iloc[0]
    ax_hex.set_xlabel(f"{args.control_label} coaccess")
    ax_hex.set_ylabel(f"{args.treatment_label} coaccess")
    ax_hex.set_title(f"All pairs — r = {r_all:.3f}", fontsize=10)
    lim = (min(ax_hex.get_xlim()[0], ax_hex.get_ylim()[0]),
           max(ax_hex.get_xlim()[1], ax_hex.get_ylim()[1]))
    ax_hex.plot(lim, lim, "--", color="gray", linewidth=0.6)

    # 2. 2x2 condition matrix
    M = np.array([[1.0, r_all], [r_all, 1.0]])
    im = ax_mat.imshow(M, vmin=-1, vmax=1, cmap="RdBu_r")
    ax_mat.set_xticks([0, 1]); ax_mat.set_yticks([0, 1])
    ax_mat.set_xticklabels([args.control_label, args.treatment_label])
    ax_mat.set_yticklabels([args.control_label, args.treatment_label])
    for i in range(2):
        for j in range(2):
            ax_mat.text(j, i, f"{M[i, j]:.3f}", ha="center", va="center",
                       color="black", fontsize=10)
    ax_mat.set_title("Condition correlation matrix")
    plt.colorbar(im, ax=ax_mat, fraction=0.05)

    # 3. Per-pair-type bar
    pt_rows = summary[summary["stratum"] == "pair_type"]
    if not pt_rows.empty:
        ax_pt.bar(pt_rows["value"], pt_rows["pearson_r"], color="#4c72b0")
        ax_pt.set_ylabel("Pearson r")
        ax_pt.set_title("Correlation by peak-pair biotype")
        ax_pt.tick_params(axis="x", rotation=30)
        for tl in ax_pt.get_xticklabels():
            tl.set_horizontalalignment("right")
        ax_pt.axhline(0, color="gray", linewidth=0.5)
    else:
        ax_pt.set_axis_off()

    # 4. Per-distance bar
    dist_rows = summary[summary["stratum"] == "distance_kb"]
    if not dist_rows.empty:
        ax_dist.bar(dist_rows["value"], dist_rows["pearson_r"], color="#55a868")
        ax_dist.set_ylabel("Pearson r")
        ax_dist.set_title("Correlation by genomic distance (kb)")
        ax_dist.axhline(0, color="gray", linewidth=0.5)
    else:
        ax_dist.set_axis_off()

    # 5. Summary table
    ax_sum.axis("off")
    table_rows = [["stratum", "value", "n_pairs", "Pearson r", "p"]]
    for _, r in summary.iterrows():
        table_rows.append([r["stratum"], str(r["value"]), f"{int(r['n_pairs']):,}",
                          f"{r['pearson_r']:.3f}", f"{r['pval']:.1e}"])
    tbl = ax_sum.table(cellText=table_rows[1:], colLabels=table_rows[0],
                      loc="center", cellLoc="left")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8)
    tbl.scale(1, 1.4)

    fig.suptitle(f"Cicero co-accessibility — {args.control_label} vs {args.treatment_label}",
                 fontsize=12)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"coacc_correlation_matrix.{ext}", dpi=200,
                   bbox_inches="tight")
    plt.close(fig)
    print(f"[2C] wrote multi-panel correlation summary "
          f"(n_strata={len(summary)})", flush=True)


def _record(stratum, value, sub):
    if len(sub) < 2:
        return {"stratum": stratum, "value": value, "n_pairs": len(sub),
                "pearson_r": np.nan, "pval": np.nan}
    r, p = stats.pearsonr(sub["coaccess_ctrl"], sub["coaccess_trt"])
    return {"stratum": stratum, "value": value, "n_pairs": int(len(sub)),
            "pearson_r": float(r), "pval": float(p)}


if __name__ == "__main__":
    main()
