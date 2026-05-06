#!/usr/bin/env python3
"""Per-(cell_type, gene) promoter MSFP heatmap + JASPAR motif overlay.

Reads a per-condition fp h5ad written by run_promoter_msfp_per_condition.py
(obs["name"]=condition labels, obsm[<chrom:start-end>]=(cond, mode, pos),
uns["scales"]=modes), clips scales to <= --max-scale (TF-binding band, not
nucleosome), runs an inline scprinter motif scan over the same window, and
stacks per-condition heatmaps above a motif-hit track.

Adapted from a one-off SREBF1 oligodendrocyte template (SDas_nf 2026-05-05);
the rendering recipe is preserved, the wiring (data-driven TF list, cycling
palette, parametrized scale ceiling) is FORGE-native.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scprinter as scp


_DEFAULT_PALETTE = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd",
                    "#ff7f0e", "#17becf", "#8c564b", "#e377c2"]


def parse_locus(locus_str):
    chrom, rng = locus_str.split(":")
    s, e = rng.split("-")
    return chrom, int(s), int(e)


def parse_tf_colors(tfs, tf_colors_arg):
    out = {}
    if tf_colors_arg:
        for token in tf_colors_arg.split(","):
            if "=" in token:
                k, v = token.split("=", 1)
                out[k.strip()] = v.strip()
    for i, tf in enumerate(tfs):
        out.setdefault(tf, _DEFAULT_PALETTE[i % len(_DEFAULT_PALETTE)])
    return out


def scan_motifs_in_window(chrom, start, end, pfm_path, cache_dir, genome, tfs):
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    genome_obj = getattr(scp.genome, genome)
    scanner = scp.motifs.Motifs(
        ref_path_motif=pfm_path,
        ref_path_fa=genome_obj.fetch_fa(),
        bg="even",
        n_jobs=1,
        motif_name_func=lambda x: x.split("\t")[-1],
    )
    scanner.prep_scanner(tf_genes=list(tfs))
    return scanner.scan_motif([(chrom, start, end)], verbose=False, split_tfs=True)


def parse_hits(hits, chrom, start, end, score_threshold=None):
    out = []
    sample_printed = False
    for h in hits:
        if not sample_printed:
            print(f"  raw hit sample: {h!r}", flush=True)
            sample_printed = True
        try:
            hit_chr = h[0]
            peak_s = int(h[1])
            tf_name = h[4]
            score = float(h[5])
            strand = int(h[6])
            motif_rel_s = int(h[7])
            motif_rel_e = int(h[8])
            if hit_chr != chrom:
                continue
            abs_s = peak_s + motif_rel_s
            abs_e = peak_s + motif_rel_e
            if abs_s < start or abs_e > end:
                continue
            if score_threshold is not None and score < score_threshold:
                continue
            out.append((abs_s - start, abs_e - start, tf_name, score, strand))
        except (IndexError, ValueError, TypeError) as exc:
            print(f"  hit parse fail: {h!r} ({exc})", flush=True)
    return out


def render(args):
    fp = ad.read_h5ad(args.fp_h5ad)
    obsm_keys = list(fp.obsm.keys())
    if not obsm_keys:
        raise ValueError(f"fp_h5ad has no obsm key: {args.fp_h5ad}")
    locus = obsm_keys[0]
    chrom, start, end = parse_locus(locus)
    width = end - start
    arr = fp.obsm[locus]
    scales = np.asarray(fp.uns["scales"])
    cond_names = fp.obs["name"].tolist()
    n_cond = arr.shape[0]
    print(f"[{args.gene}] {locus} arr={arr.shape} conditions={cond_names}",
          flush=True)

    scale_mask = scales <= args.max_scale
    scales_plot = scales[scale_mask]
    arr_plot = arr[:, scale_mask, :]
    print(f"[{args.gene}] plotting {len(scales_plot)} scales (<={args.max_scale} bp)",
          flush=True)

    tfs = tuple(t.strip() for t in args.tfs.split(",") if t.strip())
    if not tfs:
        raise ValueError("--tfs is empty after parsing; nothing to overlay")
    tf_color = parse_tf_colors(tfs, args.tf_colors)
    print(f"[{args.gene}] scanning {chrom}:{start}-{end} for {tfs}", flush=True)
    hits = scan_motifs_in_window(chrom, start, end, args.pfm, args.cache_dir,
                                 args.genome, tfs)
    print(f"[{args.gene}] raw n_hits="
          f"{len(hits) if hasattr(hits, '__len__') else '?'}", flush=True)
    motif_hits = parse_hits(hits, chrom, start, end, args.score_threshold)
    print(f"[{args.gene}] {len(motif_hits)} parsed motif hits in window", flush=True)
    for hs, he, tf, sc, st in motif_hits[:8]:
        print(f"  hit: rel {hs}-{he} ({tf} score={sc:.2f} strand={st:+d}) "
              f"abs {chrom}:{hs+start}-{he+start}", flush=True)

    vmin, vmax = 0, float(np.percentile(arr_plot, 99.5))
    fig, axes = plt.subplots(n_cond + 1, 1,
                             figsize=(12, 2.0 * n_cond + 1.6),
                             gridspec_kw={"height_ratios": [4] * n_cond + [1.2]},
                             sharex=True)
    axes = list(axes) if n_cond + 1 > 1 else [axes]
    for i in range(n_cond):
        ax = axes[i]
        im = ax.imshow(arr_plot[i], aspect="auto", origin="lower", cmap="Blues",
                       extent=[0, width, scales_plot[0], scales_plot[-1]],
                       vmin=vmin, vmax=vmax)
        ax.set_ylabel(f"{cond_names[i]}\nscale (bp)")
        ax.set_yticks([s for s in [2, 5, 10, 15, 20, 25, 30] if s <= scales_plot[-1]])
        plt.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="MSFP score")

    ax_track = axes[-1]
    ax_track.set_ylim(0, 1)
    ax_track.set_xlim(0, width)
    for hs, he, tf, sc, st in motif_hits:
        ax_track.add_patch(plt.Rectangle(
            (hs, 0.25), max(he - hs, 8), 0.5,
            color=tf_color.get(tf, "#444"), alpha=0.9, lw=0
        ))
    ax_track.set_yticks([])
    ax_track.set_xlabel(f"position (0-{width} bp) in {locus}")
    legend_label = " + ".join(f"{tf} ({tf_color.get(tf, '#444')})" for tf in tfs)
    ax_track.set_title(
        f"motif hits in window - {legend_label} - n={len(motif_hits)}",
        fontsize=9, loc="left",
    )

    title_ct = args.cell_type or fp.uns.get("cell_type", "")
    fig.suptitle(
        f"{args.gene} promoter MSFP ({title_ct}) - TF-binding scales <={args.max_scale} bp",
        fontsize=12,
    )
    fig.tight_layout()
    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"[{args.gene}] wrote {out_png}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--cell-type", default="")
    ap.add_argument("--fp-h5ad", required=True)
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--pfm", required=True)
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--genome", default="mm10")
    ap.add_argument("--max-scale", type=int, default=30)
    ap.add_argument("--tfs", required=True, help="comma-separated TF names")
    ap.add_argument("--tf-colors", default="",
                    help='"TF=#hex,TF=#hex,..." overrides; falls back to a cycling default palette')
    ap.add_argument("--score-threshold", type=float, default=None,
                    help="optional motif score floor; default = none "
                         "(matches MOTIF_SCAN_ENHANCERS default — scprinter applies its own internal cutoff)")
    args = ap.parse_args()
    render(args)


if __name__ == "__main__":
    sys.exit(main())
