#!/usr/bin/env python3
"""Of the enhancers that gained CCAN links to a gene under treatment, how many
host a TF motif of interest?

This is the cleanest spatial proxy for directionality available from
cross-sectional ATAC: for each gene whose CCAN expanded between two stratified
Cicero runs, partition the *gained* enhancer→gene links by whether the
enhancer region itself contains the TF motif site of interest. Enrichment
above the genome-wide background = direct TF binding likely drove the new
co-accessibility; at/below background = enhancer was activated indirectly.

Drop-in for any pipeline. Motif identity is set by --motif-name (axis label)
and --motif-bed (the actual coordinates). Defaults match SREBF1 in SDas but
the script is motif-agnostic — works for SOX10, AP-1, CTCF, anything you have
an enhancer-region BED for.

WARNING — coverage caveat: a motif BED derived from a single condition's
enhancer-peak universe will not contain motif hits at enhancers that exist
only in the other condition. To cover treatment-specific enhancers too, union
both per-condition peaksets, scan motifs on the union, then pass that BED.

Inputs:
  --links-control      <ccan_enhancer_gene_links_<control>.tsv>      (header: chr,start,end,CCAN,gene)
  --links-treatment    <ccan_enhancer_gene_links_<treatment>.tsv>
  --motif-bed          <bed of enhancer regions containing the motif>  (chr,start,end,...)
  --control-label      "5xFAD"
  --treatment-label    "SREBF1_OE"

Outputs:
  --out-tsv            per-gene table:
                       gene, n_<control>, n_<treatment>, delta,
                       n_gained_total, n_gained_motif, n_gained_no_motif,
                       frac_motif
  --out-png            stacked bar (top-N by delta + proportional twin panel)

Genes can be restricted with --gene-set (comma list) or --gene-file
(newline-separated). Otherwise the top-N by delta are plotted.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_links(tsv):
    """Return list of (chr, start, end, ccan, gene) rows."""
    rows = []
    with open(tsv) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, s, e, ccan, gene = parts[:5]
            rows.append((chrom, int(s), int(e), ccan, gene))
    return rows


def load_motif_bed(bed):
    """Return set of (chr, start, end) coordinates with motif."""
    motif_coords = set()
    with open(bed) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, s, e = parts[0], int(parts[1]), int(parts[2])
            motif_coords.add((chrom, s, e))
    return motif_coords


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--links-control", required=True)
    ap.add_argument("--links-treatment", required=True)
    ap.add_argument("--motif-bed", required=True)
    ap.add_argument("--control-label", default="5xFAD")
    ap.add_argument("--treatment-label", default="SREBF1_OE")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--top-n", type=int, default=25,
                    help="If no --gene-set/--gene-file given, plot top-N by delta")
    ap.add_argument("--gene-set", default=None,
                    help="Comma-separated explicit gene list to plot")
    ap.add_argument("--gene-file", default=None,
                    help="Newline-separated explicit gene list to plot")
    ap.add_argument("--min-delta", type=int, default=1,
                    help="Skip genes whose delta < this in the plot")
    ap.add_argument("--motif-name", default="SREBF1",
                    help="Used in axis labels / title only")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    fad_rows = load_links(args.links_control)
    soe_rows = load_links(args.links_treatment)
    motif_coords = load_motif_bed(args.motif_bed)

    fad_keys = {(r[0], r[1], r[2], r[4]) for r in fad_rows}
    soe_keys = {(r[0], r[1], r[2], r[4]) for r in soe_rows}

    by_gene_fad = defaultdict(int)
    by_gene_soe = defaultdict(int)
    for r in fad_rows: by_gene_fad[r[4]] += 1
    for r in soe_rows: by_gene_soe[r[4]] += 1

    # Per-gene gained = OE-only (chr,start,end,gene)
    gained_by_gene = defaultdict(list)
    for k in soe_keys - fad_keys:
        chrom, s, e, gene = k
        gained_by_gene[gene].append((chrom, s, e))

    # Tabulate
    rows = []
    universe_genes = set(by_gene_fad) | set(by_gene_soe)
    for g in universe_genes:
        a, b = by_gene_fad.get(g, 0), by_gene_soe.get(g, 0)
        delta = b - a
        gained = gained_by_gene.get(g, [])
        n_gained_total = len(gained)
        n_gained_motif = sum(1 for c in gained if c in motif_coords)
        n_gained_no_motif = n_gained_total - n_gained_motif
        frac = n_gained_motif / n_gained_total if n_gained_total else 0.0
        rows.append((g, a, b, delta, n_gained_total, n_gained_motif,
                     n_gained_no_motif, frac))

    rows.sort(key=lambda r: (-r[3], -r[5], r[0]))

    out_tsv = Path(args.out_tsv); out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as f:
        f.write(
            f"gene\tn_{args.control_label}\tn_{args.treatment_label}\tdelta\t"
            f"n_gained_total\tn_gained_with_{args.motif_name}\t"
            f"n_gained_without_{args.motif_name}\tfrac_with_{args.motif_name}\n"
        )
        for r in rows:
            f.write("\t".join(str(x) for x in r[:7]) + f"\t{r[7]:.4f}\n")
    print(f"wrote {out_tsv}")

    # Pick which genes to plot
    if args.gene_set:
        plot_genes = [g.strip() for g in args.gene_set.split(",") if g.strip()]
    elif args.gene_file:
        with open(args.gene_file) as f:
            plot_genes = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    else:
        plot_genes = [r[0] for r in rows
                      if r[3] >= args.min_delta and r[4] > 0][:args.top_n]

    g2row = {r[0]: r for r in rows}
    plot_rows = [g2row[g] for g in plot_genes if g in g2row]
    plot_rows = [r for r in plot_rows if r[4] > 0]
    if not plot_rows:
        print("No plottable genes (no gained enhancers in selected set).")
        return

    genes = [r[0] for r in plot_rows]
    n_motif = np.array([r[5] for r in plot_rows])
    n_no_motif = np.array([r[6] for r in plot_rows])
    delta = np.array([r[3] for r in plot_rows])
    frac = np.array([r[7] for r in plot_rows])

    # Background expectation: total motif-enhancer count / total enhancer count
    # — proxy for "fraction of any enhancer that contains the motif"
    total_unique_enhancers = len({(r[0], r[1], r[2]) for r in fad_rows} |
                                 {(r[0], r[1], r[2]) for r in soe_rows})
    n_motif_in_universe = sum(
        1 for c in {(r[0], r[1], r[2]) for r in fad_rows} | {(r[0], r[1], r[2]) for r in soe_rows}
        if c in motif_coords
    )
    bg_frac = n_motif_in_universe / total_unique_enhancers if total_unique_enhancers else 0
    print(f"  background fraction of enhancers with {args.motif_name} motif: "
          f"{n_motif_in_universe}/{total_unique_enhancers} = {bg_frac:.3%}")

    fig, axes = plt.subplots(1, 2, figsize=(15, max(5, 0.36 * len(genes) + 1.4)),
                             gridspec_kw={"width_ratios": [1.2, 1.0]})
    ax_abs, ax_prop = axes
    y = np.arange(len(genes))

    # Absolute stacked bar: gained-enhancers split by motif
    ax_abs.barh(y, n_motif, color="#dd8452", edgecolor="black", linewidth=0.4,
                label=f"gained, {args.motif_name}+")
    ax_abs.barh(y, n_no_motif, left=n_motif, color="#cccccc",
                edgecolor="black", linewidth=0.4,
                label=f"gained, {args.motif_name}−")
    for yi, m, nm, d in zip(y, n_motif, n_no_motif, delta):
        total = m + nm
        if total:
            ax_abs.text(total + 0.4, yi, f"Δ={d:+d}, mot={m}/{total}",
                        va="center", fontsize=8)
    ax_abs.set_yticks(y)
    ax_abs.set_yticklabels(genes, fontsize=8)
    ax_abs.invert_yaxis()
    ax_abs.set_xlabel(f"Gained enhancers ({args.treatment_label}-only links)")
    ax_abs.set_title(
        f"Gained-CCAN enhancer counts split by {args.motif_name} motif presence",
        fontsize=10
    )
    ax_abs.legend(loc="lower right", frameon=False, fontsize=8)
    ax_abs.spines["top"].set_visible(False)
    ax_abs.spines["right"].set_visible(False)

    # Proportional panel: % of gained that have motif
    ax_prop.barh(y, frac * 100, color="#dd8452", edgecolor="black",
                 linewidth=0.4, label=f"% gained with {args.motif_name}")
    ax_prop.axvline(bg_frac * 100, color="black", linestyle="--", linewidth=0.7,
                    label=f"background {bg_frac:.1%}")
    for yi, fr, m, total in zip(y, frac, n_motif, n_motif + n_no_motif):
        x_text = fr * 100 + 1.5
        ax_prop.text(x_text, yi, f"{m}/{total} ({fr:.0%})",
                     va="center", fontsize=8)
    ax_prop.set_yticks(y)
    ax_prop.set_yticklabels(genes, fontsize=8)
    ax_prop.invert_yaxis()
    ax_prop.set_xlim(0, max(110, frac.max() * 100 + 30))
    ax_prop.set_xlabel(f"% of gained enhancers containing {args.motif_name} motif")
    ax_prop.set_title(
        f"Per-gene enrichment for {args.motif_name} in newly co-accessible enhancers",
        fontsize=10
    )
    ax_prop.legend(loc="lower right", frameon=False, fontsize=8)
    ax_prop.spines["top"].set_visible(False)
    ax_prop.spines["right"].set_visible(False)

    title = (
        args.title
        or f"Directionality proxy — {args.motif_name} motif presence in CCAN-gained enhancers "
           f"({args.treatment_label} − {args.control_label})"
    )
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    out_png = Path(args.out_png); out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
