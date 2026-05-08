#!/usr/bin/env python3
"""Per-condition cis-rewiring bar plot with the OE bar stacked to show how
many of the *gained* CCAN-enhancers carry a TF motif.

Rank order is the original delta-link ranking — the user's point being that
ordering by motif matching is the wrong primary axis. Instead, the original
bar plot ordering is preserved and the OE bar is split into three layers:
  - shared with control (matching the control bar's height)
  - gained, motif−       (light grey)
  - gained, motif+       (orange)

So you read each pair as: blue = control link count; stacked OE = control +
gain_no_motif + gain_with_motif. Most genes will show very little orange —
that is the point; this preserves the visibility of the "newly accessible
enhancers don't host the TF directly" finding without sorting by it.

Usage:
  render_cis_rewiring_motif_stack.py \\
    --links-control gene_links_5xFAD.tsv \\
    --links-treatment gene_links_SREBF1_OE.tsv \\
    --motif-bed srebf1_enhancer_regions.bed \\
    --motif-name SREBF1 \\
    --panel "Myelin priority:Myrf,Cers2,Cldn11,Gjc2" \\
    --panel "Lipid pathway:Hmgcr,Hmgcs1,Insig1,Insig2,Ldlr,Scap,Sqle,Srebf1,Srebf2" \\
    --out-png out.png --out-tsv out.tsv

Drop-in for any pipeline. Pass any motif BED / name combination.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PSEUDOGENE_PATTERNS = [
    r"^Gm\d+$", r"Rik$", r"^AU\d+$", r"^Mir\d+",
    r"^[0-9]+[A-Z][a-z0-9]+", r"^LOC\d+$",
]
PSEUDOGENE_RE = re.compile("|".join(PSEUDOGENE_PATTERNS))


def is_biology_informative(gene):
    return PSEUDOGENE_RE.search(gene) is None


def load_links_keys(tsv):
    """Return list of (chr, start, end, ccan, gene) and per-gene set of (chr,start,end)."""
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
    coords = set()
    with open(bed) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                coords.add((parts[0], int(parts[1]), int(parts[2])))
    return coords


def parse_panels(panel_strs, panel_files):
    panels = []
    for s in panel_strs or []:
        if ":" not in s:
            raise SystemExit(f"--panel must be 'label:gene1,gene2,...' (got: {s})")
        label, glist = s.split(":", 1)
        genes = [g.strip() for g in glist.split(",") if g.strip()]
        if genes:
            panels.append((label.strip(), genes))
    for path in panel_files or []:
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) == 2:
                    label, glist = parts
                    genes = [g.strip() for g in glist.split(",") if g.strip()]
                    if genes:
                        panels.append((label.strip(), genes))
    return panels


def top_n_panel_by_delta(fad_rows, soe_rows, top_n, min_delta, exclude_pseudogenes,
                        label="Top-N by Δ enhancer connections"):
    fad_per = defaultdict(set)
    soe_per = defaultdict(set)
    for r in fad_rows:
        fad_per[r[4]].add((r[0], r[1], r[2]))
    for r in soe_rows:
        soe_per[r[4]].add((r[0], r[1], r[2]))
    all_genes = set(fad_per) | set(soe_per)
    deltas = []
    for g in all_genes:
        if exclude_pseudogenes and not is_biology_informative(g):
            continue
        d = len(soe_per.get(g, set())) - len(fad_per.get(g, set()))
        if d >= min_delta:
            deltas.append((g, d))
    deltas.sort(key=lambda x: x[1], reverse=True)
    genes = [g for g, _ in deltas[:top_n]]
    if not genes:
        raise SystemExit(f"No genes survived top-N selection (min-delta={min_delta}, exclude_pseudogenes={exclude_pseudogenes})")
    return [(label, genes)]


def panel_stack_horizontal(ax, genes, fad_rows, soe_rows, motif_coords,
                           control_label, treatment_label, motif_name, title):
    """Horizontal stacked bars — used for top-N agnostic mode where N is large."""
    fad_keys_per = defaultdict(set)
    soe_keys_per = defaultdict(set)
    for r in fad_rows:
        fad_keys_per[r[4]].add((r[0], r[1], r[2]))
    for r in soe_rows:
        soe_keys_per[r[4]].add((r[0], r[1], r[2]))

    n_fad, n_shared, n_gnm, n_gm, rows_out = [], [], [], [], []
    for g in genes:
        fk = fad_keys_per.get(g, set()); sk = soe_keys_per.get(g, set())
        shared = fk & sk; gained = sk - fk
        gm = {c for c in gained if c in motif_coords}; gnm = gained - gm
        n_fad.append(len(fk)); n_shared.append(len(shared))
        n_gnm.append(len(gnm)); n_gm.append(len(gm))
        rows_out.append((g, len(fk), len(sk), len(sk) - len(fk),
                         len(shared), len(gnm), len(gm)))
    n_fad = np.array(n_fad); n_shared = np.array(n_shared)
    n_gnm = np.array(n_gnm); n_gm = np.array(n_gm)

    y = np.arange(len(genes)); h = 0.36
    ax.barh(y - h/2, n_fad, h, color="#4c72b0", edgecolor="black", linewidth=0.3,
            label=control_label)
    ax.barh(y + h/2, n_shared, h, color="#4c72b0", alpha=0.55,
            edgecolor="black", linewidth=0.3,
            label=f"{treatment_label}: shared with {control_label}")
    ax.barh(y + h/2, n_gnm, h, left=n_shared, color="#cccccc",
            edgecolor="black", linewidth=0.3,
            label=f"{treatment_label}: gained, {motif_name}−")
    ax.barh(y + h/2, n_gm, h, left=n_shared + n_gnm, color="#dd8452",
            edgecolor="black", linewidth=0.3,
            label=f"{treatment_label}: gained, {motif_name}+")

    xmax = max(int((n_shared + n_gnm + n_gm).max()), int(n_fad.max()), 1)
    for yi, h_total, m in zip(y + h/2, n_shared + n_gnm + n_gm, n_gm):
        ax.text(h_total + xmax * 0.012, yi,
                f"{int(h_total)}{'  ('+str(int(m))+ 'm)' if m else ''}",
                va="center", fontsize=7)
    for yi, v in zip(y - h/2, n_fad):
        if v: ax.text(v + xmax * 0.012, yi, str(int(v)), va="center", fontsize=7)

    ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Enhancer connections (CCAN)")
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0, xmax * 1.2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    return rows_out


def panel_stack(ax, genes, fad_rows, soe_rows, motif_coords,
                control_label, treatment_label, motif_name, title):
    fad_keys_per = defaultdict(set)
    soe_keys_per = defaultdict(set)
    for r in fad_rows:
        fad_keys_per[r[4]].add((r[0], r[1], r[2]))
    for r in soe_rows:
        soe_keys_per[r[4]].add((r[0], r[1], r[2]))

    n_fad = []
    n_shared = []
    n_gained_no_motif = []
    n_gained_motif = []
    rows_out = []
    for g in genes:
        fk = fad_keys_per.get(g, set())
        sk = soe_keys_per.get(g, set())
        shared = fk & sk
        gained = sk - fk
        gained_motif = {c for c in gained if c in motif_coords}
        gained_no_motif = gained - gained_motif

        n_fad.append(len(fk))
        n_shared.append(len(shared))
        n_gained_no_motif.append(len(gained_no_motif))
        n_gained_motif.append(len(gained_motif))
        rows_out.append((g, len(fk), len(sk), len(sk) - len(fk),
                         len(shared), len(gained_no_motif), len(gained_motif)))

    n_fad = np.array(n_fad)
    n_shared = np.array(n_shared)
    n_gained_no_motif = np.array(n_gained_no_motif)
    n_gained_motif = np.array(n_gained_motif)

    x = np.arange(len(genes))
    w = 0.36

    # control bar (solid blue)
    ax.bar(x - w/2, n_fad, w, color="#4c72b0", edgecolor="black", linewidth=0.4,
           label=control_label)

    # treatment bar — stacked: shared (matched-blue), gained no motif (grey), gained motif (orange)
    ax.bar(x + w/2, n_shared, w, color="#4c72b0", edgecolor="black", linewidth=0.4,
           label=f"{treatment_label}: shared with {control_label}", alpha=0.55)
    ax.bar(x + w/2, n_gained_no_motif, w, bottom=n_shared,
           color="#cccccc", edgecolor="black", linewidth=0.4,
           label=f"{treatment_label}: gained, {motif_name}−")
    ax.bar(x + w/2, n_gained_motif, w, bottom=n_shared + n_gained_no_motif,
           color="#dd8452", edgecolor="black", linewidth=0.4,
           label=f"{treatment_label}: gained, {motif_name}+")

    ymax = max(int((n_shared + n_gained_no_motif + n_gained_motif).max()),
               int(n_fad.max()), 1) + 3
    for xi, v in zip(x - w/2, n_fad):
        if v: ax.text(xi, v + ymax * 0.015, str(int(v)), ha="center", fontsize=8)
    for xi, h_total, m in zip(x + w/2,
                               n_shared + n_gained_no_motif + n_gained_motif,
                               n_gained_motif):
        ax.text(xi, h_total + ymax * 0.015,
                f"{int(h_total)}{'  ('+str(int(m))+ 'm)' if m else ''}",
                ha="center", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Enhancer connections (CCAN)")
    ax.set_title(title, fontsize=10)
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return rows_out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--links-control", required=True)
    ap.add_argument("--links-treatment", required=True)
    ap.add_argument("--motif-bed", required=True)
    ap.add_argument("--motif-name", default="SREBF1")
    ap.add_argument("--control-label", default="5xFAD")
    ap.add_argument("--treatment-label", default="SREBF1_OE")
    ap.add_argument("--panel", action="append")
    ap.add_argument("--panel-file", action="append")
    ap.add_argument("--top-n-by-delta", type=int, default=None,
                    help="Agnostic mode: ignore --panel and rank top-N genes by Δ (n_treatment − n_control).")
    ap.add_argument("--min-delta", type=int, default=1,
                    help="Top-N mode: minimum Δ to include a gene (default 1).")
    ap.add_argument("--exclude-pseudogenes", action="store_true",
                    help="Top-N mode: drop Gm*/*Rik/Mir*/AU*/LOC*/numeric-prefixed predicted genes.")
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--title", default="Cis-rewiring with motif overlay in gained enhancers")
    ap.add_argument("--width-per-panel", type=float, default=5.4)
    ap.add_argument("--height", type=float, default=5.0)
    args = ap.parse_args()

    fad = load_links_keys(args.links_control)
    soe = load_links_keys(args.links_treatment)
    motif_coords = load_motif_bed(args.motif_bed)

    rows = []

    if args.top_n_by_delta is not None:
        panels = top_n_panel_by_delta(
            fad, soe, args.top_n_by_delta, args.min_delta, args.exclude_pseudogenes,
            label=f"Top {args.top_n_by_delta} genes by Δ enhancer connections "
                  f"({args.treatment_label} − {args.control_label})",
        )
        label, genes = panels[0]
        height = max(4.0, 0.32 * len(genes) + 1.6)
        fig, ax = plt.subplots(1, 1, figsize=(10.0, height))
        out_rows = panel_stack_horizontal(
            ax, genes, fad, soe, motif_coords,
            args.control_label, args.treatment_label, args.motif_name, label,
        )
        for r in out_rows:
            rows.append((label,) + r)
        handles, labels = ax.get_legend_handles_labels()
        seen, uniq_h, uniq_l = [], [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.append(l); uniq_h.append(h); uniq_l.append(l)
        ax.legend(uniq_h, uniq_l, loc="lower right", frameon=False, fontsize=8)
    else:
        panels = parse_panels(args.panel, args.panel_file)
        if not panels:
            raise SystemExit("Need at least one --panel/--panel-file or --top-n-by-delta")
        n_panels = len(panels)
        fig, axes = plt.subplots(1, n_panels,
                                 figsize=(args.width_per_panel * n_panels, args.height),
                                 squeeze=False)
        axes = axes.ravel()
        for ax, (label, genes) in zip(axes, panels):
            out_rows = panel_stack(ax, genes, fad, soe, motif_coords,
                                   args.control_label, args.treatment_label,
                                   args.motif_name, label)
            for r in out_rows:
                rows.append((label,) + r)
        handles, labels = axes[-1].get_legend_handles_labels()
        seen, uniq_h, uniq_l = [], [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.append(l); uniq_h.append(h); uniq_l.append(l)
        axes[-1].legend(uniq_h, uniq_l, loc="upper right", frameon=False, fontsize=8)

    fig.suptitle(args.title, fontsize=12)
    fig.tight_layout()
    out_png = Path(args.out_png); out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"wrote {out_png}")

    out_tsv = Path(args.out_tsv); out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as f:
        f.write(
            f"panel\tgene\tn_{args.control_label}\tn_{args.treatment_label}\tdelta\t"
            f"n_shared\tn_gained_without_{args.motif_name}\tn_gained_with_{args.motif_name}\n"
        )
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {out_tsv}")


if __name__ == "__main__":
    main()
