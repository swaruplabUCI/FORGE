#!/usr/bin/env python3
"""
render_cicero_lollipop.py — FORGE pipeline module

Lollipop chart: TF-motif enhancer → promoter CCAN arcs per gene,
comparing ctrl vs trt conditions for a single cell type.

Differences from the oneOff 06_cicero_lollipop.py:
  - No hardcoded path constants (RESULTS, OUT_DIR, GTF_PATH removed).
  - --gtf is REQUIRED (no default).
  - --outdir defaults to "." (callee's work directory in Nextflow).
  - PDF + PNG always written (unchanged from oneOff).

Algorithm
---------
For each condition (ctrl / trt):
  1. Load cicero_connections.tsv.gz from {ccan-base}/{condition}/.
  2. Filter to coaccess >= min-coacc.
  3. Keep arcs where one end overlaps a TF-motif peak AND the other end
     overlaps a gene promoter (TSS ± upstream/downstream bp).
  4. Annotate the TF-motif peak as Distal / Intronic / Exonic / Promoter.
  5. Count arcs per (gene, condition) → Δ = trt - ctrl.

Plot
----
  y-axis: top_n genes ranked by |Δ|
  x-axis: Δ CCAN arc count
  stem color: red = trt gain (Δ > 0), blue = ctrl gain (Δ < 0)
  dot color : green = Distal motif, blue = Intronic motif
  dot size  : scaled to max coaccessibility score

Usage:
    singularity exec scgpu_extended.sif python3 render_cicero_lollipop.py \\
        --motif-peaks   motif_peaks_Spi1.bed \\
        --ccan-base     results/cicero/per_ct/Microglia_NN \\
        --cell-type     Microglia_NN \\
        --gtf           /ref/Gencode_GRCm38/gencode.vM10.annotation.gtf \\
        --outdir        .
"""

import argparse
import bisect
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ── Color palettes ─────────────────────────────────────────────────────────────

COND_COLORS = {"TG": "#e41a1c", "WT": "#377eb8"}
ANN_COLORS  = {"Distal": "#2ca02c", "Intronic": "#1f77b4"}


# ── GTF helpers ────────────────────────────────────────────────────────────────

def parse_peak(p):
    m = re.match(r"(chr[^:]+):(\d+)-(\d+)", str(p))
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else (None, None, None)


def load_gtf(gtf_path, upstream=2000, downstream=500):
    """Parse GTF → (promoters_arr, exons_arr, gene_bodies_arr, prom_gene_map)."""
    promoters     = defaultdict(list)
    exons         = defaultdict(list)
    gene_bodies   = defaultdict(list)
    prom_gene_map = defaultdict(list)

    _open = (lambda p: gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p))
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            chrom, feat, strand = parts[0], parts[2], parts[6]
            s, e = int(parts[3]) - 1, int(parts[4])
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', parts[8]))
            gene  = attrs.get("gene_name", "")
            if feat == "gene":
                gene_bodies[chrom].append((s, e))
                tss = s if strand == "+" else e
                ps  = max(0, tss - upstream)   if strand == "+" else max(0, tss - downstream)
                pe  = tss + downstream          if strand == "+" else tss + upstream
                promoters[chrom].append((ps, pe))
                if gene:
                    prom_gene_map[chrom].append((ps, pe, gene))
            elif feat == "exon":
                exons[chrom].append((s, e))

    def _arr(d):
        return {c: np.array(sorted(v), dtype=np.int64) for c, v in d.items() if v}

    return _arr(promoters), _arr(exons), _arr(gene_bodies), prom_gene_map


def _overlaps(arr, s, e):
    if arr is None or len(arr) == 0:
        return False
    idx = bisect.bisect_left(arr[:, 0].tolist(), e)
    for i in range(idx - 1, -1, -1):
        if arr[i, 1] > s:
            return True
        if arr[i, 0] < s - 2_000_000:
            break
    return False


def annotate(chrom, s, e, promoters_arr, exons_arr, gene_bodies_arr):
    if _overlaps(promoters_arr.get(chrom), s, e):
        return "Promoter"
    if _overlaps(exons_arr.get(chrom), s, e):
        return "Exonic"
    if _overlaps(gene_bodies_arr.get(chrom), s, e):
        return "Intronic"
    return "Distal"


def peaks_to_genes(chrom, s, e, prom_gene_map):
    return {g for ps, pe, g in prom_gene_map.get(chrom, []) if ps < e and pe > s}


# ── Core analysis ──────────────────────────────────────────────────────────────

def load_motif_peaks(bed_path):
    """Load TF-motif BED → set of 'chr:start-end' peak strings."""
    peaks = set()
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split('\t')
            if len(cols) < 3:
                continue
            peaks.add(f"{cols[0].strip()}:{cols[1].strip()}-{cols[2].strip()}")
    print(f"[lollipop] {len(peaks)} motif peaks loaded from {bed_path}", flush=True)
    return peaks


def build_gene_table(motif_peaks, promoters_arr, exons_arr, gene_bodies_arr,
                     prom_gene_map, ccan_base, ctrl_cond, trt_cond, min_coacc=0.15):
    """Build per-gene arc count table comparing ctrl vs trt."""
    rows = []
    for cond in (ctrl_cond, trt_cond):
        gz = Path(ccan_base) / cond / "cicero_connections.tsv.gz"
        if not gz.exists():
            print(f"  [WARN] missing {gz}, skipping", flush=True)
            continue
        with gzip.open(gz, "rt") as f:
            df = pd.read_csv(f, sep="\t")
        df = df[df["coaccess"] >= min_coacc]
        mask_p1 = df["Peak1"].isin(motif_peaks)
        mask_p2 = df["Peak2"].isin(motif_peaks)

        # Arc: motif peak is Peak1, partner must reach a promoter
        for _, row in df[mask_p1].iterrows():
            c, s, e = parse_peak(row["Peak1"])
            ann1 = annotate(c, s, e, promoters_arr, exons_arr, gene_bodies_arr) if c else "Unknown"
            if ann1 not in ("Distal", "Intronic"):
                continue
            c2, s2, e2 = parse_peak(row["Peak2"])
            if c2 is None:
                continue
            if annotate(c2, s2, e2, promoters_arr, exons_arr, gene_bodies_arr) != "Promoter":
                continue
            for g in peaks_to_genes(c2, s2, e2, prom_gene_map):
                rows.append({"cond": cond, "motif_peak": row["Peak1"],
                             "motif_ann": ann1, "coaccess": row["coaccess"], "gene": g})

        # Arc: motif peak is Peak2, partner must reach a promoter
        for _, row in df[mask_p2 & ~mask_p1].iterrows():
            c, s, e = parse_peak(row["Peak2"])
            ann2 = annotate(c, s, e, promoters_arr, exons_arr, gene_bodies_arr) if c else "Unknown"
            if ann2 not in ("Distal", "Intronic"):
                continue
            c1, s1, e1 = parse_peak(row["Peak1"])
            if c1 is None:
                continue
            if annotate(c1, s1, e1, promoters_arr, exons_arr, gene_bodies_arr) != "Promoter":
                continue
            for g in peaks_to_genes(c1, s1, e1, prom_gene_map):
                rows.append({"cond": cond, "motif_peak": row["Peak2"],
                             "motif_ann": ann2, "coaccess": row["coaccess"], "gene": g})

    arc_df = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["cond", "motif_peak", "motif_ann", "coaccess", "gene"])

    gene_rows = []
    for gene, sub in arc_df.groupby("gene"):
        ctrl = sub[sub["cond"] == ctrl_cond]
        trt  = sub[sub["cond"] == trt_cond]
        ann_counts   = sub["motif_ann"].value_counts()
        dominant_ann = ann_counts.index[0] if len(ann_counts) else "Distal"
        gene_rows.append({
            "gene":           gene,
            "n_ctrl":         len(ctrl),
            "n_trt":          len(trt),
            "delta":          len(trt) - len(ctrl),
            "max_coacc_trt":  float(trt["coaccess"].max())  if len(trt)  else 0.0,
            "max_coacc_ctrl": float(ctrl["coaccess"].max()) if len(ctrl) else 0.0,
            "dominant_ann":   dominant_ann,
        })
    gene_df = pd.DataFrame(gene_rows).sort_values("delta", ascending=False)
    return gene_df


# ── Plotting ───────────────────────────────────────────────────────────────────

def plot_lollipop(gene_df, out_base, cell_type, ctrl_cond, trt_cond, top_n=40):
    """Render lollipop chart; writes {out_base}.pdf and {out_base}.png."""
    plot_df = (gene_df.assign(abs_delta=gene_df["delta"].abs())
               .sort_values(["abs_delta", "n_trt"], ascending=[False, False])
               .head(top_n)
               .sort_values("delta"))

    n     = len(plot_df)
    fig_h = max(4, n * 0.28 + 1.8)
    fig, ax = plt.subplots(figsize=(5.5, fig_h))
    y = np.arange(n)

    # Condition-exclusive row highlights
    for i, (_, row) in enumerate(plot_df.iterrows()):
        if row["n_ctrl"] == 0 and row["n_trt"] > 0:
            ax.axhspan(i - 0.45, i + 0.45, color="#fff0f0", zorder=0)
        elif row["n_trt"] == 0 and row["n_ctrl"] > 0:
            ax.axhspan(i - 0.45, i + 0.45, color="#f0f4ff", zorder=0)

    ax.axvline(0, color="0.6", linewidth=0.8, zorder=1)

    # Stems
    for i, (_, row) in enumerate(plot_df.iterrows()):
        d     = row["delta"]
        color = COND_COLORS.get(trt_cond, "#e41a1c") if d >= 0 \
                else COND_COLORS.get(ctrl_cond, "#377eb8")
        ax.plot([0, d], [i, i], color=color, linewidth=1.2,
                zorder=2, solid_capstyle="round")

    # Dots (size ∝ max coaccess, color = annotation class)
    for i, (_, row) in enumerate(plot_df.iterrows()):
        d    = row["delta"]
        mc   = row["max_coacc_trt"] if d >= 0 else row["max_coacc_ctrl"]
        size = 20 + 200 * mc
        color = ANN_COLORS.get(row["dominant_ann"], "#7f7f7f")
        ax.scatter(d, i, s=size, color=color, zorder=3,
                   edgecolors="white", linewidths=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"].values, fontsize=7)
    ax.set_xlabel(f"Δ CCAN arcs ({trt_cond} − {ctrl_cond})", fontsize=8)
    ax.set_title(
        f"TF-motif enhancer → promoter CCAN arcs\n{cell_type}  ·  top {n} genes by |Δ|",
        fontsize=9,
    )

    xmax = max(abs(plot_df["delta"].max()), abs(plot_df["delta"].min()), 1)
    ax.set_xlim(-xmax - 1.5, xmax + 1.5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_axisbelow(True)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color="0.88", linewidth=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    # Legends
    ann_patches = [mpatches.Patch(color=ANN_COLORS[a], label=f"{a} motif")
                   for a in ("Distal", "Intronic")]
    cond_patches = [
        mpatches.Patch(color=COND_COLORS.get(trt_cond,  "#e41a1c"),
                       label=f"{trt_cond} gain (Δ > 0)"),
        mpatches.Patch(color=COND_COLORS.get(ctrl_cond, "#377eb8"),
                       label=f"{ctrl_cond} gain (Δ < 0)"),
    ]
    size_dots = [
        plt.scatter([], [], s=20 + 200 * v, color="0.5",
                    edgecolors="white", linewidths=0.5,
                    label=f"coacc = {v:.2f}")
        for v in (0.15, 0.25, 0.40)
    ]
    leg1 = ax.legend(handles=ann_patches + cond_patches,
                     loc="lower right", fontsize=6.5, framealpha=0.85,
                     title="Motif context / direction", title_fontsize=6.5)
    ax.legend(handles=size_dots, loc="upper left", fontsize=6,
              title="Max coaccessibility", title_fontsize=6,
              framealpha=0.85, labelspacing=0.7)
    ax.add_artist(leg1)
    ax.text(1.01, 0.98,
            f"Pink bg = {trt_cond}-exclusive\nBlue bg = {ctrl_cond}-exclusive",
            transform=ax.transAxes, fontsize=5.5, va="top", ha="left", color="0.4")

    plt.tight_layout()
    for ext in ("pdf", "png"):
        p = Path(f"{out_base}.{ext}")
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=200, bbox_inches="tight")
        print(f"[lollipop] wrote {p}", flush=True)
    plt.close(fig)


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    # ── Required ──────────────────────────────────────────────────────────────
    ap.add_argument("--motif-peaks",     required=True,
                    help="BED file (chr, start, end) of TF-motif-containing peaks")
    ap.add_argument("--ccan-base",       required=True,
                    help="Base dir with {ctrl_cond}/ and {trt_cond}/ sub-dirs "
                         "each containing cicero_connections.tsv.gz")
    ap.add_argument("--cell-type",       required=True,
                    help="Cell type label for figure title")
    ap.add_argument("--gtf",             required=True,
                    help="Gencode/Ensembl GTF (gzip or plain; no default in pipeline mode)")
    # ── Optional ──────────────────────────────────────────────────────────────
    ap.add_argument("--ctrl-condition",  default="WT")
    ap.add_argument("--trt-condition",   default="TG")
    ap.add_argument("--min-coacc",       type=float, default=0.15,
                    help="Minimum coaccessibility score for arc inclusion (default 0.15)")
    ap.add_argument("--top-n",           type=int,   default=40,
                    help="Number of top genes to show by |Δ| (default 40)")
    ap.add_argument("--out-tag",         default="cicero_lollipop",
                    help="Output filename stem (default 'cicero_lollipop')")
    ap.add_argument("--outdir",          default=".",
                    help="Output directory (default: current working dir)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[lollipop] cell_type={args.cell_type}  "
          f"ctrl={args.ctrl_condition}  trt={args.trt_condition}  "
          f"min_coacc={args.min_coacc}  top_n={args.top_n}", flush=True)

    motif_peaks = load_motif_peaks(args.motif_peaks)
    if not motif_peaks:
        print("[lollipop][WARN] no motif peaks loaded — empty output.", flush=True)
        return

    print("[lollipop] parsing GTF (~30 s) ...", flush=True)
    promoters_arr, exons_arr, gene_bodies_arr, prom_gene_map = load_gtf(args.gtf)

    print("[lollipop] building gene arc table ...", flush=True)
    gene_df = build_gene_table(
        motif_peaks, promoters_arr, exons_arr, gene_bodies_arr,
        prom_gene_map, args.ccan_base,
        ctrl_cond=args.ctrl_condition, trt_cond=args.trt_condition,
        min_coacc=args.min_coacc,
    )
    print(f"[lollipop] {len(gene_df)} genes total | "
          f"{args.trt_condition}-exclusive: {(gene_df.n_ctrl == 0).sum()}  "
          f"{args.ctrl_condition}-exclusive: {(gene_df.n_trt == 0).sum()}  "
          f"both: {((gene_df.n_ctrl > 0) & (gene_df.n_trt > 0)).sum()}", flush=True)

    tsv_path = outdir / f"{args.out_tag}_genes.tsv"
    gene_df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[lollipop] gene table → {tsv_path}", flush=True)

    top_n = min(args.top_n, len(gene_df))
    if top_n == 0:
        print("[lollipop][WARN] no genes to plot — check motif peaks / ccan-base.",
              flush=True)
        return

    plot_lollipop(
        gene_df, outdir / args.out_tag,
        cell_type=args.cell_type,
        ctrl_cond=args.ctrl_condition,
        trt_cond=args.trt_condition,
        top_n=top_n,
    )
    print("[lollipop] done.", flush=True)


if __name__ == "__main__":
    main()
