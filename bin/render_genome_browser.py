#!/usr/bin/env python3
"""
render_genome_browser.py — FORGE pipeline module

Multi-track ATAC genome browser rendered with matplotlib + pyBigWig.
No pyGenomeTracks dependency, no CCAN arcs.

Modes
-----
absolute      : one fill track per cell type, colored by cell class.
differential  : one overlay track per cell type showing:
                  grey   — shared accessibility (min of ctrl and trt)
                  red    — TG-unique (trt > ctrl)
                  blue   — WT-unique (ctrl > trt)
                Requires condition-split BigWigs via --bw-manifest
                (produced by export_atac_bigwigs.py --condition-col).

Inputs
------
--bw-manifest   JSON manifest from export_atac_bigwigs.py (optional in absolute mode)
                  .bigwigs        { ct: filename }        ← absolute mode
                  .by_condition   { ct: { cond: filename } } ← differential mode
--bw-dir        BigWig base directory.
                  With --bw-manifest: manifest filenames are resolved relative here.
                  Without --bw-manifest (absolute mode only): scanned for {ct}.bw files.
--gene          Gene name
--gtf           Gencode/Ensembl GTF (gzip or plain)
--cell-types    Comma-separated ordered list to display
--mode          absolute | differential (default: absolute)
--ctrl-condition / --trt-condition   condition labels in manifest (default WT / TG)

Usage (pipeline, differential):
    singularity exec snapatac_extended.sif python3 render_genome_browser.py \\
        --bw-manifest  bigwigs/manifest.json \\
        --bw-dir       bigwigs/ \\
        --gene         Trem2 \\
        --gtf          /ref/gencode.vM10.annotation.gtf \\
        --cell-types   Micro-PVM,Astro,Oligo \\
        --mode         differential \\
        --out-png      browser_Trem2.png

Usage (oneOff / absolute, no manifest):
    singularity exec snapatac_extended.sif python3 render_genome_browser.py \\
        --bw-dir   results/pycistopic/bigwigs/.../pseudobulk_bw_files \\
        --gene     Trem2 \\
        --gtf      /ref/gencode.vM10.annotation.gtf \\
        --cell-types Micro-PVM,Astro,Oligo \\
        --out-png  browser_Trem2.png
"""

import argparse
import gzip
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec


# ── Color palettes ─────────────────────────────────────────────────────────────

# Colors sourced from SDas plot_srebf1_locus.py; space-separated names match
# the cell_type_prediction labels produced by CellTypist / FORGE annotation.
CT_COLORS = {
    # Glia (fineglia)
    'Oligo NN':              '#8c564b',
    'OPC NN':                '#e377c2',
    'Astro-OLF NN':          '#1F77B4',
    'Astro-NT NN':           '#AEC7E8',
    'Astro-TE NN':           '#6BAED6',
    'Astro-CB NN':           '#9ECAE1',
    'Astro-other NN':        '#AEC7E8',
    'Microglia NN':          '#D62728',
    'BAM NN':                '#FF7F0E',
    'DC NN':                 '#BCBD22',
    'OEC NN':                '#9467BD',
    # Non-glia broad
    'Excitatory':            '#2CA02C',
    'Inhibitory':            '#00897B',
    'Modulatory':            '#7F7F7F',
    'Vascular_Ependymal':    '#FFD700',
    'Endo NN':               '#FFD700',
    'Vascular NN':           '#DAA520',
    # Glutamatergic subtypes (AD fine-grained)
    'CA1-ProS Glut':         '#00897B',
    'CA3 Glut':              '#26A69A',
    'DG Glut':               '#80CBC4',
    'L2/3 IT CTX Glut':      '#2CA02C',
    'L4/5 IT CTX Glut':      '#41AB5D',
    'L5 IT CTX Glut':        '#17BECF',
    'L5 PPP Glut':           '#31A354',
    'L6 CT CTX Glut':        '#74C476',
    'L6 IT CTX Glut':        '#A1D99B',
    'L6b EPd Glut':          '#C7E9C0',
    'CT SUB Glut':           '#006D2C',
    # GABAergic subtypes
    'Sst Gaba':              '#9467BD',
    'Vip Gaba':              '#AD9EC0',
    'Pvalb chandelier Gaba': '#7B4F96',
    'STR D1 Gaba':           '#C9BEDE',
    'STR D2 Gaba':           '#B3A2C7',
    # Legacy underscore-key aliases (for backwards compat / direct lookups)
    'Micro-PVM':             '#D62728',
    'Microglia_NN':          '#D62728',
    'BAM_NN':                '#FF7F0E',
    'DC_NN':                 '#BCBD22',
    'Astro':                 '#1F77B4',
    'Oligo':                 '#8c564b',
    'Oligo_NN':              '#8c564b',
    'OPC':                   '#e377c2',
    'OPC_NN':                '#e377c2',
    'L2_3_IT_CTX':           '#2CA02C',
    'L5_IT_CTX':             '#17BECF',
    'CA1-ProS':              '#00897B',
    'DG':                    '#7F7F7F',
    'Endo':                  '#FFD700',
}
DEFAULT_CT_COLOR = '#888888'

# Differential fill colors
DIFF_TRT_COLOR    = '#e41a1c'    # trt-unique (e.g. TG)
DIFF_CTRL_COLOR   = '#377eb8'    # ctrl-unique (e.g. WT)
DIFF_SHARED_COLOR = '#aaaaaa'    # shared coverage


# ── GTF helpers ────────────────────────────────────────────────────────────────

def _open(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p)


def get_gene_coords(gtf_path, gene_name):
    """Return (chrom, gene_start0, gene_end, strand) for first matching gene."""
    pat = re.compile(rf'gene_name "{re.escape(gene_name)}"')
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            if pat.search(parts[8]):
                return parts[0], int(parts[3]) - 1, int(parts[4]), parts[6]
    sys.exit(f"[browser] gene '{gene_name}' not found in GTF")


def parse_merged_exons(gtf_path, gene_name, chrom):
    """Return sorted, merged exon intervals for gene_name on chrom.
    Adapted from SDas plot_srebf1_locus.py — single clean exon row."""
    pat = re.compile(rf'gene_name "{re.escape(gene_name)}"')
    exons = []
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'exon' or parts[0] != chrom:
                continue
            if pat.search(parts[8]):
                exons.append((int(parts[3]) - 1, int(parts[4])))
    if not exons:
        return []
    exons.sort()
    merged = [list(exons[0])]
    for s, e in exons[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def draw_gene_track(ax, gene_name, strand, merged_exons, gene_start, gene_end,
                    reg_start, reg_end):
    """Draw merged-exon gene structure — single clean row (SDas V2 style)."""
    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(-1.2, 0.6)
    ax.set_yticks([])
    ax.spines[['top', 'left', 'right']].set_visible(False)

    gs_clip = max(gene_start, reg_start)
    ge_clip = min(gene_end,   reg_end)

    # Backbone
    ax.plot([gs_clip, ge_clip], [0, 0],
            color='#333333', linewidth=1.5, solid_capstyle='butt', zorder=1)

    # Strand arrows along backbone
    span = ge_clip - gs_clip
    n_arrows = max(5, span // 5000)
    arrow_xs = np.linspace(gs_clip + span * 0.05, ge_clip - span * 0.05, n_arrows)
    dx = span * 0.015 * (-1 if strand == '-' else 1)
    for xp in arrow_xs:
        ax.annotate(
            '', xy=(xp + dx, 0), xytext=(xp, 0),
            arrowprops=dict(arrowstyle='-|>', color='#666666',
                            lw=0.6, mutation_scale=6),
        )

    # Merged exon blocks
    for (es, ee) in merged_exons:
        es = max(es, reg_start)
        ee = min(ee, reg_end)
        if ee > es:
            ax.add_patch(mpatches.Rectangle(
                (es, -0.38), ee - es, 0.76,
                facecolor='#2c2c2c', edgecolor='none', zorder=2,
            ))

    # Gene name + strand below backbone
    mid = (gs_clip + ge_clip) / 2
    ax.text(mid, -0.70, gene_name, ha='center', va='top',
            fontsize=8, fontstyle='italic', color='#333333')
    ax.text(mid, -1.00, f'({strand} strand)', ha='center', va='top',
            fontsize=6.5, color='#666666')

    ax.set_xlabel(f'{reg_start // 1_000_000 * 1_000_000 // 1}', fontsize=0)  # hidden
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x / 1e6:.3f}'))
    ax.tick_params(axis='x', labelsize=7)


# ── BigWig helpers ─────────────────────────────────────────────────────────────

def read_bw(path, chrom, start, end, n_bins=700):
    """Read BigWig values; returns float array length n_bins."""
    try:
        import pyBigWig
    except ImportError:
        sys.exit('[browser] pyBigWig not available in this container')
    try:
        bw   = pyBigWig.open(path)
        vals = bw.stats(chrom, start, end, type='mean', nBins=n_bins)
        bw.close()
        return np.array([v if v is not None else 0.0 for v in vals], dtype=float)
    except Exception as e:
        print(f'[browser][WARN] read_bw failed {path}: {e}', flush=True)
        return np.zeros(n_bins, dtype=float)


def compute_shared_ymax(paths, chrom, start, end, n_bins=500, percentile=99.5):
    """Global y-axis ceiling across a set of BigWig files."""
    global_max = 0.0
    for path in paths:
        vals = read_bw(path, chrom, start, end, n_bins)
        if vals.size:
            global_max = max(global_max, float(np.nanpercentile(vals, percentile)))
    result = math.ceil(global_max) if global_max > 0 else 10.0
    print(f'[browser] auto ymax={result} (raw_max={global_max:.3f})', flush=True)
    return float(result)


# ── Track renderers ────────────────────────────────────────────────────────────

def draw_absolute_track(ax, vals, reg_start, reg_end, color, ymax, label):
    """Single-condition fill track."""
    x = np.linspace(reg_start, reg_end, len(vals))
    ax.fill_between(x, 0, vals, color=color, alpha=0.85, linewidth=0)
    _style_track_ax(ax, label, reg_start, reg_end, ymax, ct_color=color)


def draw_differential_track(ax, ctrl_vals, trt_vals, reg_start, reg_end, ymax,
                             label, ctrl_label='WT', trt_label='TG'):
    """Differential fill: shared=grey, trt-unique=red, ctrl-unique=blue.
    Thin outline curves on top of fills — matching SDas plot_srebf1_locus.py V2."""
    x      = np.linspace(reg_start, reg_end, len(ctrl_vals))
    shared = np.minimum(ctrl_vals, trt_vals)

    # Shared region (grey)
    ax.fill_between(x, 0, shared, color=DIFF_SHARED_COLOR, alpha=0.80, linewidth=0,
                    zorder=1, label='shared')
    # TG-unique (red) — trt > ctrl
    mask_tg = trt_vals > ctrl_vals
    ax.fill_between(x, shared, trt_vals, where=mask_tg,
                    color=DIFF_TRT_COLOR, alpha=0.75, linewidth=0,
                    zorder=2, label=trt_label)
    # WT-unique (blue) — ctrl > trt
    mask_wt = ctrl_vals > trt_vals
    ax.fill_between(x, shared, ctrl_vals, where=mask_wt,
                    color=DIFF_CTRL_COLOR, alpha=0.75, linewidth=0,
                    zorder=2, label=ctrl_label)
    # Thin outline curves on top (SDas V2 style — help distinguish overlapping regions)
    ax.plot(x, trt_vals,  color=DIFF_TRT_COLOR,  linewidth=0.5, zorder=3)
    ax.plot(x, ctrl_vals, color=DIFF_CTRL_COLOR, linewidth=0.5, zorder=3)

    _style_track_ax(ax, label, reg_start, reg_end, ymax)


def _style_track_ax(ax, label, reg_start, reg_end, ymax, ct_color=None):
    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(0, ymax)
    # SDas V2 style: only the ymax tick, no zero tick cluttering the baseline
    ax.set_yticks([ymax])
    ax.set_yticklabels([f'{ymax:.0f}'], fontsize=5.5)
    ax.tick_params(axis='x', labelbottom=False, bottom=False)
    ax.tick_params(axis='y', length=2, pad=1)
    ax.spines[['top', 'right', 'bottom']].set_visible(False)
    color = ct_color if ct_color else CT_COLORS.get(label, DEFAULT_CT_COLOR)
    ax.set_ylabel(label, fontsize=6.5, rotation=0, labelpad=4,
                  va='center', ha='right', color=color)


def draw_xaxis(ax, reg_start, reg_end, chrom):
    """Minimal genomic x-axis."""
    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(0, 1)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=8)
    span = reg_end - reg_start
    # Choose a sensible tick interval
    for interval in [100, 500, 1000, 2000, 5000, 10000, 25000, 50000, 100000]:
        if span / interval < 10:
            break
    ticks = list(range(
        int(math.ceil(reg_start / interval)) * interval,
        reg_end, interval))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:,}' for t in ticks], rotation=30, ha='right')
    ax.text(0.5, -0.5, chrom, transform=ax.transAxes,
            ha='center', va='top', fontsize=8, color='#333333')


# ── DA peaks track ─────────────────────────────────────────────────────────────

def draw_da_peaks_track(ax, da_bed_path, reg_start, reg_end):
    """Tick marks for DA peak midpoints."""
    ax.set_xlim(reg_start, reg_end)
    ax.set_ylim(0, 1)
    ax.axis('off')
    try:
        with open(da_bed_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                ps, pe = int(parts[1]), int(parts[2])
                mid = (ps + pe) / 2
                if reg_start <= mid <= reg_end:
                    ax.axvspan(ps, pe, ymin=0.1, ymax=0.9,
                               color='#FF4500', alpha=0.4, linewidth=0)
        ax.text(-0.01, 0.5, 'DA', transform=ax.transAxes,
                ha='right', va='center', fontsize=6, color='#FF4500')
    except Exception as e:
        print(f'[browser][WARN] DA peaks load failed: {e}', flush=True)


# ── Manifest helpers ───────────────────────────────────────────────────────────

def load_manifest(manifest_path, bw_dir):
    """Return (ct_bw_map, by_condition_map) with absolute paths."""
    bw_dir = Path(bw_dir)
    with open(manifest_path) as f:
        m = json.load(f)

    ct_bw_map = {ct: str(bw_dir / fname)
                 for ct, fname in m.get('bigwigs', {}).items()}
    by_cond_raw = m.get('by_condition', {})
    by_cond_map = {
        ct: {cond: str(bw_dir / fname) for cond, fname in conds.items()}
        for ct, conds in by_cond_raw.items()
    }
    return ct_bw_map, by_cond_map


# ── Main ───────────────────────────────────────────────────────────────────────

def collect_bigwigs_from_dir(bw_dir, cell_types):
    """Absolute-mode fallback: resolve {ct}.bw files directly from a directory."""
    bw_dir = Path(bw_dir)
    ct_bw_map = {}
    for ct in cell_types:
        for stem in (ct, ct.replace(' ', '_'), ct.replace('_', ' ')):
            path = bw_dir / f'{stem}.bw'
            if path.exists():
                ct_bw_map[ct] = str(path)
                break
        else:
            print(f'[browser][WARN] no BigWig for {ct!r} in {bw_dir}', flush=True)
    return ct_bw_map


def parse_args():
    p = argparse.ArgumentParser(
        description='Render ATAC genome browser (matplotlib, no arcs)')
    # Required
    p.add_argument('--gene',             required=True)
    p.add_argument('--gtf',              required=True)
    p.add_argument('--bw-dir',           required=True,
                   help='Directory containing BigWig files. '
                        'With --bw-manifest, manifest paths are resolved here. '
                        'Without --bw-manifest (absolute mode only), scanned for {ct}.bw.')
    p.add_argument('--cell-types',       required=True,
                   help='Comma-separated ordered list of cell types to show')
    p.add_argument('--out-png',          required=True)
    # Optional
    p.add_argument('--bw-manifest',      default=None,
                   help='manifest.json from export_atac_bigwigs.py '
                        '(required for differential mode; auto-scanned for absolute)')
    p.add_argument('--mode',             default='absolute',
                   choices=['absolute', 'differential'])
    p.add_argument('--ctrl-condition',   default='WT',
                   help='Control condition label in manifest.by_condition (default WT)')
    p.add_argument('--trt-condition',    default='TG',
                   help='Treatment condition label in manifest.by_condition (default TG)')
    p.add_argument('--window',           type=int, default=15000,
                   help='bp padding around gene body (default 15000)')
    p.add_argument('--n-bins',           type=int, default=700,
                   help='BigWig resolution bins (default 700)')
    p.add_argument('--max-value',        type=float, default=None,
                   help='y-axis ceiling; auto-computed from 99.5th pct if omitted')
    p.add_argument('--da-peaks',         default=None,
                   help='Optional DA peaks BED for annotation strip')
    p.add_argument('--out-pdf',          default=None,
                   help='PDF path; auto-derived from --out-png if omitted')
    p.add_argument('--width',            type=float, default=8.0,
                   help='Figure width in inches (default 8.0)')
    return p.parse_args()


def main():
    args = parse_args()

    # ── Region ────────────────────────────────────────────────────────────────
    chrom, gs, ge, strand = get_gene_coords(args.gtf, args.gene)
    w      = args.window
    r_s    = max(0, gs - w)
    r_e    = ge + w
    print(f'[browser] {args.gene}: {chrom}:{gs}-{ge} ({strand}), '
          f'region: {chrom}:{r_s}-{r_e}', flush=True)

    cell_types = [ct.strip() for ct in args.cell_types.split(',') if ct.strip()]

    # ── Manifest ──────────────────────────────────────────────────────────────
    if args.mode == 'differential' and not args.bw_manifest:
        sys.exit('[browser] ERROR: --bw-manifest is required for --mode differential')

    if args.bw_manifest and Path(args.bw_manifest).exists():
        ct_bw_map, by_cond_map = load_manifest(args.bw_manifest, args.bw_dir)
    else:
        if args.bw_manifest:
            print(f'[browser][WARN] manifest not found: {args.bw_manifest}; '
                  f'falling back to directory scan', flush=True)
        ct_bw_map  = collect_bigwigs_from_dir(args.bw_dir, cell_types)
        by_cond_map = {}

    # ── Cell type defaulting ──────────────────────────────────────────────────
    # When --cell-types '' (empty), auto-select from the manifest rather than
    # rendering nothing. Differential mode picks CTs with both conditions present.
    if not cell_types:
        if args.mode == 'differential':
            cell_types = sorted(
                ct for ct, conds in by_cond_map.items()
                if args.ctrl_condition in conds and args.trt_condition in conds
            )
            print(f'[browser] --cell-types empty; auto-selected '
                  f'{len(cell_types)} differential CTs (both {args.ctrl_condition}'
                  f'/{args.trt_condition} present) from manifest', flush=True)
        else:
            cell_types = sorted(ct_bw_map.keys())
            print(f'[browser] --cell-types empty; auto-selected '
                  f'{len(cell_types)} CTs from manifest', flush=True)

    # Validate BigWig availability per mode
    tracks = []   # [(label, path_or_None_for_diff)]
    for ct in cell_types:
        if args.mode == 'differential':
            cond_map = by_cond_map.get(ct, {})
            ctrl_path = cond_map.get(args.ctrl_condition)
            trt_path  = cond_map.get(args.trt_condition)
            if not ctrl_path or not trt_path:
                print(f'[browser][WARN] {ct}: missing condition BigWigs '
                      f'({args.ctrl_condition}/{args.trt_condition}) in manifest — skip',
                      flush=True)
                continue
            if not Path(ctrl_path).exists() or not Path(trt_path).exists():
                print(f'[browser][WARN] {ct}: BigWig file(s) not on disk — skip', flush=True)
                continue
            tracks.append((ct, ctrl_path, trt_path))
        else:
            path = ct_bw_map.get(ct)
            if not path:
                print(f'[browser][WARN] {ct}: not in manifest — skip', flush=True)
                continue
            if not Path(path).exists():
                print(f'[browser][WARN] {ct}: file not on disk — skip', flush=True)
                continue
            tracks.append((ct, path))

    if not tracks:
        sys.exit('[browser] ERROR: no tracks available for rendering')
    print(f'[browser] rendering {len(tracks)} tracks', flush=True)

    # ── Shared y-max ──────────────────────────────────────────────────────────
    if args.max_value is not None:
        ymax = args.max_value
    else:
        all_paths = []
        for t in tracks:
            all_paths.extend(t[1:])   # absolute: (ct, path); diff: (ct, ctrl, trt)
        ymax = compute_shared_ymax(all_paths, chrom, r_s, r_e)

    # ── Figure layout ─────────────────────────────────────────────────────────
    n_ct  = len(tracks)
    has_da = args.da_peaks and Path(args.da_peaks).exists()

    row_h_track = 1.0
    row_h_da    = 0.3
    row_h_gene  = 2.5   # single merged-exon row (SDas style); was per-transcript multi-row
    row_h_xaxis = 0.6
    top_margin    = 0.7
    bottom_margin = 0.4
    left_frac  = 0.13
    right_frac = 0.97

    n_extra = (1 if has_da else 0)
    hr = ([row_h_da] if has_da else []) + [row_h_track] * n_ct + [row_h_gene, row_h_xaxis]
    fig_h = top_margin + bottom_margin + sum(hr)
    top_frac    = 1.0 - top_margin    / fig_h
    bottom_frac =       bottom_margin / fig_h

    gs_layout = gridspec.GridSpec(
        len(hr), 1, height_ratios=hr, hspace=0.03,
        left=left_frac, right=right_frac,
        top=top_frac, bottom=bottom_frac,
    )
    fig = plt.figure(figsize=(args.width, fig_h))

    row_idx = 0

    # DA peaks strip
    if has_da:
        ax_da = fig.add_subplot(gs_layout[row_idx])
        draw_da_peaks_track(ax_da, args.da_peaks, r_s, r_e)
        row_idx += 1

    # BigWig tracks
    ax_ref = None
    for t in tracks:
        ax = fig.add_subplot(gs_layout[row_idx], sharex=ax_ref)
        if ax_ref is None:
            ax_ref = ax
        ct_label = t[0]

        if args.mode == 'differential':
            _, ctrl_path, trt_path = t
            ctrl_vals = read_bw(ctrl_path, chrom, r_s, r_e, args.n_bins)
            trt_vals  = read_bw(trt_path,  chrom, r_s, r_e, args.n_bins)
            draw_differential_track(ax, ctrl_vals, trt_vals, r_s, r_e, ymax,
                                    ct_label, args.ctrl_condition, args.trt_condition)
        else:
            _, path = t
            vals  = read_bw(path, chrom, r_s, r_e, args.n_bins)
            color = CT_COLORS.get(ct_label, DEFAULT_CT_COLOR)
            draw_absolute_track(ax, vals, r_s, r_e, color, ymax, ct_label)

        row_idx += 1

    # Gene structure — merged-exon single row (SDas V2 style)
    merged_exons = parse_merged_exons(args.gtf, args.gene, chrom)
    ax_gene = fig.add_subplot(gs_layout[row_idx], sharex=ax_ref)
    draw_gene_track(ax_gene, args.gene, strand, merged_exons, gs, ge, r_s, r_e)
    row_idx += 1

    # x-axis
    ax_x = fig.add_subplot(gs_layout[row_idx], sharex=ax_ref)
    draw_xaxis(ax_x, r_s, r_e, chrom)

    # ── Title + legend ────────────────────────────────────────────────────────
    mode_tag = ' — differential (TG vs WT)' if args.mode == 'differential' else ''
    fig.suptitle(
        f'{args.gene} locus  ·  ATAC accessibility{mode_tag}  '
        f'·  {chrom}:{r_s:,}–{r_e:,}',
        fontsize=9, y=1.0 - 0.02 / fig_h,
    )

    if args.mode == 'differential':
        legend_elements = [
            mpatches.Patch(color=DIFF_SHARED_COLOR, label='shared'),
            mpatches.Patch(color=DIFF_TRT_COLOR,    label=args.trt_condition),
            mpatches.Patch(color=DIFF_CTRL_COLOR,   label=args.ctrl_condition),
        ]
        fig.legend(handles=legend_elements, loc='upper right',
                   bbox_to_anchor=(right_frac, top_frac - 0.01),
                   fontsize=7, frameon=False, ncol=3)

    # ── Output ────────────────────────────────────────────────────────────────
    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    print(f'[browser] wrote {out_png}', flush=True)

    out_pdf = Path(args.out_pdf) if args.out_pdf else out_png.with_suffix('.pdf')
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'[browser] wrote {out_pdf}', flush=True)

    plt.close(fig)


if __name__ == '__main__':
    main()
