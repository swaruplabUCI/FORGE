#!/usr/bin/env python3
"""
render_msfp_promoter_strip.py  —  FORGE pipeline module

Promoter MSFP strip: 150 bp zoom window with TF motif logo insets
and genomic coordinates. Renders one gene per row.

Differences from the oneOff 04_msfp_promoter.py:
  - No hardcoded path constants (RESULTS, OUT_DIR, DEFAULT_PFM,
    DEFAULT_CACHE, FOOTPRINTS_DIR) — all paths are passed as CLI args.
  - --pfm and --cache-dir are REQUIRED (no defaults).
  - --gtf is OPTIONAL; when omitted, TSS proximity filtering is skipped.
  - PDF is always written alongside PNG; --out-pdf overrides the auto-derived path.

Modes
-----
  absolute     : single-condition heatmap (arr[0])
  differential : Δ = arr[1] − arr[0]
  all_three    : ctrl | trt | Δ stacked per gene (requires 3-D arr with ≥2 slices)

Usage (Nextflow context):
    singularity exec snapatac_extended.sif python3 render_msfp_promoter_strip.py \\
        --scan-dir      results/scprinter/footprints/Microglia_NN \\
        --genes         Trem2,Apoe \\
        --tfs           Spi1 \\
        --pfm           /ref/scprinter/JASPAR2022_core_nonredundant.jaspar \\
        --cache-dir     /ref/scprinter \\
        --mode          all_three \\
        --out-png       msfp_promoter_Spi1_Microglia.png
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import scprinter as scp
import logomaker


# ── Shared genomic helpers ────────────────────────────────────────────────────

def parse_locus(locus_str):
    chrom, rng = locus_str.split(':')
    s, e = rng.split('-')
    return chrom, int(s), int(e)


def scan_motifs_in_window(chrom, start, end, pfm_path, cache_dir, genome, tfs):
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    genome_obj = getattr(scp.genome, genome)
    scanner = scp.motifs.Motifs(
        ref_path_motif=pfm_path,
        ref_path_fa=genome_obj.fetch_fa(),
        bg='even',
        n_jobs=1,
        motif_name_func=lambda x: x.split('\t')[-1],
    )
    try:
        scanner.prep_scanner(tf_genes=list(tfs))
        hits = scanner.scan_motif([(chrom, start, end)], verbose=False, split_tfs=True)
    finally:
        del scanner
    return hits


def parse_hits(hits, chrom, win_start):
    out = []
    sample_printed = False
    for h in hits:
        if not sample_printed:
            print(f'  raw hit sample: {h!r}', flush=True)
            sample_printed = True
        try:
            if h[0] != chrom:
                continue
            peak_s  = int(h[1])
            tf_name = h[4]
            score   = float(h[5])
            strand  = int(h[6])
            abs_s   = peak_s + int(h[7])
            abs_e   = peak_s + int(h[8])
            out.append((abs_s - win_start, abs_e - win_start, tf_name, score, strand))
        except (IndexError, ValueError, TypeError) as exc:
            print(f'  hit parse fail: {h!r} ({exc})', flush=True)
    return out


# ── JASPAR logo helpers ───────────────────────────────────────────────────────

def parse_jaspar_pfm(pfm_path, tf_name_substr):
    rows_acgt = {b: [] for b in 'ACGT'}
    inside = False
    with open(pfm_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if tf_name_substr.upper() in line.upper():
                    inside = True
                    rows_acgt = {b: [] for b in 'ACGT'}
                elif inside:
                    break
                continue
            if not inside:
                continue
            for b in 'ACGT':
                if line.strip().startswith(b):
                    nums = [float(x) for x in line.replace('[', '').replace(']', '').split()
                            if x not in (b,)]
                    rows_acgt[b] = nums
                    break
    arr = np.array([rows_acgt[b] for b in 'ACGT'], dtype=float)
    col_sums = arr.sum(axis=0)
    col_sums[col_sums == 0] = 1.0
    freq = (arr / col_sums).T
    return pd.DataFrame(freq, columns=list('ACGT'))


def draw_motif_logo(ax_parent, pfm_path, tf_name_substr, rel_start, rel_end, zoom_w, strand):
    freq_df = parse_jaspar_pfm(pfm_path, tf_name_substr)
    if freq_df.empty:
        return
    if strand == -1:
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        freq_df = (freq_df[['A', 'C', 'G', 'T']]
                   .rename(columns=comp)
                   .iloc[::-1]
                   .reset_index(drop=True))
    x0_frac = max(0.0, rel_start / zoom_w)
    x1_frac = min(1.0, rel_end   / zoom_w)
    if x1_frac <= x0_frac:
        return
    inset = ax_parent.inset_axes([x0_frac, 0.38, x1_frac - x0_frac, 0.60])
    try:
        logomaker.Logo(freq_df, ax=inset, color_scheme='classic', show_spines=False)
    except Exception as exc:
        print(f'  logomaker failed ({tf_name_substr}): {exc}', flush=True)
        return
    inset.set_xticks([])
    inset.set_yticks([])
    inset.patch.set_alpha(0)


# ── Reference sequence helpers ────────────────────────────────────────────────

_COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


def draw_seq_track(ax, seq, zoom_w):
    """Plot reference DNA sequence as per-base colored letters."""
    _COLORS = {'A': '#2ca02c', 'C': '#1f77b4', 'G': '#ff7f0e',
               'T': '#d62728', 'N': '#aaaaaa'}
    ax.set_xlim(0, zoom_w)
    ax.set_ylim(0, 1)
    ax.axis('off')
    n = len(seq)
    fs = max(3.5, min(6.0, 150.0 / n))
    for i, base in enumerate(seq):
        ax.text(i + 0.5, 0.5, base.upper(),
                ha='center', va='center',
                fontsize=fs, fontfamily='monospace', fontweight='bold',
                color=_COLORS.get(base.upper(), '#aaaaaa'),
                transform=ax.transData, clip_on=True)


# ── Gene data loading ─────────────────────────────────────────────────────────

def load_gene_promoter(gene, h5ad_path, pfm_path, cache_dir, genome, tfs,
                       max_scale, zoom_half_width, mode='absolute'):
    """Load scPrinter footprint h5ad and return plotting dict."""
    print(f'\n=== loading {gene} ({mode}) ===', flush=True)
    fp = ad.read_h5ad(h5ad_path)
    locus = list(fp.obsm.keys())[0]
    chrom, win_start, win_end = parse_locus(locus)
    win_width = win_end - win_start

    arr    = np.asarray(fp.obsm[locus])
    scales = np.asarray(fp.uns['scales'])
    print(f'  arr.shape={arr.shape}  n_scales={len(scales)}', flush=True)

    msfp_ctrl_raw = None
    msfp_trt_raw  = None

    if arr.ndim == 3:
        if mode == 'all_three':
            if arr.shape[0] >= 2:
                msfp_ctrl_raw = arr[0]
                msfp_trt_raw  = arr[1]
                msfp = arr[1] - arr[0]   # Δ stored as primary for vmax/rendering
            else:
                print(f'  {gene}: only {arr.shape[0]} condition slice(s) — '
                      f'all_three requires ≥2; falling back to absolute.', flush=True)
                msfp = arr[0]
        elif mode == 'differential' and arr.shape[0] >= 2:
            msfp = arr[1] - arr[0]
        else:
            msfp = arr[0]
    elif arr.ndim == 2:
        if mode in ('differential', 'all_three'):
            print(f'  {gene}: arr is 2-D (no condition axis) — '
                  f'{mode} mode not available; using as-is.', flush=True)
        msfp = arr
    else:
        print(f'  {gene}: unexpected arr shape {arr.shape}; skipping', flush=True)
        return None

    scale_mask       = scales <= max_scale
    scales_plot      = scales[scale_mask]
    msfp_masked      = msfp[scale_mask]
    msfp_ctrl_masked = msfp_ctrl_raw[scale_mask] if msfp_ctrl_raw is not None else None
    msfp_trt_masked  = msfp_trt_raw[scale_mask]  if msfp_trt_raw  is not None else None

    raw_hits   = scan_motifs_in_window(chrom, win_start, win_end,
                                       pfm_path, cache_dir, genome, tfs)
    all_hits   = parse_hits(raw_hits, chrom, win_start)
    tfs_upper  = [t.strip().upper() for t in tfs]
    query_hits = [(hs, he, tf, sc, st) for hs, he, tf, sc, st in all_hits
                  if tf.upper() in tfs_upper]

    if not query_hits:
        print(f'  {gene}: no motif found for {tfs} — skipping', flush=True)
        return None

    best_hit    = max(query_hits, key=lambda x: x[3])
    zoom_center = (best_hit[0] + best_hit[1]) // 2
    z0 = max(0, zoom_center - zoom_half_width)
    z1 = min(win_width, zoom_center + zoom_half_width)
    zoom_w = z1 - z0

    msfp_zoom      = msfp_masked[:, z0:z1]
    msfp_zoom_ctrl = msfp_ctrl_masked[:, z0:z1] if msfp_ctrl_masked is not None else None
    msfp_zoom_trt  = msfp_trt_masked[:, z0:z1]  if msfp_trt_masked  is not None else None
    zoom_hits = [(hs - z0, he - z0, tf, sc, st)
                 for hs, he, tf, sc, st in all_hits
                 if hs >= z0 and he <= z1]

    strand_char = '+' if best_hit[4] >= 0 else '-'
    abs_left  = win_start + z0
    abs_right = win_start + z1
    print(f'  {gene}: PASS  motif={best_hit[2]}  score={best_hit[3]:.2f}  '
          f'strand={strand_char}  window={chrom}:{abs_left}-{abs_right}',
          flush=True)

    # Fetch reference sequence in TF-frame orientation.
    try:
        genome_obj = getattr(scp.genome, genome)
        raw_seq = genome_obj.fetch_seq(chrom, abs_left, abs_right)
        ref_seq = _revcomp(raw_seq) if strand_char == '-' else raw_seq
    except Exception as e:
        print(f'  [WARN] fetch_seq failed: {e}', flush=True)
        ref_seq = 'N' * zoom_w

    return {
        'gene':           gene,
        'locus':          locus,
        'chrom':          chrom,
        'win_start':      win_start,
        'z0': z0, 'z1':  z1,
        'zoom_w':         zoom_w,
        'abs_left':       abs_left,
        'abs_right':      abs_right,
        'msfp_zoom':      msfp_zoom,       # Δ in all_three, else absolute/differential
        'msfp_zoom_ctrl': msfp_zoom_ctrl,  # ctrl slice (all_three only, else None)
        'msfp_zoom_trt':  msfp_zoom_trt,   # trt slice  (all_three only, else None)
        'scales_plot':    scales_plot,
        'zoom_hits':      zoom_hits,
        'best_hit':       best_hit,
        'strand':         strand_char,
        'ref_seq':        ref_seq,
    }


# ── Figure renderer ───────────────────────────────────────────────────────────

def render_strip(gene_data_list, args, shared_vmax, cell_type_label,
                 shared_abs_vmax=None, shared_diff_vmax=None):
    """Render promoter MSFP strip.

    Modes
    -----
    absolute / differential : 3 rows per gene (msfp | logo | seq), one colorbar.
    all_three               : 5 rows per gene (ctrl | trt | Δ | logo | seq),
                              two colorbars (abs top, diff bottom).
    """
    N    = len(gene_data_list)
    mode = getattr(args, 'mode', 'absolute')
    all_three = (mode == 'all_three')

    # Row heights (inches)
    row_h_msfp    = 1.8 if all_three else 2.8
    row_h_logo    = 0.55
    row_h_seq     = 0.32
    rows_per_gene = 5 if all_three else 3
    top_margin    = 0.9
    bottom_margin = 0.6

    msfp_rows_per_gene = 3 if all_three else 1
    fig_h = (top_margin + bottom_margin
             + N * (msfp_rows_per_gene * row_h_msfp + row_h_logo + row_h_seq))
    fig_w = 7.0

    left_col  = 0.14
    right_col = 0.85 if all_three else 0.88
    cbar_x0   = 0.87 if all_three else 0.91
    cbar_w    = 0.022

    top_frac    = 1.0 - top_margin    / fig_h
    bottom_frac =       bottom_margin / fig_h

    hr = []
    for _ in range(N):
        if all_three:
            hr.extend([row_h_msfp, row_h_msfp, row_h_msfp, row_h_logo, row_h_seq])
        else:
            hr.extend([row_h_msfp, row_h_logo, row_h_seq])

    gs = gridspec.GridSpec(
        N * rows_per_gene, 1, height_ratios=hr, hspace=0.04,
        left=left_col, right=right_col,
        top=top_frac, bottom=bottom_frac,
    )

    fig    = plt.figure(figsize=(fig_w, fig_h))
    ax_ref = None
    tfs_upper = [t.strip().upper() for t in args.tfs]

    for k, d in enumerate(gene_data_list):
        base   = k * rows_per_gene
        zoom_w = d['zoom_w']
        scales = d['scales_plot']

        # ── all_three layout: ctrl | trt | Δ ──────────────────────────────────
        if all_three:
            ax_ctrl = fig.add_subplot(gs[base],     sharex=ax_ref)
            ax_trt  = fig.add_subplot(gs[base + 1],
                                      sharex=ax_ref if ax_ref is not None else ax_ctrl)
            ax_delt = fig.add_subplot(gs[base + 2],
                                      sharex=ax_ref if ax_ref is not None else ax_ctrl)
            ax_l    = fig.add_subplot(gs[base + 3],
                                      sharex=ax_ref if ax_ref is not None else ax_ctrl)
            ax_s    = fig.add_subplot(gs[base + 4],
                                      sharex=ax_ref if ax_ref is not None else ax_ctrl)
            if ax_ref is None:
                ax_ref = ax_ctrl

            abs_vmax  = shared_abs_vmax  if shared_abs_vmax  is not None else shared_vmax
            diff_vmax = shared_diff_vmax if shared_diff_vmax is not None else shared_vmax

            ctrl_arr = (d.get('msfp_zoom_ctrl') if d.get('msfp_zoom_ctrl') is not None
                        else d['msfp_zoom'])
            trt_arr  = (d.get('msfp_zoom_trt')  if d.get('msfp_zoom_trt')  is not None
                        else d['msfp_zoom'])
            delt_arr = d['msfp_zoom']   # precomputed Δ = trt − ctrl

            _yticks = [s for s in [2, 5, 10, 15, 20, 25, 30] if s <= scales[-1]]
            ctrl_label = getattr(args, 'control_condition',   'ctrl')
            trt_label  = getattr(args, 'treatment_condition', 'trt')
            for ax_row, arr_row, vmax_row, tag in [
                (ax_ctrl, ctrl_arr, abs_vmax,  ctrl_label),
                (ax_trt,  trt_arr,  abs_vmax,  trt_label),
                (ax_delt, delt_arr, diff_vmax, f'Δ {trt_label}−{ctrl_label}'),
            ]:
                ax_row.imshow(
                    arr_row, aspect='auto', origin='lower',
                    cmap='RdBu_r', vmin=-vmax_row, vmax=vmax_row,
                    interpolation='bilinear',
                    extent=[0, zoom_w, scales[0], scales[-1]],
                )
                ax_row.set_yticks(_yticks)
                ax_row.tick_params(axis='y', labelsize=6)
                ax_row.tick_params(axis='x', labelbottom=False)
                ax_row.text(0.01, 0.88, tag,
                            transform=ax_row.transAxes,
                            fontsize=7, va='top', ha='left',
                            color='#222222', style='italic', zorder=5)

            ax_ctrl.set_ylabel(d['gene'], fontsize=10, rotation=0,
                               labelpad=60, va='center', fontweight='bold')
            if k == N - 1:
                ax_delt.text(1.01, 0.5, 'scale (bp)', fontsize=7, rotation=90,
                             ha='left', va='center', transform=ax_delt.transAxes)
            if k == 0:
                ax_ctrl.set_title(
                    f"Promoter MSFP  ·  {cell_type_label}  ·  ctrl | trt | Δ",
                    fontsize=9, loc='center', pad=3)

        # ── absolute / differential layout ────────────────────────────────────
        else:
            ax_m = fig.add_subplot(gs[base], sharex=ax_ref)
            ax_l = fig.add_subplot(gs[base + 1],
                                   sharex=ax_ref if ax_ref is not None else ax_m)
            ax_s = fig.add_subplot(gs[base + 2],
                                   sharex=ax_ref if ax_ref is not None else ax_m)
            if ax_ref is None:
                ax_ref = ax_m

            ax_m.imshow(
                d['msfp_zoom'], aspect='auto', origin='lower',
                cmap='RdBu_r', vmin=-shared_vmax, vmax=shared_vmax,
                interpolation='bilinear',
                extent=[0, zoom_w, scales[0], scales[-1]],
            )
            ax_m.set_yticks([s for s in [2, 5, 10, 15, 20, 25, 30] if s <= scales[-1]])
            ax_m.tick_params(axis='y', labelsize=7)
            ax_m.tick_params(axis='x', labelbottom=False)
            ax_m.set_ylabel(d['gene'], fontsize=10, rotation=0,
                            labelpad=60, va='center', fontweight='bold')
            if k == 0:
                ax_m.set_title(
                    f"Promoter MSFP  ·  {cell_type_label}"
                    + ('  [Δ TG−WT]' if mode == 'differential' else ''),
                    fontsize=9, loc='center', pad=3)
            if k == N - 1:
                ax_m.text(1.01, 0.5, 'scale (bp)', fontsize=7, rotation=90,
                          ha='left', va='center', transform=ax_m.transAxes)

        # ── Logo row (shared across modes) ────────────────────────────────────
        ax_l.set_ylim(0, 1)
        ax_l.set_xlim(0, zoom_w)
        ax_l.set_yticks([])
        ax_l.tick_params(axis='x', bottom=False, labelbottom=False)
        ax_l.spines[:].set_visible(False)

        for hs, he, tf, sc_score, st in d['zoom_hits']:
            tf_key = tf.upper() if tf.upper() in tfs_upper else tf
            draw_motif_logo(ax_l, args.pfm, tf_key, hs, he, zoom_w, st)

        chev = '▶▶' if d['strand'] == '+' else '◀◀'
        ax_l.text(0.01, 0.10, f"{d['chrom']}:{d['abs_left']:,} {chev}",
                  ha='left', va='bottom', fontsize=7, color='#333333',
                  transform=ax_l.transAxes, zorder=3)
        ax_l.text(0.99, 0.10, f"{chev} {d['abs_right']:,}",
                  ha='right', va='bottom', fontsize=7, color='#333333',
                  transform=ax_l.transAxes, zorder=3)

        # ── Reference sequence track ──────────────────────────────────────────
        draw_seq_track(ax_s, d.get('ref_seq', ''), zoom_w)
        strand_label = f"ref seq  ({d['strand']} strand)"
        ax_s.text(-0.01, 0.5, strand_label, ha='right', va='center',
                  fontsize=6, color='#555555', transform=ax_s.transAxes)

    # ── Colorbars ─────────────────────────────────────────────────────────────
    total_h = top_frac - bottom_frac - 0.10

    if all_three:
        abs_vmax  = shared_abs_vmax  if shared_abs_vmax  is not None else shared_vmax
        diff_vmax = shared_diff_vmax if shared_diff_vmax is not None else shared_vmax
        half_h    = total_h / 2 - 0.02

        sm_abs = plt.cm.ScalarMappable(
            cmap=plt.get_cmap('RdBu_r'),
            norm=mcolors.Normalize(vmin=-abs_vmax, vmax=abs_vmax),
        )
        sm_abs.set_array([])
        cbar_ax_abs = fig.add_axes(
            [cbar_x0, bottom_frac + 0.05 + half_h + 0.04, cbar_w, half_h])
        fig.colorbar(sm_abs, cax=cbar_ax_abs, label='MSFP score')
        cbar_ax_abs.tick_params(labelsize=7)

        sm_diff = plt.cm.ScalarMappable(
            cmap=plt.get_cmap('RdBu_r'),
            norm=mcolors.Normalize(vmin=-diff_vmax, vmax=diff_vmax),
        )
        sm_diff.set_array([])
        cbar_ax_diff = fig.add_axes(
            [cbar_x0, bottom_frac + 0.05, cbar_w, half_h])
        fig.colorbar(sm_diff, cax=cbar_ax_diff, label='ΔMSFP (TG − WT)')
        cbar_ax_diff.tick_params(labelsize=7)
    else:
        sm = plt.cm.ScalarMappable(
            cmap=plt.get_cmap('RdBu_r'),
            norm=mcolors.Normalize(vmin=-shared_vmax, vmax=shared_vmax),
        )
        sm.set_array([])
        cbar_ax = fig.add_axes([cbar_x0, bottom_frac + 0.05, cbar_w, total_h])
        cbar_label = 'ΔMSFP (TG − WT)' if mode == 'differential' else 'MSFP score'
        fig.colorbar(sm, cax=cbar_ax, label=cbar_label)
        cbar_ax.tick_params(labelsize=7)

    mode_tag = {
        'differential': '  ·  Δ TG − WT',
        'all_three':    '  ·  ctrl | trt | Δ',
    }.get(mode, '')
    fig.suptitle(
        f"Promoter MSFP  ·  {cell_type_label}  ·  {N} genes  "
        f"·  {args.zoom_half_width * 2} bp window{mode_tag}",
        fontsize=10, y=1.0 - 0.04 / fig_h,
    )
    return fig


# ── h5ad discovery ────────────────────────────────────────────────────────────

def _find_h5ads(scan_dir, target_genes):
    d = Path(scan_dir)
    result = []
    for gene in target_genes:
        matches = sorted(d.glob(f'footprints_*{gene}*.h5ad'))
        exact = [m for m in matches if m.stem.endswith(f'_{gene}')]
        hits  = exact if exact else matches
        if not hits:
            print(f'[WARN] no h5ad for {gene} in {scan_dir}', flush=True)
        else:
            result.append((gene, str(hits[0])))
            if len(hits) > 1:
                print(f'  {gene}: {len(hits)} matches, using {hits[0].name}', flush=True)
    return result


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Render promoter MSFP strip (pipeline module)')
    # ── Required ──────────────────────────────────────────────────────────────
    ap.add_argument('--scan-dir',            required=True,
                    help='Directory containing footprints_*.h5ad files for one CT')
    ap.add_argument('--genes',               required=True,
                    help='Comma-separated gene names')
    ap.add_argument('--tfs',                 required=True,
                    help='Comma-separated TF motif names (e.g. Spi1,SPIB)')
    ap.add_argument('--pfm',                 required=True,
                    help='JASPAR .jaspar PFM file (no default in pipeline mode)')
    ap.add_argument('--cache-dir',           required=True,
                    help='scPrinter dataset cache directory (no default in pipeline mode)')
    ap.add_argument('--out-png',             required=True)
    # ── Optional ──────────────────────────────────────────────────────────────
    ap.add_argument('--cell-type',           default='',
                    help='Cell type label for figure title')
    ap.add_argument('--mode',                default='absolute',
                    choices=['absolute', 'differential', 'all_three'])
    ap.add_argument('--control-condition',   default='WT',
                    help='Label for condition slice 0 (default: WT)')
    ap.add_argument('--treatment-condition', default='TG',
                    help='Label for condition slice 1 (default: TG)')
    ap.add_argument('--genome',              default='mm10')
    ap.add_argument('--max-scale',           type=int, default=30)
    ap.add_argument('--zoom-half-width',     type=int, default=75)
    ap.add_argument('--out-pdf',             default=None,
                    help='PDF output path; auto-derived from --out-png if omitted')
    args = ap.parse_args()

    target_genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    tfs          = tuple(t.strip() for t in args.tfs.split(',') if t.strip())
    args.tfs     = tfs

    print(f'[setup] genes={target_genes}  tfs={tfs}  '
          f'zoom={args.zoom_half_width*2}bp  max_scale={args.max_scale}  '
          f'mode={args.mode}', flush=True)

    gene_pairs = _find_h5ads(args.scan_dir, target_genes)
    if not gene_pairs:
        sys.exit('ERROR: no footprint h5ad files found')

    passed, skipped = [], []
    for gene, h5ad_path in gene_pairs:
        try:
            d = load_gene_promoter(
                gene=gene, h5ad_path=h5ad_path,
                pfm_path=args.pfm, cache_dir=args.cache_dir,
                genome=args.genome, tfs=tfs,
                max_scale=args.max_scale,
                zoom_half_width=args.zoom_half_width,
                mode=args.mode,
            )
            (passed if d is not None else skipped).append(d if d is not None else gene)
        except Exception as exc:
            import traceback
            print(f'[WARN] {gene} error: {exc}', flush=True)
            traceback.print_exc()
            skipped.append(gene)

    print(f'\n[screen] PASSED  {len(passed)}: {[d["gene"] for d in passed]}', flush=True)
    print(f'[screen] SKIPPED {len(skipped)}: {skipped}', flush=True)

    if not passed:
        sys.exit('ERROR: no genes passed motif screening.')

    # ── Shared vmax computation ───────────────────────────────────────────────
    if args.mode == 'all_three':
        abs_arrays  = ([d['msfp_zoom_ctrl'] for d in passed
                        if d.get('msfp_zoom_ctrl') is not None]
                       + [d['msfp_zoom_trt'] for d in passed
                          if d.get('msfp_zoom_trt') is not None])
        diff_arrays = [d['msfp_zoom'] for d in passed]   # Δ stored here for all_three

        shared_abs_vmax  = (max(float(np.percentile(np.abs(a), 99)) for a in abs_arrays)
                            if abs_arrays else 1.0)
        shared_diff_vmax = max(float(np.percentile(np.abs(a), 99)) for a in diff_arrays)
        shared_abs_vmax  = max(shared_abs_vmax,  1e-6)
        shared_diff_vmax = max(shared_diff_vmax, 1e-6)
        shared_vmax      = shared_abs_vmax
        print(f'[render] shared_abs_vmax={shared_abs_vmax:.4f}  '
              f'shared_diff_vmax={shared_diff_vmax:.4f}', flush=True)
    else:
        shared_abs_vmax  = None
        shared_diff_vmax = None
        shared_vmax = max(float(np.percentile(np.abs(d['msfp_zoom']), 99))
                          for d in passed)
        shared_vmax = max(shared_vmax, 1e-6)
        print(f'[render] shared_vmax={shared_vmax:.4f}', flush=True)

    ct_label = args.cell_type or Path(args.scan_dir).name
    fig = render_strip(passed, args, shared_vmax, ct_label,
                       shared_abs_vmax=shared_abs_vmax,
                       shared_diff_vmax=shared_diff_vmax)

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    print(f'[render] wrote {out_png}', flush=True)

    # PDF is always written; --out-pdf overrides the auto-derived path
    out_pdf = Path(args.out_pdf) if args.out_pdf else out_png.with_suffix('.pdf')
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'[render] wrote {out_pdf}', flush=True)

    plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
