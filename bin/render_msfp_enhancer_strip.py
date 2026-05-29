#!/usr/bin/env python3
"""
render_msfp_enhancer_strip.py — pipeline-ready enhancer MSFP strip renderer.

Promotes oneOff/AD_Mm_10X_r5/05_msfp_enhancer.py to src_FORGE/bin.
Generates a publication-quality strip figure: one row per enhancer region,
each row showing a 150 bp zoom window centred on the best TF-motif hit,
with a JASPAR PWM logo inset and a per-base reference-sequence track.

Input h5ad format (from TARGETED_ENHANCER_FOOTPRINTING or
ENHANCER_FOOTPRINTING_PER_CT):
  obsm  — {region_id: ndarray}
            ndarray shape (n_scales, n_positions)       → absolute only
                    shape (n_conditions, n_scales, n_positions) → supports differential
  uns   — {'scales': array_of_bp_values}

Region selection filters (at least one must be satisfied per region):
  --context-bp        TSS proximity window around --target-gene
  --cicero-connections  Cicero co-accessibility links to --target-gene

Modes (--mode):
  absolute      First condition slice (arr[0]).  Default.
  differential  arr[1] - arr[0] (treatment minus control delta-MSFP).
                Requires ≥2 condition slices; falls back to absolute with a
                warning if the h5ad has only one slice.

Output:
  Both PNG (--out-png) and PDF are always written.  PDF path defaults to
  --out-png with the extension replaced (.pdf); pass --out-pdf to override.

Usage (pipeline — called by render_msfp_enhancer_strip.nf):
    python render_msfp_enhancer_strip.py \\
        --enhancer-h5ad enhancer_footprints_Microglia_NN_Spi1.h5ad \\
        --tfs Spi1 --target-gene Trem2 \\
        --gtf gencode.vM10.annotation.gtf \\
        --pfm JASPAR2022_core_nonredundant.jaspar \\
        --cache-dir /path/to/scprinter_cache \\
        --genome mm10 --mode differential \\
        --out-png msfp_enhancer_Spi1_Microglia_NN_diff.png
"""

import argparse
import gzip
import re
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


# ── Motif scan + logo helpers ─────────────────────────────────────────────────

def parse_locus(locus_str):
    chrom, rng = locus_str.split(':')
    s, e = rng.split('-')
    return chrom, int(s), int(e)


def scan_motifs_in_window(chrom, start, end, pfm_path, cache_dir, genome, tfs):
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    genome_obj = getattr(scp.genome, genome)
    # scPrinter Motifs scanner holds open file handles to the FASTA and index.
    # Always delete the scanner object explicitly so those handles are released
    # before the function returns — leaving them open risks file corruption on
    # shared DFS mounts if the process is killed mid-flight.
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
    for h in hits:
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
        except (IndexError, ValueError, TypeError):
            pass
    return out


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


def draw_motif_logo(ax_parent, pfm_path, tf_name_substr, rel_start, rel_end,
                    zoom_w, strand):
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


# ── GTF helper ────────────────────────────────────────────────────────────────

def get_gene_tss(gtf_path, gene_name):
    """Return (chrom, tss) for gene, or (None, None)."""
    pat = re.compile(rf'gene_name "{re.escape(gene_name)}"')
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            if pat.search(parts[8]):
                chrom  = parts[0]
                strand = parts[6]
                tss = int(parts[3]) - 1 if strand == '+' else int(parts[4])
                return chrom, tss
    return None, None


# ── Cicero proximity filter ───────────────────────────────────────────────────

def get_cicero_peaks_for_gene(cicero_gz, target_gene, gtf_path, slop=2000):
    """Return set of peak IDs linked to target_gene's promoter in Cicero."""
    gene_chrom, tss = get_gene_tss(gtf_path, target_gene)
    if not tss:
        return set()

    prom_s = tss - slop
    prom_e = tss + slop
    pat    = re.compile(r'(chr[^:]+):(\d+)-(\d+)')
    linked = set()

    try:
        with gzip.open(cicero_gz, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                for pk in (parts[0], parts[1]):
                    m = pat.match(pk)
                    if m and m.group(1) == gene_chrom:
                        ps, pe = int(m.group(2)), int(m.group(3))
                        if ps < prom_e and pe > prom_s:
                            other = parts[1] if pk == parts[0] else parts[0]
                            linked.add(other)
    except Exception as e:
        print(f'[msfp_enh][WARN] Cicero load failed: {e}', flush=True)

    return linked


# ── Reference sequence helpers ────────────────────────────────────────────────

_COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

def _revcomp(seq):
    return seq.translate(_COMP)[::-1]


def draw_seq_track(ax, seq, zoom_w):
    """Plot reference DNA sequence as per-base colored letters.

    Sequence is already in TF-frame orientation (revcomped for − strand before
    calling).  Each base is centered at x = i + 0.5 so columns align exactly
    with the MSFP heatmap and logo rows (shared x-axis, range [0, zoom_w]).

    Color convention: A=green, C=blue, G=orange, T=red, N=grey.
    """
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


# ── Core loading ──────────────────────────────────────────────────────────────

def load_enhancer_regions(h5ad_path, tfs, pfm_path, cache_dir, genome,
                          max_scale, zoom_half_width,
                          gene_chrom, gene_tss, context_bp,
                          cicero_linked_peaks, mode='absolute'):
    """Load enhancer h5ad and return list of plot dicts for passing regions.

    Parameters
    ----------
    mode : {'absolute', 'differential'}
        'absolute'     — use the first condition slice (arr[0]).
        'differential' — compute arr[1] - arr[0] (treatment minus control).
                         Falls back to arr[0] with a warning when the h5ad
                         carries fewer than two condition slices.
    """
    print(f'[msfp_enh] loading {h5ad_path}  mode={mode}', flush=True)
    fp = ad.read_h5ad(h5ad_path)

    region_keys = list(fp.obsm.keys())
    print(f'  {len(region_keys)} enhancer regions', flush=True)

    scales     = np.asarray(fp.uns.get('scales', []))
    scale_mask = scales <= max_scale
    if not scale_mask.any():
        print('[msfp_enh][WARN] no scales <= max_scale', flush=True)
        return []

    scales_plot = scales[scale_mask]
    tfs_upper   = [t.strip().upper() for t in tfs]
    passed      = []

    for region_key in region_keys:
        try:
            chrom, win_start, win_end = parse_locus(region_key)
        except Exception:
            continue

        # Filter: TSS proximity OR Cicero linkage
        in_context = (gene_chrom and chrom == gene_chrom and
                      gene_tss is not None and
                      abs((win_start + win_end) // 2 - gene_tss) <= context_bp)
        in_cicero  = region_key in cicero_linked_peaks

        if not in_context and not in_cicero:
            continue

        arr = np.asarray(fp.obsm[region_key])
        if arr.ndim == 3:
            if mode == 'differential':
                if arr.shape[0] >= 2:
                    msfp = arr[1] - arr[0]
                else:
                    print(f'  [WARN] {region_key}: only {arr.shape[0]} condition '
                          f'slice(s) in h5ad — differential requires ≥2; '
                          f'falling back to arr[0].', flush=True)
                    msfp = arr[0]
            else:
                msfp = arr[0]
        elif arr.ndim == 2:
            if mode == 'differential':
                print(f'  [WARN] {region_key}: arr is 2-D (no condition axis) — '
                      f'differential mode not available; using as-is.', flush=True)
            msfp = arr
        else:
            continue

        msfp_masked = msfp[scale_mask]

        # Motif scan
        try:
            raw_hits   = scan_motifs_in_window(chrom, win_start, win_end,
                                               pfm_path, cache_dir, genome, tfs)
            all_hits   = parse_hits(raw_hits, chrom, win_start)
            query_hits = [(hs, he, tf, sc, st) for hs, he, tf, sc, st in all_hits
                          if tf.upper() in tfs_upper]
        except Exception as e:
            print(f'  [WARN] motif scan failed for {region_key}: {e}', flush=True)
            continue

        if not query_hits:
            continue

        best_hit    = max(query_hits, key=lambda x: x[3])
        win_width   = win_end - win_start
        zoom_center = (best_hit[0] + best_hit[1]) // 2
        z0 = max(0, zoom_center - zoom_half_width)
        z1 = min(win_width, zoom_center + zoom_half_width)
        zoom_w = z1 - z0

        msfp_zoom = msfp_masked[:, z0:z1]
        zoom_hits = [(hs - z0, he - z0, tf, sc, st)
                     for hs, he, tf, sc, st in all_hits
                     if hs >= z0 and he <= z1]

        strand_char = '+' if best_hit[4] >= 0 else '-'
        abs_left  = win_start + z0
        abs_right = win_start + z1
        source = 'Cicero' if in_cicero else 'proximity'
        print(f'  PASS [{source}] {region_key}  motif={best_hit[2]}  '
              f'score={best_hit[3]:.2f}  strand={strand_char}', flush=True)

        try:
            genome_obj = getattr(scp.genome, genome)
            raw_seq = genome_obj.fetch_seq(chrom, abs_left, abs_right)
            ref_seq = _revcomp(raw_seq) if strand_char == '-' else raw_seq
        except Exception as e:
            print(f'  [WARN] fetch_seq failed for {region_key}: {e}', flush=True)
            ref_seq = 'N' * zoom_w

        passed.append({
            'region':      region_key,
            'chrom':       chrom,
            'win_start':   win_start,
            'z0': z0, 'z1': z1,
            'zoom_w':      zoom_w,
            'abs_left':    abs_left,
            'abs_right':   abs_right,
            'msfp_zoom':   msfp_zoom,
            'scales_plot': scales_plot,
            'zoom_hits':   zoom_hits,
            'best_hit':    best_hit,
            'strand':      strand_char,
            'ref_seq':     ref_seq,
            'source':      source,
        })

        if len(passed) >= 20:
            print('  [INFO] reached 20 enhancer limit', flush=True)
            break

    return passed


# ── Figure renderer ───────────────────────────────────────────────────────────

def render_enhancer_strip(region_data_list, args, shared_vmax, tf_label, ct_label,
                          mode='absolute'):
    N = len(region_data_list)
    tfs_upper = [t.strip().upper() for t in args.tfs]

    row_h_msfp    = 2.5
    row_h_logo    = 0.55
    row_h_seq     = 0.32
    top_margin    = 0.9
    bottom_margin = 0.6
    fig_h = top_margin + bottom_margin + N * (row_h_msfp + row_h_logo + row_h_seq)
    fig_w = 7.0

    left_col  = 0.14
    right_col = 0.88
    cbar_x0   = 0.91

    top_frac    = 1.0 - top_margin    / fig_h
    bottom_frac =       bottom_margin / fig_h

    hr = []
    for _ in range(N):
        hr.extend([row_h_msfp, row_h_logo, row_h_seq])

    gs = gridspec.GridSpec(
        N * 3, 1, height_ratios=hr, hspace=0.04,
        left=left_col, right=right_col,
        top=top_frac, bottom=bottom_frac,
    )

    fig    = plt.figure(figsize=(fig_w, fig_h))
    ax_ref = None
    ax_m_list, ax_l_list, ax_s_list = [], [], []

    for k in range(N):
        ax_m = fig.add_subplot(gs[3 * k], sharex=ax_ref)
        ax_l = fig.add_subplot(gs[3 * k + 1],
                               sharex=ax_ref if ax_ref is not None else ax_m)
        ax_s = fig.add_subplot(gs[3 * k + 2],
                               sharex=ax_ref if ax_ref is not None else ax_m)
        if ax_ref is None:
            ax_ref = ax_m
        ax_m_list.append(ax_m)
        ax_l_list.append(ax_l)
        ax_s_list.append(ax_s)

    sm = plt.cm.ScalarMappable(
        cmap=plt.get_cmap('RdBu_r'),
        norm=mcolors.Normalize(vmin=-shared_vmax, vmax=shared_vmax),
    )
    sm.set_array([])

    for k, d in enumerate(region_data_list):
        ax_m   = ax_m_list[k]
        ax_l   = ax_l_list[k]
        zoom_w = d['zoom_w']
        scales = d['scales_plot']

        ax_m.imshow(
            d['msfp_zoom'], aspect='auto', origin='lower',
            cmap='RdBu_r', vmin=-shared_vmax, vmax=shared_vmax,
            interpolation='bilinear',
            extent=[0, zoom_w, scales[0], scales[-1]],
        )
        ax_m.set_yticks([s for s in [2, 5, 10, 15, 20, 25, 30] if s <= scales[-1]])
        ax_m.tick_params(axis='y', labelsize=6)
        ax_m.tick_params(axis='x', labelbottom=False)

        region_label = d['region']
        src_tag      = f"[{d['source']}]"
        ax_m.set_ylabel(f"{region_label}\n{src_tag}", fontsize=6, rotation=0,
                        labelpad=90, va='center')

        if k == 0:
            mode_tag = '  [Δ TG − WT]' if mode == 'differential' else ''
            ax_m.set_title(f'Enhancer MSFP  ·  {tf_label}  ·  {ct_label}{mode_tag}',
                           fontsize=9, loc='center', pad=3)
        if k == N - 1:
            ax_m.text(1.01, 0.5, 'scale (bp)', fontsize=7, rotation=90,
                      ha='left', va='center', transform=ax_m.transAxes)

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
                  ha='left', va='bottom', fontsize=6, color='#333333',
                  transform=ax_l.transAxes)
        ax_l.text(0.99, 0.10, f"{chev} {d['abs_right']:,}",
                  ha='right', va='bottom', fontsize=6, color='#333333',
                  transform=ax_l.transAxes)

        ax_s = ax_s_list[k]
        draw_seq_track(ax_s, d.get('ref_seq', ''), zoom_w)
        strand_label = f"ref seq  ({d['strand']} strand)"
        ax_s.text(-0.01, 0.5, strand_label, ha='right', va='center',
                  fontsize=6, color='#555555', transform=ax_s.transAxes)

    cbar_ax = fig.add_axes([cbar_x0, bottom_frac + 0.05,
                            0.025, top_frac - bottom_frac - 0.10])
    cbar_label = 'ΔMSFP (TG − WT)' if mode == 'differential' else 'MSFP score'
    fig.colorbar(sm, cax=cbar_ax, label=cbar_label)
    cbar_ax.tick_params(labelsize=6)

    mode_tag = '  ·  Δ TG − WT' if mode == 'differential' else ''
    fig.suptitle(
        f'Enhancer MSFP  ·  {tf_label}  ·  {ct_label}  ·  '
        f'{N} regions  ·  {args.zoom_half_width * 2} bp window{mode_tag}',
        fontsize=10, y=1.0 - 0.04 / fig_h,
    )
    return fig


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Render enhancer MSFP strip (pipeline-ready version).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required
    p.add_argument('--enhancer-h5ad', required=True,
                   help='enhancer_footprints_{CT}_{TF}.h5ad from '
                        'TARGETED_ENHANCER_FOOTPRINTING or ENHANCER_FOOTPRINTING_PER_CT')
    p.add_argument('--tfs',           required=True,
                   help='Comma-separated TF names to scan (e.g. Spi1,Irf8)')
    p.add_argument('--pfm',           required=True,
                   help='JASPAR .jaspar PFM file (e.g. JASPAR2022_core_nonredundant.jaspar)')
    p.add_argument('--cache-dir',     required=True,
                   help='scPrinter cache directory (contains genome FASTA index)')
    p.add_argument('--out-png',       required=True,
                   help='Output PNG path')
    # Optional: region filtering
    p.add_argument('--target-gene',      default=None,
                   help='Gene name for TSS proximity filter')
    p.add_argument('--context-bp',       type=int, default=500_000,
                   help='TSS proximity window (bp)')
    p.add_argument('--cicero-connections', default=None,
                   help='cicero_connections.tsv.gz for Cicero linkage filter')
    p.add_argument('--gtf',              default=None,
                   help='GTF for TSS lookup and Cicero proximity filter. '
                        'Both filters are skipped if not provided.')
    # Optional: rendering
    p.add_argument('--genome',           default='mm10',
                   help='Genome key for scPrinter (e.g. mm10, hg38)')
    p.add_argument('--max-scale',        type=int, default=30,
                   help='Scale ceiling in bp (TF-binding band)')
    p.add_argument('--zoom-half-width',  type=int, default=75,
                   help='Half-width of zoom window centred on best motif hit')
    p.add_argument('--cell-type',        default='',
                   help='Cell-type label for figure title (falls back to h5ad stem)')
    p.add_argument('--mode',             default='absolute',
                   choices=['absolute', 'differential'],
                   help='absolute: first condition slice; '
                        'differential: arr[1]-arr[0] (treatment minus control)')
    # Optional: output
    p.add_argument('--out-pdf',          default=None,
                   help='PDF output path. Defaults to --out-png with .pdf extension.')
    return p.parse_args()


def main():
    args = parse_args()
    tfs      = tuple(t.strip() for t in args.tfs.split(',') if t.strip())
    args.tfs = tfs

    # Resolve gene TSS (only when --gtf is provided)
    gene_chrom, gene_tss = None, None
    if args.target_gene and args.gtf and Path(args.gtf).exists():
        gene_chrom, gene_tss = get_gene_tss(args.gtf, args.target_gene)
        if gene_tss:
            print(f'[msfp_enh] target gene {args.target_gene}: '
                  f'{gene_chrom}:{gene_tss}', flush=True)
        else:
            print(f'[msfp_enh][WARN] {args.target_gene} not found in GTF', flush=True)
    elif args.target_gene and not args.gtf:
        print('[msfp_enh][WARN] --target-gene given but --gtf not provided; '
              'TSS proximity filter disabled.', flush=True)

    # Cicero-linked peaks (requires both --cicero-connections and --gtf)
    cicero_linked = set()
    if args.cicero_connections and args.target_gene and args.gtf:
        cicero_linked = get_cicero_peaks_for_gene(
            args.cicero_connections, args.target_gene, args.gtf)
        print(f'[msfp_enh] {len(cicero_linked)} Cicero-linked peaks for '
              f'{args.target_gene}', flush=True)

    region_data = load_enhancer_regions(
        h5ad_path=args.enhancer_h5ad,
        tfs=tfs,
        pfm_path=args.pfm,
        cache_dir=args.cache_dir,
        genome=args.genome,
        max_scale=args.max_scale,
        zoom_half_width=args.zoom_half_width,
        gene_chrom=gene_chrom,
        gene_tss=gene_tss,
        context_bp=args.context_bp,
        cicero_linked_peaks=cicero_linked,
        mode=args.mode,
    )

    if not region_data:
        print('[msfp_enh] WARN: no enhancer regions passed filters — skipping', flush=True)
        sys.exit(77)

    shared_vmax = max(float(np.percentile(np.abs(d['msfp_zoom']), 99))
                      for d in region_data)
    shared_vmax = max(shared_vmax, 1e-6)
    print(f'[msfp_enh] {len(region_data)} regions, shared_vmax={shared_vmax:.4f}',
          flush=True)

    tf_label = ','.join(tfs)
    ct_label = args.cell_type or Path(args.enhancer_h5ad).stem

    fig = render_enhancer_strip(region_data, args, shared_vmax, tf_label, ct_label,
                                mode=args.mode)

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    print(f'[msfp_enh] wrote {out_png}', flush=True)

    # Always write PDF — derive path from PNG if --out-pdf not given.
    out_pdf = Path(args.out_pdf) if args.out_pdf else out_png.with_suffix('.pdf')
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'[msfp_enh] wrote {out_pdf}', flush=True)

    plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
