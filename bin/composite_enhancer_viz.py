#!/usr/bin/env python3
"""
composite_enhancer_viz.py — Recipe D Steps D3-D7

FIX-97: Rearchitected to remove scPRINTER on-the-fly panel generation (which hung
loading 1.9GB printer h5ad and never worked). Now assembles from:
  Panel A: Genome browser tracks (pyGenomeTracks — gene models, CCANs, ATAC coverage)
  Panel B: Pre-computed enhancer footprint composites from ENHANCER_FOOTPRINTING
  Panel C: TF motif logo (logomaker)

Inputs:
  --track-manifest       JSON from prepare_enhancer_viz_tracks.py
  --motif-scan           Enhancer motif scan TSV (from MOTIF_SCAN_ENHANCERS)
  --footprints-dir       Directory with enhancer footprint PNGs (from ENHANCER_FOOTPRINTING)
  --gene                 Gene name to render
  --tf-name              TF name
  --dpi                  Output DPI (default 200)
  --outdir               Output directory (default .)
"""
import argparse
import json
import os
import re
import sys
import subprocess
import warnings

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import numpy as np
import pandas as pd


def render_pygenometracks(config_path, region_str, output_png, dpi=150, width=40):
    """Render genome browser tracks for a given region via pyGenomeTracks CLI."""
    cmd = [
        'pyGenomeTracks',
        '--tracks', config_path,
        '--region', region_str,
        '--outFileName', output_png,
        '--dpi', str(dpi),
        '--width', str(width),
        '--fontSize', '10',
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"[D3] Browser tracks rendered to {output_png}")
        return output_png
    except subprocess.CalledProcessError as e:
        print(f"[WARN] pyGenomeTracks failed: {e.stderr[:500]}")
        return None
    except FileNotFoundError:
        print("[WARN] pyGenomeTracks not found — skipping browser panel")
        return None


def find_promoter_pngs(promoter_dir, gene, cell_type=None):
    """Find pre-computed promoter MSFP composite PNGs from SCPRINTER_FOOTPRINTING.

    Restored after FIX-97: the original promoter panel loaded the 1.9GB printer
    h5ad on-the-fly and hung. The pre-rendered PNGs already exist at
    `{promoter_dir}/{ct_safe}/plots/{ct}_{gene}_composite.png` so we can stitch
    them in without any heavy load.

    Pattern (per scprinter/footprints publishDir layout):
        {promoter_dir}/{ct_safe}/plots/{ct_label}_{gene}_composite.png

    Returns dict {{cell_type: composite_png_path}}, restricted to one ct if
    `cell_type` is provided.
    """
    import glob
    results = {}
    if not promoter_dir or not os.path.isdir(promoter_dir):
        return results

    # The composite_png contains all three sub-panels (insertion / footprint
    # / TFBS) — same pre-rendering convention as the enhancer outputs.
    pattern = os.path.join(promoter_dir, f'**/plots/*_{gene}_composite.png')
    for png in glob.glob(pattern, recursive=True):
        basename = os.path.basename(png)
        suffix = f'_{gene}_composite.png'
        if not basename.endswith(suffix):
            continue
        ct = basename[:-len(suffix)]
        if not ct:
            continue
        if cell_type and ct != cell_type:
            continue
        results[ct] = png

    if results:
        print(f"[D4b] Found {len(results)} pre-computed promoter composites for "
              f"{gene}{' / ' + cell_type if cell_type else ''}: {list(results.keys())}")
    return results


def find_footprint_pngs(footprints_dir, tf_name, cell_type=None):
    """Find pre-computed footprint PNGs from ENHANCER_FOOTPRINTING for this TF.

    FIX-97: Uses pre-computed PNGs instead of loading 1.9GB printer h5ad.
    Returns dict of {cell_type: composite_png_path}.

    If cell_type is given, restricts the result to that single cell type so the
    per-(ct) split COMPOSITE only assembles its own row.
    """
    import glob
    results = {}
    if not footprints_dir or not os.path.isdir(footprints_dir):
        return results

    # ENHANCER_FOOTPRINTING outputs: {cell_type}/{TF}/enhancer_plots/{cell_type}_{TF}_composite.png
    # Substring match on tf_name; resolve correct ct via prefix-strip.
    pattern = os.path.join(footprints_dir, f'**/*{tf_name}*_composite.png')
    for png in glob.glob(pattern, recursive=True):
        basename = os.path.basename(png)
        suffix = f'_{tf_name}_composite.png'
        if not basename.endswith(suffix):
            continue
        ct = basename[:-len(suffix)]
        if not ct:
            continue
        if cell_type and ct != cell_type:
            continue
        results[ct] = png

    if results:
        print(f"[D4] Found {len(results)} pre-computed footprint composites for "
              f"{tf_name}{' / ' + cell_type if cell_type else ''}: {list(results.keys())}")
    return results


def filter_ini_to_celltype(ini_path, target_ct, out_path):
    """Rewrite a multi-track .ini to retain only the fixed sections + one ct's bigwig.

    The ini built by prepare_enhancer_viz_tracks.py has 4 fixed sections
    (`[genes]`, `[spacer]`, `[links]`, `[enhancers]`), one section per cell-type
    bigwig (with `title = <raw cell_type>`), and `[x-axis]`. Per-ct COMPOSITE
    only needs one bigwig — this slims the ini accordingly so pyGenomeTracks
    renders a single coverage track instead of 31, dropping runtime ~30x.

    Returns the path written, or None if the target_ct section was not found.
    """
    keep_fixed = {'genes', 'spacer', 'links', 'enhancers', 'x-axis'}
    sections = []
    cur_name = None
    cur_lines = []
    with open(ini_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('[') and stripped.endswith(']'):
                if cur_name is not None:
                    sections.append((cur_name, cur_lines))
                cur_name = stripped[1:-1]
                cur_lines = [line]
            else:
                if cur_name is None:
                    continue
                cur_lines.append(line)
    if cur_name is not None:
        sections.append((cur_name, cur_lines))

    target_ct_lower = target_ct.lower()
    kept = []
    found_target = False
    for name, body in sections:
        if name in keep_fixed:
            kept.append((name, body))
            continue
        title_match = False
        for ln in body:
            stripped = ln.strip()
            if stripped.startswith('title'):
                _, _, val = stripped.partition('=')
                if val.strip().lower() == target_ct_lower:
                    title_match = True
                    break
        if title_match:
            kept.append((name, body))
            found_target = True

    if not found_target:
        return None

    with open(out_path, 'w') as f:
        for _, body in kept:
            for ln in body:
                f.write(ln)
            if not body[-1].endswith('\n'):
                f.write('\n')
    return out_path


def render_motif_logo(motif_pwm, tf_name, ax=None):
    """Render a sequence logo from a position frequency matrix.

    `motif_pwm` is a DataFrame with columns {A,C,G,T} where each row sums to
    1 (probabilities), as produced by load_motif_pwm. logomaker's
    transform_matrix is applied to convert probabilities → information
    content (bits) so the logo height encodes per-position conservation.
    """
    try:
        import logomaker
    except ImportError:
        print("[WARN] logomaker not available — skipping motif logo")
        return ax

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 1.5))

    # Convert probability matrix → information matrix (bits per base)
    try:
        info_matrix = logomaker.transform_matrix(
            motif_pwm,
            from_type='probability',
            to_type='information',
        )
    except Exception as e:
        # If the matrix is already counts (legacy callers), let logomaker
        # auto-handle it.
        print(f"[D5] transform_matrix fallback ({e}); rendering raw matrix")
        info_matrix = motif_pwm

    logomaker.Logo(
        info_matrix,
        ax=ax,
        color_scheme='classic',
        font_name='DejaVu Sans',
    )
    ax.set_title(tf_name, fontsize=12, fontweight='bold')
    ax.set_ylabel('bits')
    ax.set_xlim([-0.5, len(info_matrix) - 0.5])
    return ax


_JASPAR_PFM_CACHE = {}


def _parse_jaspar_pfm(pfm_path):
    """Parse a JASPAR .jaspar PFM file into {tf_name_upper: (motif_id, df_freq)}.

    JASPAR2022 format (one record per motif):
        >MA0867.2\tSOX4
        A  [   12   34   ...   ]
        C  [    1    9   ...   ]
        G  [    8   17   ...   ]
        T  [   23   11   ...   ]

    For each motif, counts are normalized column-wise to frequencies
    (each position sums to 1). The dict is keyed on the upper-case TF name
    (case-insensitive lookup downstream) and additionally on the motif ID
    (e.g. 'MA0867.2'). When a TF appears in multiple motif versions, the
    LATEST record wins (file order = ascending MA-ID, so this picks the
    newest version naturally).
    """
    if pfm_path in _JASPAR_PFM_CACHE:
        return _JASPAR_PFM_CACHE[pfm_path]

    out = {}
    cur_id = None
    cur_tf = None
    cur_rows = {}
    with open(pfm_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_tf and cur_rows.keys() == {'A', 'C', 'G', 'T'}:
                    df = pd.DataFrame(cur_rows)
                    col_sums = df.sum(axis=1).replace(0, np.nan)
                    freq = df.div(col_sums, axis=0).fillna(0.25)
                    out[cur_tf.upper()] = (cur_id, freq)
                    out[cur_id] = (cur_id, freq)
                # Header parse
                parts = line[1:].split('\t')
                cur_id = parts[0].strip()
                cur_tf = parts[1].strip() if len(parts) > 1 else cur_id
                cur_rows = {}
                continue
            # Body row: "A  [ 12 34 ... ]"
            stripped = line.strip()
            if not stripped:
                continue
            for base in ('A', 'C', 'G', 'T'):
                if stripped.startswith(base):
                    inner = stripped[len(base):].strip()
                    if inner.startswith('['):
                        inner = inner[1:]
                    if inner.endswith(']'):
                        inner = inner[:-1]
                    try:
                        vals = [float(x) for x in inner.split()]
                    except ValueError:
                        vals = []
                    cur_rows[base] = vals
                    break
        # Flush trailing record
        if cur_tf and cur_rows.keys() == {'A', 'C', 'G', 'T'}:
            df = pd.DataFrame(cur_rows)
            col_sums = df.sum(axis=1).replace(0, np.nan)
            freq = df.div(col_sums, axis=0).fillna(0.25)
            out[cur_tf.upper()] = (cur_id, freq)
            out[cur_id] = (cur_id, freq)

    _JASPAR_PFM_CACHE[pfm_path] = out
    return out


def load_motif_pwm(tf_name, motif_scan_path, jaspar_pfms_path=''):
    """Load a motif PWM for `tf_name` from the JASPAR PFM file used upstream.

    Strategy:
      1. If `jaspar_pfms_path` is given (the canonical `params.scprinter.pfms`
         file), parse it and look up by case-insensitive TF name.
      2. Otherwise, fall back to scprinter's bundled JASPAR2022_core() (returns
         a path in the user's scprinter cache; same md5 as the canonical file
         on this system).
      3. If neither lookup succeeds, return None (logo panel is dropped).

    Returns a frequency DataFrame[positions x {A,C,G,T}] suitable for
    logomaker, or None on failure. logomaker auto-transforms frequencies
    into information-content (bits) before rendering, so callers can label
    the y-axis 'bits' as expected.
    """
    pfms_path = jaspar_pfms_path
    if not pfms_path or not os.path.exists(pfms_path):
        try:
            import scprinter as scp
            pfms_path = scp.motifs.JASPAR2022_core()
        except Exception as e:
            print(f"[WARN] Could not resolve JASPAR PFMs ({e}); skipping logo")
            return None

    if not pfms_path or not os.path.exists(pfms_path):
        print(f"[WARN] JASPAR PFM file not found at {pfms_path!r}; skipping logo")
        return None

    try:
        pfm_index = _parse_jaspar_pfm(pfms_path)
    except Exception as e:
        print(f"[WARN] Failed to parse JASPAR PFMs at {pfms_path}: {e}")
        return None

    # Case-insensitive TF lookup, then motif-ID lookup
    hit = pfm_index.get(tf_name.upper()) or pfm_index.get(tf_name)
    if hit is None:
        print(f"[WARN] PWM not found for TF '{tf_name}' in JASPAR PFMs "
              f"({len([k for k in pfm_index if not k.startswith('MA')])} TF names indexed)")
        return None

    motif_id, freq_df = hit
    print(f"[D5] Loaded PWM for '{tf_name}' (JASPAR {motif_id}, {len(freq_df)} positions)")
    return freq_df


def _crop_inner_whitespace(img, top_search_frac=0.30, blank_thresh=2.0,
                           keep_pad_px=24):
    """ITER-v4: collapse the long blank-row band baked into the source
    enhancer-footprint PNG between its header lines and the first sub-panel.

    The source PNG (from ENHANCER_FOOTPRINTING) is ~9500 px tall with a few
    text-header rows at top, then ~400 px of pure white, then the first
    heatmap. matplotlib's imshow can't address that gap because it's
    embedded in the pixels themselves. We scan for the longest run of
    near-white rows in the upper `top_search_frac` of the image and excise
    it (keeping `keep_pad_px` rows of margin on each side so the header
    isn't pressed against the heatmap).

    Args:
        img: ndarray HxWx{3,4} as returned by mpimg.imread (float in [0,1]).
        top_search_frac: only look in the top fraction of the image so we
            don't accidentally collapse legitimate gaps lower down.
        blank_thresh: per-row darkness threshold (in 0-255 scale) below
            which a row is treated as blank. Default 2.0 = ~99.2% white.
        keep_pad_px: rows of blank margin to preserve at each side of the
            removed band.

    Returns the (possibly cropped) image. Falls through to original image
    on any failure to keep the pipeline robust.
    """
    try:
        if img.ndim < 2:
            return img
        h, w = img.shape[:2]
        if h < 200:
            return img

        rgb = img[..., :3] if img.shape[-1] >= 3 else img
        scale = 255.0 if rgb.dtype.kind == 'f' else 1.0
        # Per-row darkness = 255 - mean luminance. Larger = more content.
        row_dark = 255.0 - rgb.mean(axis=(1, 2)) * scale

        upper = int(h * top_search_frac)
        is_blank = row_dark[:upper] < blank_thresh

        # Find longest run of blank rows in upper band.
        best_start, best_len = 0, 0
        cur_start, cur_len = 0, 0
        for i, b in enumerate(is_blank):
            if b:
                if cur_len == 0:
                    cur_start = i
                cur_len += 1
                if cur_len > best_len:
                    best_start, best_len = cur_start, cur_len
            else:
                cur_len = 0

        # Only excise if the gap is meaningful (>~0.5" at 200 dpi).
        if best_len < 100:
            return img

        cut_start = best_start + keep_pad_px
        cut_end = best_start + best_len - keep_pad_px
        if cut_end <= cut_start:
            return img

        kept = np.concatenate([img[:cut_start], img[cut_end:]], axis=0)
        return kept
    except Exception as e:
        print(f"[ITER-v4] _crop_inner_whitespace skipped ({e}); using raw image")
        return img


def assemble_composite_figure(
    browser_png, footprint_pngs, motif_pwm, tf_name, gene,
    title='', save_path=None, dpi=200, promoter_pngs=None,
):
    """FIX-97: Assemble composite from browser panel + pre-computed footprint PNGs.

    Args:
        browser_png: Path to pyGenomeTracks browser panel PNG (or None)
        footprint_pngs: Dict of {cell_type: composite_png_path} from ENHANCER_FOOTPRINTING
        promoter_pngs: Dict of {cell_type: composite_png_path} from
            SCPRINTER_FOOTPRINTING (per-(ct, gene) promoter MSFP). May be None
            or empty if scPRINTER didn't run for this gene — panel just gets
            skipped. Restores the promoter panel removed by FIX-97 using
            pre-rendered PNGs (no h5ad load → no hang).
        motif_pwm: DataFrame PWM for motif logo (or None)
        tf_name: TF name for labels
        gene: Gene name for labels
    """
    panel_imgs = []
    panel_labels = []

    # Panel A: Browser tracks
    if browser_png and os.path.exists(browser_png):
        panel_imgs.append(mpimg.imread(browser_png))
        panel_labels.append('A. Genomic Context: Gene Model, CCANs, ATAC Coverage')

    # Panel B (Tier-1 promoter binding): pre-rendered per-(ct, gene) MSFP
    # composite. Inserted before enhancer-footprint panels so the visual
    # hierarchy mirrors the filter cascade (promoter > enhancer-TF >
    # differential). Skipped silently when promoter_pngs is empty.
    if promoter_pngs:
        for ct, png_path in sorted(promoter_pngs.items()):
            if os.path.exists(png_path):
                panel_imgs.append(mpimg.imread(png_path))
                letter = chr(ord('A') + len(panel_imgs) - 1)
                panel_labels.append(f'{letter}. {ct} — {gene} Promoter MSFP (any-TF binding)')

    # Next panels: pre-computed enhancer-TF footprint composites
    # ITER-v4: collapse the inner whitespace gap baked into the source PNG
    # (between its header lines and the first sub-panel).
    base_offset = len(panel_imgs)
    for i, (ct, png_path) in enumerate(sorted(footprint_pngs.items())):
        if os.path.exists(png_path):
            img = _crop_inner_whitespace(mpimg.imread(png_path))
            panel_imgs.append(img)
            # ITER-v4: source PNG carries its own internal heading
            # ("OPC-Oligo - TF SOX10 Enhancer Footprinting | ..."), so the
            # matplotlib panel title was redundant. Suppress with empty label.
            panel_labels.append('')

    if not panel_imgs:
        print("[WARN] No panels available — cannot assemble composite figure")
        return None

    # Determine layout. Floor + cap each panel height so multi-panel stacks
    # (e.g. browser + 5 footprint cell types) do not squeeze each panel to
    # ~3 in, and single-panel browser-only cases do not blow up to ~78 in.
    # aspect='auto' lets each imshow fill its row cleanly at the allocated
    # height; source PNGs carry their own labels so mild aspect distortion
    # is preferable to letterboxing or cramming.
    n_panels = len(panel_imgs)
    has_motif = motif_pwm is not None
    # ITER-v1: vertical expansion. fig_width 24 -> 16 (less squat) and
    # MAX_PANEL_HEIGHT 6 -> 60 so the 9527px-tall enhancer footprint
    # composite (aspect 0.21) renders close to its natural height.
    fig_width = 16
    MIN_PANEL_HEIGHT = 3.0
    MAX_PANEL_HEIGHT = 60.0
    panel_heights = []
    for img in panel_imgs:
        h, w = img.shape[:2]
        natural = fig_width * (h / w)
        panel_heights.append(max(MIN_PANEL_HEIGHT, min(MAX_PANEL_HEIGHT, natural)))
    motif_height = 2.0
    total_height = sum(panel_heights) + (motif_height if has_motif else 0)
    height_ratios = panel_heights + ([motif_height] if has_motif else [])
    n_rows = n_panels + (1 if has_motif else 0)

    # ITER-v3: tighten whitespace.
    #   - hspace 0.3 -> 0.02 (kill large gridspec inter-panel gap)
    #   - set_title pad 6 (default) -> 2 (tighten label-to-content)
    #   - suptitle y 1.01 -> 0.999 + GridSpec top=0.985 to hug header
    fig = plt.figure(figsize=(fig_width, total_height))
    gs = gridspec.GridSpec(n_rows, 1, figure=fig, height_ratios=height_ratios,
                            hspace=0.02, top=0.985, bottom=0.005)

    for row, (img, label) in enumerate(zip(panel_imgs, panel_labels)):
        ax = fig.add_subplot(gs[row])
        ax.imshow(img, aspect='auto')
        ax.axis('off')
        # ITER-v4: skip set_title when label is empty (source PNG already
        # carries an internal heading; matplotlib title would duplicate it).
        if label:
            ax.set_title(label, fontsize=14, fontweight='bold', loc='left', pad=2)

    # Motif logo at bottom
    if has_motif:
        ax_logo = fig.add_subplot(gs[n_panels])
        render_motif_logo(motif_pwm, tf_name, ax=ax_logo)
        # ITER-v4: drop the letter prefix — source PNG already labels its
        # internal sub-panels A..J, so adding "B. TF Motif" collides
        # confusingly with that scheme.
        ax_logo.set_title(f'TF Motif: {tf_name}',
                         fontsize=14, fontweight='bold', loc='left', pad=2)

    if title:
        fig.suptitle(title, fontsize=18, fontweight='bold', y=0.999)

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        pdf_path = save_path.replace('.png', '.pdf')
        fig.savefig(pdf_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        print(f"[D6] Composite figure saved to {save_path} and {pdf_path}")

    return fig


def main():
    parser = argparse.ArgumentParser(description='Composite enhancer visualization (Recipe D — FIX-97)')
    parser.add_argument('--track-manifest', required=True, help='JSON from prepare_enhancer_viz_tracks.py')
    parser.add_argument('--motif-scan', default='', help='Enhancer motif scan TSV')
    parser.add_argument('--footprints-dir', default='', help='Directory with pre-computed footprint PNGs')
    parser.add_argument('--promoter-dir', default='',
                        help='Directory with pre-computed scPRINTER per-(ct, gene) '
                             'promoter MSFP composites '
                             '(results/scprinter/footprints). Restores the '
                             'promoter panel removed by FIX-97; only renders '
                             'when this gene was in scprinter.target_genes.')
    parser.add_argument('--gene', required=True, help='Gene name to render')
    parser.add_argument('--tf-name', required=True, help='TF name')
    parser.add_argument('--cell-type', default='',
                        help='If set, render a per-cell-type composite — slim the .ini to this '
                             'ct\'s bigwig only and keep just that ct\'s footprint composite. '
                             'When empty, falls back to legacy multi-ct stack (all bigwigs + '
                             'all footprint cell types).')
    parser.add_argument('--gene-safe', default='', help='Filesystem-safe gene name (FIX-49). Defaults to raw.')
    parser.add_argument('--tf-name-safe', default='', help='Filesystem-safe TF name (FIX-49). Defaults to raw.')
    parser.add_argument('--cell-type-safe', default='',
                        help='Filesystem-safe cell-type name. Defaults to sanitized --cell-type.')
    parser.add_argument('--jaspar-pfms', default='',
                        help='Path to the JASPAR .jaspar PFM file used upstream by '
                             'MOTIF_SCAN_ENHANCERS (e.g. params.scprinter.pfms). When '
                             'empty, falls back to scprinter.motifs.JASPAR2022_core() '
                             'cache. Required for the TF motif logo panel.')
    parser.add_argument('--dpi', type=int, default=200)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load track manifest
    with open(args.track_manifest) as f:
        manifest = json.load(f)

    gene = args.gene
    gene_original = gene  # FIX-95: Preserve original case for output filenames (Nextflow expects this)
    tf_name = args.tf_name
    cell_type = args.cell_type or ''
    # FIX-49: filesystem-safe variants for path construction. Keep raw gene_original/
    # tf_name for manifest lookups, motif PWM matching, JSON content, plot titles.
    gene_fs = args.gene_safe or gene_original
    tf_fs   = getattr(args, 'tf_name_safe', '') or tf_name
    ct_fs   = args.cell_type_safe or (re.sub(r'[/\s\(\)]+', '_', cell_type) if cell_type else '')

    # Suffix used in output filenames — adds _{ct} only in per-ct mode.
    out_suffix = f'{gene_fs}_{tf_fs}' + (f'_{ct_fs}' if ct_fs else '')

    gene_regions = manifest.get('gene_regions', {})
    gene_info = gene_regions.get(gene)
    # Case-insensitive fallback (JASPAR names may differ from GTF case)
    if not gene_info:
        gene_upper_map = {k.upper(): k for k in gene_regions}
        matched = gene_upper_map.get(gene.upper())
        if matched:
            print(f"[D3] Normalized gene '{gene}' -> '{matched}' (manifest case)")
            gene = matched
            gene_info = gene_regions[gene]
    if not gene_info:
        print(f"[SKIP] Gene '{gene}' not found in track manifest (likely JASPAR "
              f"composite/ortholog name with no human gene match). Available: "
              f"{list(gene_regions.keys())[:10]}")
        # Write stub summary: Nextflow declares `summary` non-optional, so a
        # graceful skip with no file is treated as a missing output and fails the run.
        stub = {
            'gene': gene_original,
            'tf': tf_name,
            'cell_type': cell_type or None,
            'skipped': True,
            'reason': 'gene_not_in_track_manifest',
        }
        with open(os.path.join(args.outdir, f'summary_{out_suffix}.json'), 'w') as f:
            json.dump(stub, f, indent=2)
        sys.exit(0)

    ini_path = gene_info['ini']
    # FIX-95: .ini files staged into track_inis/ subdir by Nextflow
    if not os.path.exists(ini_path):
        ini_in_subdir = os.path.join('track_inis', os.path.basename(ini_path))
        if os.path.exists(ini_in_subdir):
            ini_path = ini_in_subdir
    region_str = gene_info['region']
    chrom = gene_info.get('chrom', '')
    tss = gene_info.get('tss', 0)

    print(f"[D3-D7] Rendering composite for {gene} / {tf_name}"
          + (f" / {cell_type}" if cell_type else ""))
    print(f"        Region: {region_str}")

    # In per-ct mode, slim the multi-track .ini down to one bigwig section.
    # Keeps `[genes]`, `[spacer]`, `[links]`, `[enhancers]`, the matching ct
    # section, `[x-axis]` — drops the other 30 ATAC bigwig tracks. This is the
    # main runtime win (pyGenomeTracks goes from rendering 31 coverage tracks
    # to 1).
    ini_to_render = ini_path
    if cell_type:
        slim_path = os.path.join(args.outdir, f'tracks_{out_suffix}.ini')
        result = filter_ini_to_celltype(ini_path, cell_type, slim_path)
        if result:
            ini_to_render = result
            print(f"[D3] Slim ini built for ct={cell_type}: {slim_path}")
        else:
            print(f"[WARN] No bigwig section matched cell_type='{cell_type}' "
                  f"in {ini_path}; using full ini (slow path).")

    # Step D3: Render pyGenomeTracks browser panel
    # Two outputs: PNG (used downstream by COMPOSITE if panel A is re-enabled)
    # and a standalone vector PDF that is the PRIMARY consumer-facing browser
    # asset (the composite no longer embeds Panel A — the multi-panel layout
    # was compressing it to ~3-6 in tall, illegible for gene tracks).
    browser_png = None
    browser_pdf = None
    if region_str:
        browser_png = os.path.join(args.outdir, f'browser_{out_suffix}.png')
        result = render_pygenometracks(ini_to_render, region_str, browser_png, dpi=args.dpi)
        if result is None:
            browser_png = None

        browser_pdf = os.path.join(args.outdir, f'browser_{out_suffix}.pdf')
        pdf_result = render_pygenometracks(ini_to_render, region_str, browser_pdf,
                                           dpi=args.dpi, width=40)
        if pdf_result is None:
            browser_pdf = None

    # Step D4 (FIX-97): Find pre-computed footprint PNGs from ENHANCER_FOOTPRINTING
    footprint_pngs = find_footprint_pngs(
        args.footprints_dir, tf_name,
        cell_type=cell_type if cell_type else None,
    )

    # Step D4b (Tier-1 panel restoration, post-FIX-97): per-(ct, gene) scPRINTER
    # promoter MSFP composite. Empty dict when this gene wasn't in
    # scprinter.target_genes — assemble_composite_figure skips the panel cleanly.
    promoter_pngs = find_promoter_pngs(
        args.promoter_dir, gene_original,
        cell_type=cell_type if cell_type else None,
    )

    # Step D5: Load TF motif PWM (from canonical JASPAR PFMs file)
    motif_pwm = load_motif_pwm(tf_name, args.motif_scan, args.jaspar_pfms)

    # Step D6: Assemble composite figure
    # Panel A (browser tracks) is intentionally omitted — it lives as a
    # standalone vector PDF (browser_*.pdf) instead. Including it here forced
    # the multi-panel composite to compress every panel to 3-6 in tall, which
    # made gene tracks illegible. Composite now starts at panel A = promoter
    # MSFP (or first enhancer footprint when scPRINTER didn't run for the
    # gene), giving each remaining panel ~50% more vertical room.
    composite_path = os.path.join(args.outdir, f'composite_{out_suffix}.png')
    title = (f'Regulatory Architecture: {gene} — {tf_name}'
             + (f' [{cell_type}]' if cell_type else ''))
    fig = assemble_composite_figure(
        browser_png=None,
        footprint_pngs=footprint_pngs,
        promoter_pngs=promoter_pngs,
        motif_pwm=motif_pwm,
        tf_name=tf_name,
        gene=gene,
        title=title,
        save_path=composite_path,
        dpi=args.dpi,
    )

    if fig:
        plt.close(fig)

    # Write summary entry
    summary = {
        'gene': gene,
        'tf': tf_name,
        'cell_type': cell_type or None,
        'region': region_str,
        'composite_png': composite_path if fig else None,
        'browser_png': browser_png,
        'browser_pdf': browser_pdf,
        'footprint_cell_types': list(footprint_pngs.keys()),
        'has_motif': motif_pwm is not None,
    }
    summary_path = os.path.join(args.outdir, f'summary_{out_suffix}.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"[composite_enhancer_viz] Done — {gene} / {tf_name}"
          + (f" / {cell_type}" if cell_type else ""))


if __name__ == '__main__':
    main()
