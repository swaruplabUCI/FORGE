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


def find_footprint_pngs(footprints_dir, tf_name):
    """Find pre-computed footprint PNGs from ENHANCER_FOOTPRINTING for this TF.

    FIX-97: Uses pre-computed PNGs instead of loading 1.9GB printer h5ad.
    Returns dict of {cell_type: composite_png_path}.
    """
    import glob
    results = {}
    if not footprints_dir or not os.path.isdir(footprints_dir):
        return results

    # ENHANCER_FOOTPRINTING outputs: {cell_type}/{TF}/enhancer_plots/{cell_type}_{TF}_composite.png
    pattern = os.path.join(footprints_dir, f'**/*{tf_name}*_composite.png')
    for png in glob.glob(pattern, recursive=True):
        basename = os.path.basename(png)
        # Extract cell type from filename: {cell_type}_{TF}_composite.png
        ct = basename.replace(f'_{tf_name}_composite.png', '')
        if ct and ct != basename:
            results[ct] = png

    if results:
        print(f"[D4] Found {len(results)} pre-computed footprint composites for {tf_name}: {list(results.keys())}")
    return results


def render_motif_logo(motif_pwm, tf_name, ax=None):
    """Render a sequence logo from a position weight matrix."""
    try:
        import logomaker
    except ImportError:
        print("[WARN] logomaker not available — skipping motif logo")
        return ax

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 1.5))

    logo = logomaker.Logo(
        motif_pwm,
        ax=ax,
        color_scheme='classic',
        font_name='DejaVu Sans',
    )
    ax.set_title(tf_name, fontsize=12, fontweight='bold')
    ax.set_ylabel('bits')
    ax.set_xlim([-0.5, len(motif_pwm) - 0.5])
    return ax


def load_motif_pwm(tf_name, motif_scan_path):
    """Attempt to load a motif PWM for the given TF from the motif scan output."""
    # Try to extract from scprinter motif database
    try:
        import scprinter as scp
        # scprinter stores motifs internally; try to access them
        pwm = scp.motifs.get_motif_pwm(tf_name)
        if pwm is not None:
            return pd.DataFrame(pwm, columns=['A', 'C', 'G', 'T'])
    except Exception:
        pass

    # Fallback: generate a placeholder PWM
    print(f"[WARN] Could not load PWM for {tf_name} — using placeholder")
    return None


def assemble_composite_figure(
    browser_png, footprint_pngs, motif_pwm, tf_name, gene,
    title='', save_path=None, dpi=200,
):
    """FIX-97: Assemble composite from browser panel + pre-computed footprint PNGs.

    Args:
        browser_png: Path to pyGenomeTracks browser panel PNG (or None)
        footprint_pngs: Dict of {cell_type: composite_png_path} from ENHANCER_FOOTPRINTING
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

    # Panels B+: Pre-computed footprint composites per cell type
    for i, (ct, png_path) in enumerate(sorted(footprint_pngs.items())):
        if os.path.exists(png_path):
            panel_imgs.append(mpimg.imread(png_path))
            letter = chr(ord('B') + i)
            panel_labels.append(f'{letter}. {ct} — {tf_name} Enhancer Footprint')

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
    fig_width = 24
    MIN_PANEL_HEIGHT = 3.0
    MAX_PANEL_HEIGHT = 6.0
    panel_heights = []
    for img in panel_imgs:
        h, w = img.shape[:2]
        natural = fig_width * (h / w)
        panel_heights.append(max(MIN_PANEL_HEIGHT, min(MAX_PANEL_HEIGHT, natural)))
    motif_height = 2.0
    total_height = sum(panel_heights) + (motif_height if has_motif else 0)
    height_ratios = panel_heights + ([motif_height] if has_motif else [])
    n_rows = n_panels + (1 if has_motif else 0)

    fig = plt.figure(figsize=(fig_width, total_height))
    gs = gridspec.GridSpec(n_rows, 1, figure=fig, height_ratios=height_ratios, hspace=0.3)

    for row, (img, label) in enumerate(zip(panel_imgs, panel_labels)):
        ax = fig.add_subplot(gs[row])
        ax.imshow(img, aspect='auto')
        ax.axis('off')
        ax.set_title(label, fontsize=14, fontweight='bold', loc='left')

    # Motif logo at bottom
    if has_motif:
        ax_logo = fig.add_subplot(gs[n_panels])
        render_motif_logo(motif_pwm, tf_name, ax=ax_logo)
        letter = chr(ord('A') + n_panels)
        ax_logo.set_title(f'{letter}. TF Motif: {tf_name}',
                         fontsize=14, fontweight='bold', loc='left')

    if title:
        fig.suptitle(title, fontsize=18, fontweight='bold', y=1.01)

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
    parser.add_argument('--gene', required=True, help='Gene name to render')
    parser.add_argument('--tf-name', required=True, help='TF name')
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
            'skipped': True,
            'reason': 'gene_not_in_track_manifest',
        }
        with open(os.path.join(args.outdir, f'summary_{gene_original}_{tf_name}.json'), 'w') as f:
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

    print(f"[D3-D7] Rendering composite for {gene} / {tf_name}")
    print(f"        Region: {region_str}")

    # Step D3: Render pyGenomeTracks browser panel
    browser_png = None
    if region_str:
        browser_png = os.path.join(args.outdir, f'browser_{gene_original}.png')
        result = render_pygenometracks(ini_path, region_str, browser_png, dpi=args.dpi)
        if result is None:
            browser_png = None

    # Step D4 (FIX-97): Find pre-computed footprint PNGs from ENHANCER_FOOTPRINTING
    footprint_pngs = find_footprint_pngs(args.footprints_dir, tf_name)

    # Step D5: Load TF motif PWM
    motif_pwm = load_motif_pwm(tf_name, args.motif_scan)

    # Step D6: Assemble composite figure
    composite_path = os.path.join(args.outdir, f'composite_{gene_original}_{tf_name}.png')
    fig = assemble_composite_figure(
        browser_png=browser_png,
        footprint_pngs=footprint_pngs,
        motif_pwm=motif_pwm,
        tf_name=tf_name,
        gene=gene,
        title=f'Regulatory Architecture: {gene} — {tf_name}',
        save_path=composite_path,
        dpi=args.dpi,
    )

    if fig:
        plt.close(fig)

    # Write summary entry
    summary = {
        'gene': gene,
        'tf': tf_name,
        'region': region_str,
        'composite_png': composite_path if fig else None,
        'browser_png': browser_png,
        'footprint_cell_types': list(footprint_pngs.keys()),
        'has_motif': motif_pwm is not None,
    }
    summary_path = os.path.join(args.outdir, f'summary_{gene_original}_{tf_name}.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"[composite_enhancer_viz] Done — {gene} / {tf_name}")


if __name__ == '__main__':
    main()
