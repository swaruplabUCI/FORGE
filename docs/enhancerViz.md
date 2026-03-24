# Recipe D: Composite Enhancer Visualization — Agent Recipe

## Purpose

Standalone visualization protocol for assembling genome browser tracks, scprinter footprints, CCAN/DORC arcs, and TF motif logos into a single publication-quality composite figure per locus. Designed to render the outputs of Recipes A–C (enhancer footprinting, multiome regulatory inference, CellChat signaling) into an interpretable visual summary.

This recipe assumes that upstream analyses have already been completed and that the following outputs are available on disk or in memory:

- **scprinter printer object** with computed TFBS binding scores and multiscale footprints
- **Cicero CCAN connections** (CSV) and/or **DORC peak-gene associations** (from scprinter)
- **Cell-type pseudobulk ATAC coverage** (as BAMs or bigWig files)
- **SCENIC+ eRegulon region sets** (optional, for enhancer annotations)
- **TF motif PWMs** (from JASPAR, HOCOMOCO, FigR, or scprinter's motif database)

---

## Prerequisites

```bash
pip install pyGenomeTracks logomaker pyBigWig matplotlib pybedtools --break-system-packages
# pyGenomeTracks also requires: pysam, hic2cool (optional for Hi-C)
# On headless HPC: matplotlib.use('Agg') must be set before importing pyplot
```

```
# Required data files:
- GTF gene annotation (e.g., gencode.v38.annotation.gtf.gz)
- Per-cell-type bigWig coverage files (from deeptools bamCoverage or equivalent)
- Cicero connections CSV or DORC significant associations table
- scprinter h5ad with pre-computed binding scores and footprints
- TF motif PWMs as DataFrames or JASPAR/MEME format files
```

---

## Architecture Overview

The composite figure is built from three rendering engines stitched into a single matplotlib figure:

1. **pyGenomeTracks** — genome browser tracks (bigWig coverage, gene models, CCAN/DORC arcs, peak annotations). Renders into matplotlib axes with genomic coordinate alignment.
2. **scprinter plots** — Tn5 insertion profiles, TFBS binding scores, multiscale footprint heatmaps. These are matplotlib-native and can be captured into subplot axes.
3. **logomaker** — TF motif sequence logos. Renders directly into matplotlib Axes objects.

```
┌─────────────────────────────────────────────────────┐
│  Panel A: Gene model track (pyGenomeTracks / GTF)   │
├─────────────────────────────────────────────────────┤
│  Panel B: CCAN / DORC arc links                     │
│           (pyGenomeTracks links track)              │
├─────────────────────────────────────────────────────┤
│  Panel C: Pseudobulk ATAC coverage per cell type    │
│           (pyGenomeTracks bigWig or scprinter Tn5)  │
├─────────────────────────────────────────────────────┤
│  Panel D: TFBS binding score heatmap                │
│           (scprinter plot_binding_score)             │
├─────────────────────────────────────────────────────┤
│  Panel E: Multiscale footprint heatmap              │
│           (scprinter plot_footprints)               │
├──────────────────────┬──────────────────────────────┤
│  Panel F: TF motif   │  Panel G: Enhancer footprint │
│  logo (logomaker)    │  at distal element (scprinter)│
└──────────────────────┴──────────────────────────────┘
```

---

## Step D1: Prepare Track Data Files

All tracks need to be exported as standard genomic file formats so pyGenomeTracks can render them.

```python
import pandas as pd
import numpy as np
import pyBigWig

# --- D1a: Export pseudobulk ATAC coverage as bigWig files ---
# scprinter stores Tn5 insertions internally. For pyGenomeTracks,
# export per-cell-type coverage as bigWig files.
#
# Option 1: Use scprinter's internal data to build bigWigs
# Option 2: Use cellranger-arc pseudobulk BAMs + deeptools bamCoverage
# Option 3: If you already have per-group bigWigs, skip this step.
#
# Using deeptools (bash, most common on HPC):
# for group in celltype1 celltype2 celltype3; do
#     bamCoverage -b ${group}_pseudobulk.bam -o ${group}_coverage.bw \
#         --normalizeUsing RPKM --binSize 10 --smoothLength 30
# done

# --- D1b: Export Cicero CCANs as links file for arc tracks ---
# pyGenomeTracks links format: chr1 start1 end1 chr2 start2 end2 score
cicero_conns = pd.read_csv('cicero_connections.csv')

def peak_to_coords(peak_str):
    """Parse 'chr1-1000-2000' or 'chr1:1000-2000' to (chr, start, end)."""
    parts = peak_str.replace(':', '-').split('-')
    return parts[0], int(parts[1]), int(parts[2])

links_rows = []
for _, row in cicero_conns.iterrows():
    if row['coaccess'] < 0.25:
        continue
    chr1, s1, e1 = peak_to_coords(row['Peak1'])
    chr2, s2, e2 = peak_to_coords(row['Peak2'])
    if chr1 == chr2:  # pyGenomeTracks links must be same chromosome
        links_rows.append([chr1, s1, e1, chr2, s2, e2, row['coaccess']])

links_df = pd.DataFrame(links_rows,
    columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'score'])
links_df.to_csv('cicero_links.bedpe', sep='\t', header=False, index=False)

# --- D1c: Export DORC peak-gene links (if using multiome) ---
# Same format but scored by expression correlation instead of co-accessibility
dorc_links_rows = []
for _, row in dorc_sig.iterrows():
    chr_p, s_p, e_p = peak_to_coords(row['peak'])
    # Get TSS of linked gene
    gene_tss = tss_reference[tss_reference['gene'] == row['gene']].iloc[0]
    chr_g, s_g, e_g = gene_tss['chr'], gene_tss['start'], gene_tss['end']
    if chr_p == chr_g:
        dorc_links_rows.append([chr_p, s_p, e_p, chr_g, s_g, e_g, row['corr']])

dorc_links_df = pd.DataFrame(dorc_links_rows,
    columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'score'])
dorc_links_df.to_csv('dorc_links.bedpe', sep='\t', header=False, index=False)

# --- D1d: Export enhancer peak annotations as BED ---
# Highlight enhancer regions with binding evidence from footprinting
enhancer_bed = enhancer_peaks[['chr', 'start', 'end']].copy()
enhancer_bed['name'] = 'enhancer'
enhancer_bed['score'] = 1000
enhancer_bed['strand'] = '.'
enhancer_bed.to_csv('enhancer_annotations.bed', sep='\t', header=False, index=False)
```

---

## Step D2: pyGenomeTracks Configuration

pyGenomeTracks uses an INI-format config file to define track layout. Build this programmatically for the target locus.

```python
def write_pygenometracks_ini(
    config_path,
    gene_gtf,
    bigwig_files,       # dict: {cell_type_name: path_to_bigwig}
    links_file,         # path to cicero or DORC links bedpe
    enhancer_bed,       # path to enhancer annotation BED
    links_label='CCANs',
    colors=None,        # dict: {cell_type_name: color}
):
    """Generate a pyGenomeTracks INI config file."""

    lines = []

    # Gene model track
    lines.append('[genes]')
    lines.append(f'file = {gene_gtf}')
    lines.append('title = Genes')
    lines.append('height = 3')
    lines.append('fontsize = 10')
    lines.append('style = UCSC')
    lines.append('merge_transcripts = true')
    lines.append('')

    # Spacer
    lines.append('[spacer]')
    lines.append('height = 0.5')
    lines.append('')

    # CCAN / DORC arc links
    lines.append('[links]')
    lines.append(f'file = {links_file}')
    lines.append(f'title = {links_label}')
    lines.append('height = 3')
    lines.append('color = darkblue')
    lines.append('line_width = 1')
    lines.append('links_type = arcs')
    lines.append('orientation = inverted')
    lines.append('')

    # Enhancer annotation
    lines.append('[enhancers]')
    lines.append(f'file = {enhancer_bed}')
    lines.append('title = Enhancers')
    lines.append('height = 0.5')
    lines.append('color = #ff6b35')
    lines.append('border_color = none')
    lines.append('display = collapsed')
    lines.append('')

    # BigWig coverage tracks per cell type
    if colors is None:
        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab10')
        colors = {k: f'#{int(c[0]*255):02x}{int(c[1]*255):02x}{int(c[2]*255):02x}'
                  for k, c in zip(bigwig_files.keys(),
                                  [cmap(i) for i in range(len(bigwig_files))])}

    for cell_type, bw_path in bigwig_files.items():
        section_name = cell_type.replace(' ', '_').replace('/', '_')
        lines.append(f'[{section_name}]')
        lines.append(f'file = {bw_path}')
        lines.append(f'title = {cell_type}')
        lines.append('height = 2')
        lines.append(f'color = {colors.get(cell_type, "#333333")}')
        lines.append('min_value = 0')
        lines.append('number_of_bins = 500')
        lines.append('')

    # X-axis
    lines.append('[x-axis]')
    lines.append('fontsize = 10')
    lines.append('')

    with open(config_path, 'w') as f:
        f.write('\n'.join(lines))

    return config_path


# Example usage:
write_pygenometracks_ini(
    config_path='tracks.ini',
    gene_gtf='genes.gtf.gz',
    bigwig_files={
        'Monocyte':     'mono_coverage.bw',
        'T cell':       'tcell_coverage.bw',
        'B cell':       'bcell_coverage.bw',
    },
    links_file='cicero_links.bedpe',
    enhancer_bed='enhancer_annotations.bed',
    links_label='Cicero CCANs',
)
```

---

## Step D3: Render pyGenomeTracks Panel

```python
import subprocess

def render_pygenometracks(config_path, region_str, output_png, dpi=150, width=40):
    """
    Render genome browser tracks for a given region.
    region_str: e.g. 'chr11:5,200,000-5,350,000'
    """
    cmd = [
        'pyGenomeTracks',
        '--tracks', config_path,
        '--region', region_str,
        '--outFileName', output_png,
        '--dpi', str(dpi),
        '--width', str(width),
        '--fontSize', '10',
    ]
    subprocess.run(cmd, check=True)
    return output_png

# Render the browser tracks panel
browser_png = render_pygenometracks(
    'tracks.ini',
    'chr11:5,200,000-5,350,000',  # example: HBB locus
    'browser_tracks.png'
)
```

**Alternative: pyGenomeTracks Python API (for direct matplotlib integration)**

```python
# For tighter matplotlib integration, use the Python API directly
from pygenometracks.tracksClass import PlotTracks

# This returns a matplotlib figure you can composite with other panels
tracks = PlotTracks('tracks.ini', fig_width=40, dpi=150)
fig = tracks.plot('chr11', 5200000, 5350000)
# fig is a matplotlib Figure object — you can extract axes or save
```

---

## Step D4: Generate scprinter Panels

Capture scprinter plots as matplotlib figures. Two strategies depending on whether scprinter functions accept an `ax=` parameter:

```python
import scprinter as scp
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

# --- Strategy A: If scprinter functions accept ax= parameter ---
# (Check your scprinter version — newer versions may support this)
# fig_scp, axes_scp = plt.subplots(2, 1, figsize=(20, 8))
# scp.pl.plot_binding_score(printer, ..., ax=axes_scp[0])
# scp.pl.plot_footprints(printer, ..., ax=axes_scp[1])

# --- Strategy B: Capture scprinter output as images (always works) ---
def capture_scprinter_plot(plot_func, *args, **kwargs):
    """
    Call an scprinter plotting function and capture the result as an image.
    Returns a numpy array suitable for plt.imshow().
    """
    # scprinter creates its own figure; capture it
    plot_func(*args, **kwargs)
    fig = plt.gcf()
    fig.savefig('_temp_scp_panel.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return mpimg.imread('_temp_scp_panel.png')

# Capture TFBS heatmap
tfbs_img = capture_scprinter_plot(
    scp.pl.plot_binding_score,
    printer, save_key='enhancer_GATA1_TFBS',
    group_names=uniq_groups[order],
    kind='heatmap', region=target_region,
    row_label=cell_type_labels[order]
)

# Capture multiscale footprint
fp_img = capture_scprinter_plot(
    scp.pl.plot_footprints,
    printer, save_key='enhancer_GATA1_footprint',
    group_names=uniq_groups[order],
    row_label=cell_type_labels[order],
    region=target_region,
    stack=True, scales=[10, 20, 50, 100],
    cmap='Blues', vmin=0.5, vmax=2.0
)

# Capture Tn5 insertion profile
tn5_img = capture_scprinter_plot(
    scp.pl.plot_group_atac,
    printer, grouping[0:5], np.arange(5),
    target_region, smooth=5
)
```

---

## Step D5: Generate TF Motif Logo

```python
import logomaker
import pandas as pd

def render_motif_logo(motif_pwm, tf_name, ax=None):
    """
    Render a sequence logo from a position weight matrix.

    motif_pwm: pd.DataFrame with columns ['A', 'C', 'G', 'T']
               and rows = positions. Values = frequencies or information content.
    tf_name: string label for the motif.
    ax: matplotlib Axes to render into (if None, creates new figure).
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 1.5))

    # Convert frequencies to information content if needed
    # (logomaker can accept raw frequencies and type='information')
    logo = logomaker.Logo(
        motif_pwm,
        ax=ax,
        color_scheme='classic',    # A=green, C=blue, G=orange, T=red
        font_name='DejaVu Sans',
    )
    ax.set_title(tf_name, fontsize=12, fontweight='bold')
    ax.set_ylabel('bits')
    ax.set_xlim([-0.5, len(motif_pwm) - 0.5])
    return ax


# Example: load a motif PWM (from JASPAR, HOCOMOCO, or scprinter's motif set)
# JASPAR format or similar:
motif_pwm = pd.DataFrame({
    'A': [0.05, 0.01, 0.95, 0.01, 0.05, 0.80, 0.01, 0.01],
    'C': [0.05, 0.01, 0.01, 0.01, 0.80, 0.05, 0.01, 0.01],
    'G': [0.80, 0.01, 0.01, 0.95, 0.10, 0.10, 0.01, 0.95],
    'T': [0.10, 0.97, 0.03, 0.03, 0.05, 0.05, 0.97, 0.03],
})  # Example: GATA motif WGATAA
```

---

## Step D6: Assemble Composite Figure

```python
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

def assemble_composite_figure(
    browser_png,         # path to pyGenomeTracks output
    tfbs_img,            # numpy array from capture_scprinter_plot
    footprint_img,       # numpy array from capture_scprinter_plot
    tn5_img,             # numpy array from capture_scprinter_plot
    motif_pwm,           # pd.DataFrame for logomaker
    tf_name,             # string
    enhancer_fp_img=None, # optional: footprint at distal enhancer
    title='',
    save_path=None,
    figsize=(24, 28),
):
    """
    Assemble all panels into a single composite figure.
    """
    # Load browser tracks image
    browser_img = mpimg.imread(browser_png)

    # Define grid layout
    fig = plt.figure(figsize=figsize, constrained_layout=True)

    if enhancer_fp_img is not None:
        gs = gridspec.GridSpec(6, 2, figure=fig,
            height_ratios=[3, 2, 2, 2, 2, 1.5],
            width_ratios=[1, 1])
    else:
        gs = gridspec.GridSpec(5, 2, figure=fig,
            height_ratios=[3, 2, 2, 2, 1.5],
            width_ratios=[1, 1])

    row = 0

    # Panel A: Browser tracks (gene model + CCANs + coverage) — full width
    ax_browser = fig.add_subplot(gs[row, :])
    ax_browser.imshow(browser_img, aspect='auto')
    ax_browser.axis('off')
    ax_browser.set_title('A. Genomic Context: Gene Model, CCANs, ATAC Coverage',
                         fontsize=14, fontweight='bold', loc='left')
    row += 1

    # Panel B: Tn5 insertion profile — full width
    ax_tn5 = fig.add_subplot(gs[row, :])
    ax_tn5.imshow(tn5_img, aspect='auto')
    ax_tn5.axis('off')
    ax_tn5.set_title('B. Tn5 Insertion Profile by Cell Type',
                     fontsize=14, fontweight='bold', loc='left')
    row += 1

    # Panel C: TFBS binding score heatmap — full width
    ax_tfbs = fig.add_subplot(gs[row, :])
    ax_tfbs.imshow(tfbs_img, aspect='auto')
    ax_tfbs.axis('off')
    ax_tfbs.set_title('C. TF Binding Score (Promoter Region)',
                     fontsize=14, fontweight='bold', loc='left')
    row += 1

    # Panel D: Multiscale footprint heatmap — full width
    ax_fp = fig.add_subplot(gs[row, :])
    ax_fp.imshow(footprint_img, aspect='auto')
    ax_fp.axis('off')
    ax_fp.set_title('D. Multiscale Footprint',
                     fontsize=14, fontweight='bold', loc='left')
    row += 1

    # Panel E (left): TF motif logo
    # Panel F (right): Enhancer footprint (if available)
    if enhancer_fp_img is not None:
        # Bottom row split: logo left, enhancer footprint right
        ax_logo = fig.add_subplot(gs[row, 0])
        render_motif_logo(motif_pwm, tf_name, ax=ax_logo)
        ax_logo.set_title('E. TF Motif', fontsize=14, fontweight='bold', loc='left')

        ax_enh_fp = fig.add_subplot(gs[row, 1])
        ax_enh_fp.imshow(enhancer_fp_img, aspect='auto')
        ax_enh_fp.axis('off')
        ax_enh_fp.set_title('F. Enhancer Footprint (Distal)',
                           fontsize=14, fontweight='bold', loc='left')
    else:
        # Just the logo, centered
        ax_logo = fig.add_subplot(gs[row, :])
        render_motif_logo(motif_pwm, tf_name, ax=ax_logo)
        ax_logo.set_title('E. TF Motif', fontsize=14, fontweight='bold', loc='left')

    if title:
        fig.suptitle(title, fontsize=18, fontweight='bold', y=1.01)

    if save_path:
        fig.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f'Composite figure saved to: {save_path}')

    return fig


# --- Example invocation ---
fig = assemble_composite_figure(
    browser_png='browser_tracks.png',
    tfbs_img=tfbs_img,
    footprint_img=fp_img,
    tn5_img=tn5_img,
    motif_pwm=motif_pwm,
    tf_name='GATA1',
    enhancer_fp_img=enhancer_fp_img,  # from footprinting a distal DORC enhancer
    title='Regulatory Architecture: HBB Locus — GATA1 Binding',
    save_path='composite_HBB_GATA1.png',
)
```

---

## Step D7: Batch Generation Across Loci (Agent Loop)

For systematic analysis, generate composite figures for all DORC genes or eRegulon targets:

```python
def generate_composite_for_gene(
    gene_name, tf_name, printer, adata_rna,
    dorc_sig, eregulon_region_sets,
    grouping, uniq_groups, order, cell_type_labels,
    genome, tracks_ini_template,
    output_dir='composite_figures'
):
    """
    Generate a composite figure for one gene locus.
    Orchestrates: pyGenomeTracks rendering, scprinter panel capture,
    motif logo generation, and final assembly.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    # 1. Get gene coordinates (TSS +/- window)
    gene_info = tss_reference[tss_reference['gene'] == gene_name].iloc[0]
    chrom = gene_info['chr']
    tss = gene_info['start']
    window = 100000  # 100kb window
    region_str = f"{chrom}:{max(0, tss - window)}-{tss + window}"

    # 2. Render browser tracks
    browser_png = f'{output_dir}/{gene_name}_browser.png'
    render_pygenometracks(tracks_ini_template, region_str, browser_png)

    # 3. Define the promoter region for scprinter
    promoter_region = pd.DataFrame([{
        'Chromosome': chrom,
        'Start': max(0, tss - 2000),
        'End': tss + 2000
    }])

    # 4. Capture scprinter panels (Tn5, TFBS, footprint)
    # ... (use capture_scprinter_plot as in Step D4) ...

    # 5. Get enhancer regions for this TF from eRegulon/DORC
    # ... (lookup from existing region sets) ...

    # 6. Capture enhancer footprint panel
    # ... (scprinter footprint at distal enhancer region) ...

    # 7. Load TF motif PWM
    # ... (from JASPAR or scprinter motif database) ...

    # 8. Assemble
    fig = assemble_composite_figure(
        browser_png=browser_png,
        tfbs_img=tfbs_img,
        footprint_img=fp_img,
        tn5_img=tn5_img,
        motif_pwm=motif_pwm,
        tf_name=tf_name,
        enhancer_fp_img=enhancer_fp_img,
        title=f'Regulatory Architecture: {gene_name} — {tf_name}',
        save_path=f'{output_dir}/composite_{gene_name}_{tf_name}.png',
    )
    plt.close(fig)

# Batch: generate for all DORC genes with footprinted TFs
for gene in dorc_gene_list:
    for tf in footprinted_tfs_for_gene.get(gene, []):
        generate_composite_for_gene(gene, tf, ...)
```

---

## Summary Table

| Step | Tool | Input | Output |
|------|------|-------|--------|
| D1 | deeptools / pandas | BAMs, Cicero conns, DORC links | bigWig, BEDPE, BED files |
| D2 | pyGenomeTracks | Track files | INI config |
| D3 | pyGenomeTracks | INI + region | Browser tracks PNG or matplotlib Figure |
| D4 | scprinter.pl + matplotlib | printer + regions | Captured panel images |
| D5 | logomaker | PWM DataFrame | Sequence logo in matplotlib Axes |
| D6 | matplotlib gridspec | All panels | Composite figure (PNG/PDF/SVG) |
| D7 | Batch loop | Gene/TF lists | Directory of composite figures |

---

## Technical Notes

1. **Coordinate alignment:** pyGenomeTracks handles coordinate alignment natively for its own tracks. For scprinter panels (captured as images), alignment with the browser coordinate system is visual rather than pixel-exact. For publication figures, ensure the region window matches between pyGenomeTracks and scprinter calls.
2. **Vector vs. raster:** pyGenomeTracks can output SVG/PDF for vector graphics. scprinter captured panels are raster (PNG). For the highest quality, use the pyGenomeTracks Python API to render directly into matplotlib axes, then overlay scprinter panels. Alternatively, export all panels as high-DPI PNG (300+) and composite.
3. **HPC considerations:** pyGenomeTracks requires a display backend. On headless HPC, set `matplotlib.use('Agg')` before import. All rendering is file-based (no interactive display needed).
4. **IGV alternative:** If pyGenomeTracks is not available, IGV batch mode can replace Steps D2-D3. Use: `echo "goto chr11:5,200,000-5,350,000\nsnapshot browser_tracks.png\nexit" | igv -b /dev/stdin`. However, pyGenomeTracks is preferred for programmatic reproducibility and matplotlib integration.

---

## Cross-References to Other Recipes

- **Recipe A outputs used here:** Cicero CCANs (Step D1b), chromVAR clusters for pseudobulk grouping, scprinter TFBS/footprint scores at enhancer regions
- **Recipe B outputs used here:** DORC peak-gene links (Step D1c), eRegulon enhancer annotations (Step D1d), RNA-informed cell groupings
- **Recipe C outputs used here:** Signaling-validated TF targets can be highlighted in the composite figure as a distinct annotation track or color-coded in the TFBS heatmap
- **scprinter API:** `scp.pl.plot_binding_score()`, `scp.pl.plot_footprints()`, `scp.pl.plot_group_atac()`, `scp.pl.plot_region_atac()` — all documented in the main recipe (enhancer_footprinting_recipe.md)
