#!/usr/bin/env python3
"""
run_scprinter_footprinting.py — FIX-38/42 Rewrite

Runs scPrinter TF footprinting in two modes:

  GLOBAL MODE (no --control-condition / --treatment-condition):
    FIX-38: Rewritten for per-cell-type fan-out.  Each Nextflow task receives
    ONE cell type via --cell-type-label.  Groups barcodes for the FOCAL cell
    type only, computes multiscale footprints + TFBS binding scores at each
    target gene's promoter, and produces three publication-quality plots per
    gene:
      Plot A — Tn5 insertion profile (raw chromatin accessibility)
      Plot B — Multiscale footprint heatmap (scale x position, the key
               scprinter visualization)
      Plot C — TFBS binding score heatmap (TF occupancy prediction)

    The old run_global_mode() grouped ALL cell types and plotted a stacked
    heatmap, which was nonsensical in the per-cell-type fan-out context
    (cell-type-specific TFs were being displayed across unrelated cell types).

  DIFFERENTIAL MODE (--control-condition AND --treatment-condition):
    FIX-42: All plots emphasize the contrast between conditions:
      - Differential MSFP heatmap (treatment − control, RdBu_r) with TSS line
      - Single-scale line plot (overlay at --diff-scale, same axes)
      - Differential TFBS binding score (treatment − control, fill_between)
    Groups barcodes by condition within a single --cell-type.
    Plus enhancer footprinting stub for future module integration.

FIX-42: Integrated SDas refactor improvements:
    - Precise TSS from gene_coordinate_resolver 'tss' field (was midpoint)
    - Auto-detection of barcode normalization mode (choose_normalization)
    - Differential heatmap (RdBu_r, treatment-control) + TSS vertical line
    - Differential single-scale line plot at user-selected scale
    - sample_id_regex parameter for flexible barcode-to-sample mapping

Outputs:
    - footprints_<label>_<gene>.h5ad
    - tfbs_<label>_<gene>.h5ad  (global mode, FIX-38)
    - enhancer_footprints_<label>_<gene>.h5ad  (differential only)
    - plots/<label>_<gene>_*.png
    - footprint_summary.csv
    - enhancer_footprint_summary.csv  (differential only)
"""

import os
import re
import sys
import argparse
import builtins
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

# FIX-44: scPRINTER's get_binding_score() uses Python 2's `long` type internally,
# which was removed in Python 3.  Inject `long = int` into builtins so the
# library finds it at runtime.  This is the standard Py2→3 compatibility shim.
builtins.long = int

import scprinter as scp
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import logomaker
except ImportError:
    logomaker = None


# =============================================================================
# Helpers
# =============================================================================
def init_scprinter(cache_dir: str):
    """Set scprinter dataset cache path."""
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp


def parse_args():
    p = argparse.ArgumentParser(
        description="Run scPrinter footprinting (FIX-38: focal cell-type + TFBS scoring)"
    )
    p.add_argument("--peak-matrix", required=True,
                   help="Path to filtered peak matrix .h5ad")
    p.add_argument("--cell-type", required=True,
                   help="Cell type name (differential mode) or 'all' (global mode)")
    p.add_argument("--target-genes", required=True,
                   help="Comma-separated list of target genes")
    p.add_argument("--cache-dir", required=True,
                   help="Path to scPrinter cache directory")
    p.add_argument("--genome", required=True,
                   help="Genome key for scPrinter (e.g. 'hg38', 'mm10')")
    p.add_argument("--printer-path", required=True,
                   help="Path to existing scPrinter printer .h5ad")
    p.add_argument("--coordinates-file",
                   help="Path to gene coordinates JSON from gene_coordinate_resolver.py")

    # FIX-23b: Label for output filenames (e.g. cell type name in discovery mode)
    p.add_argument("--cell-type-label", default=None,
                   help="Label for output file prefixes (e.g. 'Astrocytes' in discovery mode)")

    # FIX-23: Made optional — omit both for global mode
    p.add_argument("--control-condition", default=None,
                   help="Control condition label (omit for global mode)")
    p.add_argument("--treatment-condition", default=None,
                   help="Treatment condition label (omit for global mode)")

    p.add_argument("--cpus", type=int, default=1,
                   help="Number of jobs/threads for footprinting")
    p.add_argument("--scale", type=int, default=50,
                   help="Scale for single-scale insertion profile plots")

    # FIX-38: TFBS binding score toggle
    p.add_argument("--skip-tfbs", action="store_true", default=False,
                   help="Skip TFBS binding score computation (faster, fewer plots)")

    # FIX-42: SDas-style parameters
    p.add_argument("--sample-id-regex", default=None,
                   help="Regex with named group 'sample' for extracting sample ID "
                        "from barcode (e.g. '(?P<sample>.+?)_\\d+$')")
    p.add_argument("--diff-scale", type=int, default=50,
                   help="Scale index for differential single-scale line plot (default: 50)")
    p.add_argument("--tf-mapping", default=None,
                   help="Path to tf_target_genes.json (maps genes back to their regulating TFs)")
    p.add_argument("--cicero-connections", default="",
                   help="Cicero connections TSV.gz for CCAN arc panel")
    p.add_argument("--pfm-path", default="",
                   help="JASPAR PFM file for TF motif logo panel")
    # FIX-46: Configurable cell type column name (was hardcoded 'cell_type')
    p.add_argument("--cell-type-col", default="celltypist_prediction",
                   help="obs column in peak matrix containing cell type labels "
                        "(default: celltypist_prediction)")
    return p.parse_args()


def normalize_peak_barcode(bc: str) -> str:
    """Normalize peak-matrix barcode to match printer.obs_names.

    FIX-37: Peak matrix uses colon delimiter (e.g. 'L1_Donor_0_july:1000602')
    while the printer uses underscore (e.g. 'L1_Donor_0_july_1000602').
    Simply replace the first colon with underscore.
    """
    if ":" in bc:
        return bc.replace(":", "_", 1)
    return bc


def tail_only(bc: str) -> str:
    """Strip sample prefix, keeping only the cell barcode portion.
    Check colon first — our pipeline uses sample_id:barcode format,
    and sample_id itself may contain underscores (e.g. Mouse_Kidney_BD).
    """
    if ":" in bc:
        return bc.split(":", 1)[-1]
    if "_" in bc:
        return bc.split("_", 1)[-1]
    return bc


# FIX-37e: 10x barcode regex — matches ACGTN+-\d+ at the end of a string,
# preceded by a colon or underscore separator.  Handles sample names that
# contain underscores (e.g. '10k_PBMC:ACGTACGT-1' or '10k_PBMC_ACGTACGT-1').
_10X_BC_PATTERN = re.compile(r'[_:]([ACGTN]+-\d+)$')


def strip_10x(bc: str) -> str:
    """FIX-37e: Extract the 10x-style barcode from a prefixed barcode string.

    Examples:
        '10k_PBMC:ACGTACGT-1'    → 'ACGTACGT-1'
        '10k_PBMC_ACGTACGT-1'    → 'ACGTACGT-1'
        'L1_Donor_0_july:100602' → '100602'  (BD style, less relevant)
        'ACGTACGT-1'             → 'ACGTACGT-1'  (no prefix, returned as-is)
    """
    m = _10X_BC_PATTERN.search(bc)
    if m:
        return m.group(1)
    return bc


def choose_normalization(peak_barcodes: np.ndarray, printer_barcodes: set) -> str:
    """FIX-42: Auto-detect which barcode normalization maps peak→printer barcodes.

    Probes a sample of peak-matrix barcodes with three strategies:
      - 'head_tail': replace first colon with underscore (FIX-37 default for BD)
      - 'tail_only': strip everything before the first colon/underscore
      - '10x_strip': regex-extract the 10x barcode (ACGTN+-\\d+) from the end
                     (FIX-37e: handles sample names with underscores like '10k_PBMC')

    Returns the strategy name that produces the most matches.  Falls back to
    'head_tail' on ties (the safer default for BD Rhapsody data).
    """
    sample = peak_barcodes[:min(500, len(peak_barcodes))]

    n_head_tail = sum(1 for bc in sample if normalize_peak_barcode(bc) in printer_barcodes)
    n_tail_only = sum(1 for bc in sample if tail_only(bc) in printer_barcodes)
    n_10x_strip = sum(1 for bc in sample if strip_10x(bc) in printer_barcodes)

    counts = {'head_tail': n_head_tail, 'tail_only': n_tail_only, '10x_strip': n_10x_strip}
    chosen = max(counts, key=counts.get)
    # On all-zero ties, fall back to head_tail (BD default)
    if counts[chosen] == 0:
        chosen = 'head_tail'
    print(f"  FIX-42 choose_normalization: head_tail={n_head_tail}, tail_only={n_tail_only}, "
          f"10x_strip={n_10x_strip} (of {len(sample)} sampled) → using '{chosen}'")
    return chosen


def apply_normalization(bc: str, mode: str) -> str:
    """Apply the chosen barcode normalization strategy."""
    if mode == "tail_only":
        return tail_only(bc)
    if mode == "10x_strip":
        return strip_10x(bc)
    return normalize_peak_barcode(bc)


def load_gene_coordinates(coord_file: str) -> Dict:
    """Load gene coordinates from resolver output.

    FIX-42: Now extracts the precise TSS position from the 'tss' field
    output by gene_coordinate_resolver.py, instead of approximating as
    the midpoint of the promoter window.  The TSS field is strand-aware
    (start for + strand, end for - strand).
    """
    with open(coord_file, 'r') as f:
        coords = json.load(f)
    gene_regions = {}
    for gene, info in coords.items():
        if info.get('found') and 'promoter' in info:
            p = info['promoter']
            region = {
                "chr": p['chr'],
                "start": p['start'],
                "end": p['end'],
            }
            # FIX-42: Use precise TSS from resolver (strand-aware)
            if 'tss' in info:
                region["tss"] = int(info['tss'])
            else:
                # Fallback: midpoint of promoter window
                region["tss"] = (int(p['start']) + int(p['end'])) // 2
            gene_regions[gene] = region
    return gene_regions


def _get_tss(region: Dict) -> int:
    """FIX-42: Get TSS position from region dict.

    Uses the precise 'tss' field from gene_coordinate_resolver when available,
    otherwise falls back to midpoint of the promoter window.
    """
    if "tss" in region:
        return int(region["tss"])
    return (int(region["start"]) + int(region["end"])) // 2


def _format_msfp_axes(ax, region: Dict, gene: str, modes=None,
                      xlabel=True, ylabel=True):
    """FIX-41/42/44b: Post-process multiscale footprint heatmap axes.

    scprinter's plot_footprints renders raw array indices on both axes.
    This function relabels them to be publication-quality:
      - X-axis: TSS-centered distance in bp (e.g. -1000 to +1000)
      - Y-axis: Footprint scale in bp (2, 10, 20, 50, 100)
      - TSS vertical line at the exact TSS position

    FIX-42: Uses precise TSS from gene_coordinate_resolver instead of
    midpoint approximation.
    FIX-44b: 7 ticks, proper pixel-to-genomic coordinate mapping,
    TSS vertical line, warning messages instead of silent pass.

    Parameters:
        ax: matplotlib Axes from plot_footprints
        region: dict with 'chr', 'start', 'end', optionally 'tss'
        gene: gene name for title
        modes: np.array of scales used (default: np.arange(2, 101))
        xlabel: whether to set x-axis label
        ylabel: whether to set y-axis label
    """
    if modes is None:
        modes = np.arange(2, 101)

    start = int(region["start"])
    end = int(region["end"])
    tss = _get_tss(region)

    # X-axis: TSS-centered
    try:
        xlim = ax.get_xlim()
        n_ticks = 7  # FIX-44b: 7 ticks for finer resolution
        # Generate tick positions evenly across the x-axis pixel extent
        tick_positions = np.linspace(xlim[0], xlim[1], n_ticks)
        # Map pixel positions to genomic coordinates via linear interpolation
        tick_genomic = start + (tick_positions - xlim[0]) / (xlim[1] - xlim[0]) * (end - start)
        tick_labels = [f"{int(g - tss):+d}" for g in tick_genomic]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        if xlabel:
            ax.set_xlabel(f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                          fontsize=9)

        # FIX-44b: TSS vertical line
        tss_pixel = xlim[0] + (tss - start) / (end - start) * (xlim[1] - xlim[0])
        if xlim[0] <= tss_pixel <= xlim[1]:
            ax.axvline(tss_pixel, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
    except Exception as e:
        print(f"  [WARN] _format_msfp_axes x-axis formatting failed: {e}")

    # Y-axis: Scale labels
    try:
        y_lo, y_hi = ax.get_ylim()  # FIX-44b: clearer variable names
        # modes go from 2 to 100; map representative scales to axis positions
        scale_ticks = [2, 10, 20, 50, 100]
        n_modes = len(modes)
        # Positions in axis coordinates (scale index within modes array)
        tick_positions_y = []
        tick_labels_y = []
        for s in scale_ticks:
            if s in modes:
                idx = np.where(modes == s)[0][0]
                # Map index to axis coordinates
                frac = idx / (n_modes - 1) if n_modes > 1 else 0.5
                pos = y_lo + frac * (y_hi - y_lo)
                tick_positions_y.append(pos)
                tick_labels_y.append(f"{s} bp")
        if tick_positions_y:
            ax.set_yticks(tick_positions_y)
            ax.set_yticklabels(tick_labels_y, fontsize=8)
        if ylabel:
            ax.set_ylabel("Footprint Scale", fontsize=9)
    except Exception as e:
        print(f"  [WARN] _format_msfp_axes y-axis formatting failed: {e}")


def _format_tss_xaxis(ax, region: Dict, gene: str, add_vline: bool = True):
    """FIX-41/42/44b: Add TSS-centered x-axis labels to any plot over a promoter region.

    FIX-42: Uses precise TSS from gene_coordinate_resolver.
    FIX-44b: Optional TSS vertical line, warning messages instead of silent pass.

    Parameters:
        ax: matplotlib Axes
        region: dict with 'chr', 'start', 'end', optionally 'tss'
        gene: gene name for axis label
        add_vline: if True, draw a vertical dashed line at TSS position
    """
    start = int(region["start"])
    end = int(region["end"])
    tss = _get_tss(region)

    try:
        xlim = ax.get_xlim()
        n_ticks = 7  # FIX-44b: 7 ticks for consistency with _format_msfp_axes
        tick_positions = np.linspace(xlim[0], xlim[1], n_ticks)
        tick_genomic = start + (tick_positions - xlim[0]) / (xlim[1] - xlim[0]) * (end - start)
        tick_labels = [f"{int(g - tss):+d}" for g in tick_genomic]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel(f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                      fontsize=9)

        # FIX-44b: TSS vertical line
        if add_vline:
            tss_pixel = xlim[0] + (tss - start) / (end - start) * (xlim[1] - xlim[0])
            if xlim[0] <= tss_pixel <= xlim[1]:
                ax.axvline(tss_pixel, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
    except Exception as e:
        print(f"  [WARN] _format_tss_xaxis formatting failed: {e}")


def _init_uns_dict(printer, key: str):
    """FIX-39: Convert an HDF5-backed printer.uns entry to a plain Python dict.

    The Rust anndata backend (anndata-0.4.2) keeps internal HDF5 group handle
    caches.  When scprinter writes a new footprint/binding-score entry to
    printer.uns[key] during gene N, the cache becomes stale.  On gene N+1 the
    check ``save_key in printer.uns[key]`` calls H5Gopen2() on the stale
    handle, triggering a Rust panic.

    Converting to a plain dict before the gene loop makes all subsequent
    reads/writes go through Python, sidestepping the HDF5 cache entirely.
    """
    try:
        if key in printer.uns:
            printer.uns[key] = dict(printer.uns[key])
        else:
            printer.uns[key] = {}
        print(f"  FIX-39: printer.uns['{key}'] → Python dict ({len(printer.uns[key])} entries)")
    except Exception as e:
        # If even reading the existing dict fails (cache already corrupt from
        # a previous run), replace with empty dict so the loop can proceed.
        print(f"  FIX-39: Could not read printer.uns['{key}'] ({e}); resetting to empty dict")
        printer.uns[key] = {}


def build_focal_grouping(
    adata: ad.AnnData,
    printer,
    cell_type_label: str,
    bc_norm_mode: Optional[str] = None,
    cell_type_col: str = "celltypist_prediction",
) -> Tuple[List, np.ndarray]:
    """Build barcode grouping for a single focal cell type.

    FIX-38: Filters peak matrix to the focal cell type, normalizes barcodes,
    intersects with printer.obs_names, and returns a single-group grouping.

    FIX-42: Uses choose_normalization() auto-detection when bc_norm_mode is
    None.  Supports 'head_tail' (FIX-37 default) and 'tail_only' (SDas style).

    FIX-46: cell_type_col is configurable (was hardcoded 'cell_type').

    Returns:
        grouping: list with one element (list of barcode strings)
        uniq_groups: np.array([cell_type_label])
    """
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"adata.obs does not contain a '{cell_type_col}' column. "
                         f"Available: {list(adata.obs.columns)}")

    # Filter to focal cell type
    cell_mask = adata.obs[cell_type_col].astype(str) == cell_type_label
    n_cells = int(cell_mask.sum())
    print(f"  Focal cell type '{cell_type_label}': {n_cells} cells in peak matrix")

    if n_cells == 0:
        available = sorted(adata.obs[cell_type_col].astype(str).unique())
        raise ValueError(
            f"No cells found for cell_type '{cell_type_label}'. "
            f"Available: {available}"
        )

    # FIX-42: Auto-detect normalization mode if not specified
    bcs_raw = adata.obs_names[cell_mask].astype(str)
    printer_bcs = set(map(str, printer.obs_names))

    if bc_norm_mode is None:
        bc_norm_mode = choose_normalization(np.asarray(bcs_raw), printer_bcs)

    bcs_normalized = np.array(
        [apply_normalization(b, bc_norm_mode) for b in bcs_raw], dtype=str
    )
    bcs_common = [b for b in bcs_normalized if b in printer_bcs]

    if not bcs_common:
        # FIX-37e: Debug output to diagnose barcode format mismatches
        sample_peak = list(bcs_raw[:3])
        sample_printer = list(printer_bcs)[:3]
        sample_normalized = list(bcs_normalized[:3])
        raise ValueError(
            f"No barcodes for cell_type '{cell_type_label}' found in printer. "
            f"Peak matrix has {n_cells} cells but none match printer.obs_names. "
            f"Normalization mode: {bc_norm_mode}. "
            f"Sample peak barcodes: {sample_peak}. "
            f"After normalization: {sample_normalized}. "
            f"Sample printer barcodes: {sample_printer}"
        )

    print(f"  Matched {len(bcs_common)} / {n_cells} barcodes in printer")

    grouping = [bcs_common]
    uniq_groups = np.array([cell_type_label])
    return grouping, uniq_groups


# =============================================================================
# Helper functions for Cicero CCAN arcs and TF motif logos
# =============================================================================
def load_jaspar_pwm(pfm_path, tf_name):
    """Load a PWM from a JASPAR PFM file for the given TF name.

    Returns a pandas DataFrame with columns A, C, G, T (information content)
    or None if TF not found.
    """
    if not pfm_path or not os.path.isfile(pfm_path):
        return None

    try:
        from Bio import motifs
        # Read all motifs from JASPAR file
        with open(pfm_path) as f:
            all_motifs = motifs.parse(f, 'jaspar')

        # Try exact match on name first, then case-insensitive
        tf_upper = tf_name.upper()
        # Handle composite names like "FOS::JUND" — try both the full name and parts
        candidates = [tf_name, tf_upper]
        if '::' in tf_name:
            candidates.extend(tf_name.split('::'))
            candidates.extend(tf_upper.split('::'))

        for m in all_motifs:
            m_name = m.name.upper() if m.name else ''
            for cand in candidates:
                if m_name == cand.upper():
                    # Convert counts to PWM (position probability matrix)
                    pwm = m.counts.normalize(pseudocounts=0.5)
                    # Convert to information content for logomaker
                    rows = []
                    for i in range(len(pwm['A'])):
                        row = {'A': pwm['A'][i], 'C': pwm['C'][i],
                               'G': pwm['G'][i], 'T': pwm['T'][i]}
                        rows.append(row)
                    pwm_df = pd.DataFrame(rows)
                    # Convert to information content (bits)
                    ic = pwm_df.copy()
                    for col in ['A', 'C', 'G', 'T']:
                        ic[col] = pwm_df[col] * np.log2(pwm_df[col] / 0.25 + 1e-10)
                        ic[col] = ic[col].clip(lower=0)
                    return ic

        print(f"  [WARN] TF '{tf_name}' not found in JASPAR PFM file")
        return None
    except ImportError:
        print("  [WARN] Biopython not available — trying manual JASPAR parse")
        return _parse_jaspar_manual(pfm_path, tf_name)
    except Exception as e:
        print(f"  [WARN] Failed to load PWM for {tf_name}: {e}")
        return None


def _parse_jaspar_manual(pfm_path, tf_name):
    """Fallback manual parser for JASPAR format PFMs."""
    tf_upper = tf_name.upper()
    candidates = {tf_upper}
    if '::' in tf_name:
        candidates.update(p.upper() for p in tf_name.split('::'))

    try:
        with open(pfm_path) as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith('>'):
                # Header line: >MA0099.3  FOS::JUND
                parts = line[1:].split()
                name = parts[1] if len(parts) > 1 else parts[0]
                if name.upper() in candidates:
                    # Read next 4 lines (A, C, G, T counts)
                    counts = {}
                    for nuc in ['A', 'C', 'G', 'T']:
                        i += 1
                        if i >= len(lines):
                            break
                        row = lines[i].strip()
                        # Format: A [ 1.0  2.0  3.0 ... ]
                        vals_str = row.split('[')[1].split(']')[0] if '[' in row else row
                        vals = [float(v) for v in vals_str.split()]
                        counts[nuc] = vals

                    if len(counts) == 4:
                        n_pos = len(counts['A'])
                        rows = []
                        for pos in range(n_pos):
                            total = sum(counts[nuc][pos] for nuc in 'ACGT') + 2.0  # pseudocount
                            row = {}
                            for nuc in 'ACGT':
                                p_val = (counts[nuc][pos] + 0.5) / total
                                row[nuc] = p_val * np.log2(p_val / 0.25 + 1e-10)
                                row[nuc] = max(0, row[nuc])
                            rows.append(row)
                        return pd.DataFrame(rows)
            i += 1

        return None
    except Exception as e:
        print(f"  [WARN] Manual JASPAR parse failed: {e}")
        return None


def render_motif_logo(pwm_df, tf_name, ax):
    """Render a sequence logo from information content DataFrame."""
    try:
        import logomaker
        logomaker.Logo(pwm_df, ax=ax, color_scheme='classic', font_name='DejaVu Sans')
        ax.set_title(f'TF Motif: {tf_name}', fontsize=10, fontweight='bold')
        ax.set_ylabel('bits', fontsize=8)
        ax.set_xlim([-0.5, len(pwm_df) - 0.5])
        return True
    except ImportError:
        ax.text(0.5, 0.5, f'{tf_name}\n(logomaker not available)',
                ha='center', va='center', fontsize=10, transform=ax.transAxes)
        ax.set_title(f'TF Motif: {tf_name}', fontsize=10, fontweight='bold')
        ax.axis('off')
        return False
    except Exception as e:
        print(f"  [WARN] Logo rendering failed: {e}")
        ax.text(0.5, 0.5, f'{tf_name}\n(logo failed)',
                ha='center', va='center', fontsize=10, transform=ax.transAxes)
        ax.axis('off')
        return False


def plot_cicero_arcs(cicero_conns, best_region_str, regions_df, ax, window_kb=500):
    """Plot Cicero connection arcs centered on the best enhancer region.

    Shows co-accessibility arcs within +/-window_kb of the best region.
    """
    from matplotlib.patches import Arc as MplArc

    if cicero_conns is None or cicero_conns.empty:
        ax.text(0.5, 0.5, 'No Cicero connections available',
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')
        return

    # Parse best region
    chrom, rest = best_region_str.split(':')
    center_start, center_end = [int(x) for x in rest.split('-')]
    center = (center_start + center_end) // 2
    window = window_kb * 1000

    view_start = center - window
    view_end = center + window

    # Filter connections to those within the window on the same chromosome
    mask = (
        (cicero_conns['p1_chrom'] == chrom) &
        (cicero_conns['p2_chrom'] == chrom) &
        (cicero_conns['p1_start'] >= view_start) &
        (cicero_conns['p1_end'] <= view_end) &
        (cicero_conns['p2_start'] >= view_start) &
        (cicero_conns['p2_end'] <= view_end)
    )
    local = cicero_conns[mask].copy()

    if local.empty:
        ax.text(0.5, 0.5, f'No Cicero connections within +/-{window_kb}kb',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'CCAN arcs at {best_region_str}', fontsize=10)
        ax.axis('off')
        return

    # Build set of region midpoints for highlighting
    enh_mids = set()
    for _, row in regions_df.iterrows():
        if row['Chromosome'] == chrom:
            mid = (row['Start'] + row['End']) // 2
            if view_start <= mid <= view_end:
                enh_mids.add(mid)

    # Normalize co-accessibility for color mapping
    max_coaccess = local['coaccess'].max()
    min_coaccess = local['coaccess'].min()

    for _, row in local.iterrows():
        mid1 = (row['p1_start'] + row['p1_end']) // 2
        mid2 = (row['p2_start'] + row['p2_end']) // 2
        x_center = (mid1 + mid2) / 2
        width = abs(mid2 - mid1)
        height = width * 0.4  # arc height proportional to distance

        # Color by co-accessibility score
        if max_coaccess > min_coaccess:
            norm = (row['coaccess'] - min_coaccess) / (max_coaccess - min_coaccess)
        else:
            norm = 0.5
        color = plt.cm.RdYlBu_r(norm)
        alpha = 0.3 + 0.7 * norm

        arc = MplArc((x_center, 0), width, height, angle=0,
                      theta1=0, theta2=180, color=color, alpha=alpha, linewidth=0.8)
        ax.add_patch(arc)

    # Mark region midpoints as ticks on the baseline
    for mid in enh_mids:
        ax.plot(mid, 0, '|', color='red', markersize=6, markeredgewidth=1.5)

    # Mark the best region
    ax.axvspan(center_start, center_end, alpha=0.15, color='red', zorder=0)

    ax.set_xlim(view_start, view_end)
    max_arc_height = max(abs(r['p2_start'] - r['p1_start']) * 0.2 for _, r in local.iterrows()) if len(local) > 0 else window * 0.1
    ax.set_ylim(-max_arc_height * 0.1, max_arc_height * 1.2)
    ax.set_xlabel(f'{chrom} position', fontsize=8)
    ax.set_title(f'CCAN arcs at {best_region_str} (+/-{window_kb}kb, {len(local)} connections)', fontsize=10)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_yticks([])

    # FIX-88: Add colorbar legend for co-accessibility score (match enhancer style)
    import matplotlib.colors as mcolors
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.RdYlBu_r,
        norm=mcolors.Normalize(vmin=min_coaccess, vmax=max_coaccess))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, aspect=15)
    cbar.set_label('Cicero co-accessibility', fontsize=8)
    cbar.ax.tick_params(labelsize=7)


def load_all_cicero_connections(conns_path):
    """Load all Cicero connections without enhancer filtering."""
    try:
        conns = pd.read_csv(conns_path, sep='\t', compression='infer')
        for prefix, col in [('p1', 'Peak1'), ('p2', 'Peak2')]:
            parts = conns[col].str.replace(':', '-').str.split('-', expand=True)
            conns[f'{prefix}_chrom'] = parts[0]
            conns[f'{prefix}_start'] = parts[1].astype(int)
            conns[f'{prefix}_end'] = parts[2].astype(int)
        print(f"  Loaded {len(conns)} Cicero connections")
        return conns
    except Exception as e:
        print(f"  [WARN] Failed to load Cicero connections: {e}")
        return None


def filter_connections_near_region(cicero_conns, region, window_kb=500):
    """Filter Cicero connections to those within +/-window_kb of a region."""
    chrom = region['chr']
    center = (int(region['start']) + int(region['end'])) // 2
    window = window_kb * 1000
    view_start = center - window
    view_end = center + window
    mask = (
        (cicero_conns['p1_chrom'] == chrom) &
        (cicero_conns['p2_chrom'] == chrom) &
        (cicero_conns['p1_start'] >= view_start) &
        (cicero_conns['p1_end'] <= view_end) &
        (cicero_conns['p2_start'] >= view_start) &
        (cicero_conns['p2_end'] <= view_end)
    )
    local = cicero_conns[mask].copy()
    return local if len(local) > 0 else None


# =============================================================================
# GLOBAL MODE — focal cell type, footprint at each gene TSS
# FIX-38: Complete rewrite for per-cell-type fan-out with TFBS scoring
# =============================================================================
def run_global_mode(args, printer, genome_obj, gene_regions, target_genes):
    """
    FIX-38: Focal cell-type footprinting with TFBS binding scores.

    Each Nextflow task receives ONE cell type. This function:
    1. Groups barcodes for the focal cell type only
    2. Loads TFBS + nucleosome binding score models (unless --skip-tfbs)
    3. For each target gene's promoter:
       a. Computes get_footprint_score (multiscale, scales 2-100)
       b. Computes get_binding_score (TFBS model, contextRadius=100)
       c. Generates Plot A: Tn5 insertion profile (plot_region_atac)
       d. Generates Plot B: Multiscale footprint heatmap (plot_footprints, stack=False)
       e. Generates Plot C: TFBS binding score visualization
       f. Generates composite 3-panel figure (A + B + C stacked vertically)
    """
    label = args.cell_type_label or "global"
    # FIX-43: Sanitize label for filesystem — cell types like "Megakaryocytes/platelets"
    # contain '/' which breaks filenames (interpreted as directory separator)
    safe_label = label.replace("/", "_")
    print(f"FIX-43: scPRINTER footprinting — focal cell type mode (label={label}, safe_label={safe_label})", flush=True)

    # Load peak matrix and build focal cell type grouping
    adata = ad.read_h5ad(args.peak_matrix)
    grouping, uniq_groups = build_focal_grouping(adata, printer, label,
                                                  cell_type_col=args.cell_type_col)

    # FIX-38: Load TFBS binding score model for integrated visualization
    has_tfbs = False
    if not args.skip_tfbs:
        try:
            print("  Loading TF binding score model (pretrained_TFBS_model)...")
            printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
            has_tfbs = True
            print("  [OK] TFBS model loaded")
        except Exception as e:
            print(f"  [WARN] Could not load TFBS model: {e}")
            print("  Continuing without TFBS binding scores")

    plots_dir = Path("plots")
    plots_dir.mkdir(exist_ok=True)
    promoter_results = []

    # FIX-39: Pre-initialize printer.uns tracking dicts as plain Python dicts.
    # The Rust anndata backend (anndata-0.4.2) caches HDF5 group handles
    # internally.  When get_footprint_score() writes a new entry to
    # printer.uns["footprints"] for gene 1, the HDF5 metadata cache becomes
    # stale.  On gene 2, the check `save_key in printer.uns["footprints"]`
    # triggers H5Gopen2() on the stale cache, causing a Rust panic:
    #   "called `Result::unwrap()` on an `Err` value: H5Gopen2(): unable to
    #    open group: incorrect cache entry type"
    # By converting these to plain Python dicts BEFORE the loop, all
    # subsequent `in` checks and mutations operate on Python objects, not
    # HDF5-backed groups.
    _init_uns_dict(printer, "footprints")
    _init_uns_dict(printer, "bindingscores")

    # Load TF→gene mapping if available (maps target genes back to regulating TFs)
    gene_to_tfs = {}
    if args.tf_mapping and os.path.exists(args.tf_mapping):
        try:
            with open(args.tf_mapping) as _f:
                tf_map_data = json.load(_f)
            ct_data = tf_map_data.get('cell_types', {}).get(label, {})
            for tf_name, tf_info in ct_data.items():
                for tgt in tf_info.get('target_genes', []):
                    if tgt not in gene_to_tfs:
                        gene_to_tfs[tgt] = []
                    gene_to_tfs[tgt].append(tf_name)
            print(f"  Loaded TF mapping: {len(gene_to_tfs)} genes with TF annotations")
        except Exception as e:
            print(f"  [WARN] Could not load TF mapping: {e}")

    # Load Cicero connections once (shared across all genes)
    all_cicero_conns = None
    if args.cicero_connections and os.path.exists(args.cicero_connections):
        all_cicero_conns = load_all_cicero_connections(args.cicero_connections)

    for gene in target_genes:
        if gene not in gene_regions:
            print(f"[WARN] Gene not found in coordinates: {gene}; skipping")
            continue

        region = gene_regions[gene]
        region_str = f"{region['chr']}:{region['start']}-{region['end']}"
        # TF annotation: which TFs are predicted to bind at this gene
        regulating_tfs = gene_to_tfs.get(gene, [])
        tf_label = f"TFs: {', '.join(regulating_tfs)}" if regulating_tfs else ""
        print(f"\n{'='*60}")
        print(f"[Global/Focal] Processing {gene}: {region_str}")
        if tf_label:
            print(f"  Regulating TFs: {', '.join(regulating_tfs)}")
        print(f"  Cell type: {label}, Cells: {len(grouping[0])}")
        print(f"{'='*60}")

        regions_df = pd.DataFrame([{
            "Chromosome": region["chr"],
            "Start": int(region["start"]),
            "End": int(region["end"]),
        }])

        # --- Step 1: Compute multiscale footprint scores ---
        fp_key = f"footprints_{safe_label}_{gene}"
        print(f"  Computing multiscale footprint scores (scales 2-100)...")
        scp.tl.get_footprint_score(
            printer,
            grouping,
            uniq_groups,
            regions_df,
            modes=np.arange(2, 101),
            n_jobs=args.cpus,
            save_key=fp_key,
            backed=False,
            overwrite=True,
        )

        fp_ok = hasattr(printer, "footprintsadata") and fp_key in printer.footprintsadata
        if fp_ok:
            fp_data = printer.footprintsadata[fp_key]
            fp_h5ad_name = f"footprints_{safe_label}_{gene}.h5ad"
            print(f"  Saving footprint data to {fp_h5ad_name}")
            try:
                fp_data = fp_data.copy()
            except Exception:
                pass
            fp_data.write(fp_h5ad_name)
        else:
            print(f"  [WARN] No footprintsadata entry for key '{fp_key}'")

        # --- Step 2: Compute TFBS binding scores ---
        bs_key = f"tfbs_{safe_label}_{gene}"
        bs_ok = False
        if has_tfbs:
            print(f"  Computing TFBS binding scores (contextRadius=100)...")
            try:
                scp.tl.get_binding_score(
                    printer,
                    grouping,
                    uniq_groups,
                    regions_df,
                    model_key="TF",
                    n_jobs=args.cpus,
                    contextRadius=100,
                    save_key=bs_key,
                    backed=False,
                    overwrite=True,
                )
                bs_ok = (
                    hasattr(printer, "bindingscoreadata")
                    and bs_key in printer.bindingscoreadata
                )
                if bs_ok:
                    bs_data = printer.bindingscoreadata[bs_key]
                    bs_h5ad_name = f"tfbs_{safe_label}_{gene}.h5ad"
                    print(f"  Saving TFBS data to {bs_h5ad_name}")
                    try:
                        bs_data = bs_data.copy()
                    except Exception:
                        pass
                    bs_data.write(bs_h5ad_name)
            except Exception as e:
                print(f"  [WARN] TFBS binding score failed for {gene}: {e}")

        # --- Step 3: Generate individual plots ---

        # Plot A: Tn5 insertion profile (raw chromatin accessibility)
        try:
            fig_a, ax_a = plt.subplots(1, 1, figsize=(12, 2.5))
            scp.pl.plot_region_atac(
                printer,
                grouping[0],
                regions_df.iloc[0],
                ax=ax_a,
                color="red",
            )
            title_a = f"{gene} — {label} — Tn5 Insertion Profile"
            if tf_label:
                title_a += f"\n{tf_label}"
            ax_a.set_title(title_a, fontsize=11)
            ax_a.set_ylabel("Tn5 Insertions", fontsize=9)
            _format_tss_xaxis(ax_a, region, gene)
            plt.tight_layout()
            figpath_a = plots_dir / f"{safe_label}_{gene}_A_insertion_profile.png"
            plt.savefig(figpath_a, dpi=200, bbox_inches="tight")
            print(f"  [OK] Plot A (insertion profile): {figpath_a}")
        except Exception as e:
            print(f"  [WARN] Plot A failed for {gene}: {e}")
            figpath_a = None
        plt.close("all")

        # Plot B: Multiscale footprint heatmap (THE key scprinter visualization)
        figpath_b = None
        if fp_ok:
            try:
                fig_b, ax_b = plt.subplots(1, 1, figsize=(12, 5))
                scp.pl.plot_footprints(
                    printer,
                    save_key=fp_key,
                    group_names=[label],
                    region=region_str,
                    ax=ax_b,
                    stack=False,
                    cmap="Blues",
                    vmin=0.5,
                    vmax=2.0,
                )
                title_b = f"{gene} — {label} — Multiscale Footprint"
                if tf_label:
                    title_b += f"\n{tf_label}"
                ax_b.set_title(title_b, fontsize=11)
                # FIX-38: Ensure y-axis has small scales (fine TF footprints) at
                # bottom and large scales (nucleosome-level) at top.
                ymin, ymax = ax_b.get_ylim()
                if ymin < ymax:
                    ax_b.invert_yaxis()
                # FIX-41: TSS-centered x-axis + scale y-axis labels
                _format_msfp_axes(ax_b, region, gene)
                plt.tight_layout()
                figpath_b = plots_dir / f"{safe_label}_{gene}_B_multiscale_footprint.png"
                plt.savefig(figpath_b, dpi=200, bbox_inches="tight")
                print(f"  [OK] Plot B (multiscale footprint): {figpath_b}")
            except Exception as e:
                print(f"  [WARN] Plot B failed for {gene}: {e}")
            plt.close("all")

        # Plot C: TFBS binding score — use scPrinter's built-in renderer
        figpath_c = None
        if bs_ok:
            try:
                fig_c, ax_c = plt.subplots(1, 1, figsize=(12, 3))
                scp.pl.plot_binding_score(
                    printer,
                    save_key=bs_key,
                    group_names=[label],
                    region=region_str,
                    ax=ax_c,
                    cmap="Greens",
                    vmin=0.0,
                    vmax=1.0,
                )
                title_c = f"{gene} — {label} — TF Binding Score"
                if tf_label:
                    title_c += f"\n{tf_label}"
                ax_c.set_title(title_c, fontsize=11)
                _format_tss_xaxis(ax_c, region, gene)
                ax_c.set_ylabel("TF Binding Score", fontsize=9)

                plt.tight_layout()
                figpath_c = plots_dir / f"{safe_label}_{gene}_C_tfbs_binding.png"
                plt.savefig(figpath_c, dpi=200, bbox_inches="tight")
                print(f"  [OK] Plot C (TFBS binding score): {figpath_c}")
            except Exception as e:
                import traceback as tb
                print(f"  [WARN] Plot C failed for {gene}: {e}\n{tb.format_exc()}")
            plt.close("all")

        # --- Step 4: Composite figure (A + B + C + D + E stacked vertically) ---
        # FIX-41: Always generate composite if at least footprint succeeded.
        # Panels: Insertion (A) -> MSFP heatmap (B) -> TFBS (C) -> Arcs (D) -> Logo (E).
        # All panels share TSS-centered x-axis; only bottom panel gets xlabel.
        if fp_ok:
            try:
                height_ratios = [1.0, 2.0]  # A: insertion, B: footprint
                n_panels = 2
                if bs_ok:
                    height_ratios.append(1.2)
                    n_panels += 1

                # Check for Cicero arcs panel
                has_arcs = False
                local_conns = None
                if all_cicero_conns is not None:
                    local_conns = filter_connections_near_region(all_cicero_conns, region, window_kb=500)
                    if local_conns is not None:
                        height_ratios.append(1.5)
                        n_panels += 1
                        has_arcs = True

                # Check for TF motif logo panel
                has_logo = False
                pwm_df = None
                primary_tf = regulating_tfs[0] if regulating_tfs else None
                if primary_tf and args.pfm_path and os.path.exists(args.pfm_path):
                    pwm_df = load_jaspar_pwm(args.pfm_path, primary_tf)
                    if pwm_df is not None:
                        height_ratios.append(0.8)
                        n_panels += 1
                        has_logo = True

                fig_comp, axes_comp = plt.subplots(
                    n_panels, 1,
                    figsize=(12, sum(height_ratios) * 2.5),
                    gridspec_kw={"height_ratios": height_ratios},
                )
                if n_panels == 1:
                    axes_comp = [axes_comp]
                elif n_panels > 1:
                    axes_comp = list(axes_comp)

                ax_idx = 0

                # Panel A: Insertion profile
                try:
                    scp.pl.plot_region_atac(
                        printer,
                        grouping[0],
                        regions_df.iloc[0],
                        ax=axes_comp[ax_idx],
                        color="red",
                    )
                    comp_title = f"{gene} — {label}"
                    if tf_label:
                        comp_title += f" | {tf_label}"
                    axes_comp[ax_idx].set_title(
                        comp_title, fontsize=12, fontweight="bold"
                    )
                    axes_comp[ax_idx].set_ylabel("Tn5 Insertions", fontsize=9)
                    # TSS-centered x ticks but no xlabel (not the bottom panel)
                    _format_tss_xaxis(axes_comp[ax_idx], region, gene)
                    axes_comp[ax_idx].set_xlabel("")
                except Exception as e:
                    print(f"  [WARN] Composite Panel A failed: {e}")
                ax_idx += 1

                # Panel B: Multiscale footprint heatmap
                try:
                    scp.pl.plot_footprints(
                        printer,
                        save_key=fp_key,
                        group_names=[label],
                        region=region_str,
                        ax=axes_comp[ax_idx],
                        stack=False,
                        cmap="Blues",
                        vmin=0.5,
                        vmax=2.0,
                        add_ticks=True,
                    )
                    # Fix y-axis orientation (small scales at bottom)
                    ymin, ymax = axes_comp[ax_idx].get_ylim()
                    if ymin < ymax:
                        axes_comp[ax_idx].invert_yaxis()
                    _format_msfp_axes(
                        axes_comp[ax_idx], region, gene,
                        xlabel=False, ylabel=True,
                    )
                except Exception as e:
                    print(f"  [WARN] Composite Panel B failed: {e}")
                ax_idx += 1

                # Panel C: TFBS binding score — use scPrinter's built-in renderer
                if bs_ok:
                    try:
                        scp.pl.plot_binding_score(
                            printer,
                            save_key=bs_key,
                            group_names=[label],
                            region=region_str,
                            ax=axes_comp[ax_idx],
                            cmap="Greens",
                            vmin=0.0,
                            vmax=1.0,
                        )
                        axes_comp[ax_idx].set_ylabel("TF Binding", fontsize=9)
                        _format_tss_xaxis(axes_comp[ax_idx], region, gene)
                    except Exception as e:
                        import traceback as tb
                        print(f"  [WARN] Composite Panel C failed: {e}\n{tb.format_exc()}")
                    ax_idx += 1

                # Panel D: Cicero CCAN arcs (NEW)
                if has_arcs:
                    try:
                        plot_cicero_arcs(
                            local_conns, region_str,
                            regions_df,  # pass regions_df for region marking
                            axes_comp[ax_idx], window_kb=500)
                        axes_comp[ax_idx].set_title(
                            f'D. Cicero CCANs at {gene} promoter (\u00b1500kb)',
                            fontsize=10, fontweight='bold', loc='left', pad=15)
                    except Exception as e:
                        print(f"  [WARN] Composite Panel D (arcs) failed: {e}")
                    ax_idx += 1

                # Panel E: TF motif logo (NEW)
                if has_logo:
                    try:
                        render_motif_logo(pwm_df, primary_tf, axes_comp[ax_idx])
                        axes_comp[ax_idx].set_title(
                            f'E. Binding motif: {primary_tf}',
                            fontsize=10, fontweight='bold', loc='left', pad=15)
                    except Exception as e:
                        print(f"  [WARN] Composite Panel E (logo) failed: {e}")
                    ax_idx += 1

                # Ensure only bottom panel has xlabel
                for i in range(ax_idx - 1):
                    axes_comp[i].set_xlabel("")

                plt.tight_layout()
                figpath_comp = plots_dir / f"{safe_label}_{gene}_composite.png"
                plt.savefig(figpath_comp, dpi=200, bbox_inches="tight")
                print(f"  [OK] Composite figure: {figpath_comp}")
            except Exception as e:
                print(f"  [WARN] Composite figure failed for {gene}: {e}")
            plt.close("all")

        promoter_results.append({
            "mode": "global",
            "label": label,
            "gene": gene,
            "region": region_str,
            "n_cells": len(grouping[0]),
            "has_footprint": fp_ok,
            "has_tfbs": bs_ok,
        })

    # Write summary
    if promoter_results:
        summary_df = pd.DataFrame(promoter_results)
        summary_df.to_csv("footprint_summary.csv", index=False)
        print(f"\nSaved global footprint summary ({len(promoter_results)} genes)")
    else:
        pd.DataFrame(columns=["mode", "gene", "region"]).to_csv(
            "footprint_summary.csv", index=False
        )
        print("No global footprint results produced")

    return promoter_results


# =============================================================================
# DIFFERENTIAL MODE — group by condition within cell type
# =============================================================================
def run_differential_mode(args, printer, genome_obj, gene_regions, target_genes):
    """
    Differential footprinting: groups by control/treatment condition
    within a single cell type.

    FIX-38: TFBS binding scores + per-condition multiscale heatmaps.
    FIX-42: All plots focus on the DIFFERENCE between conditions:
      1. Differential MSFP heatmap (treatment − control, RdBu_r) with TSS line
      2. Single-scale line plot at --diff-scale (overlay, same axes)
      3. Differential TFBS (treatment − control, fill_between)
    """
    label = args.cell_type_label or args.cell_type
    safe_label = label.replace("/", "_")
    print(f"FIX-38/42: scPRINTER footprinting — differential mode", flush=True)
    print(f"  Cell type:  {args.cell_type}")
    print(f"  Control:    {args.control_condition}")
    print(f"  Treatment:  {args.treatment_condition}")

    adata = ad.read_h5ad(args.peak_matrix)

    ct_col = args.cell_type_col
    if ct_col not in adata.obs.columns:
        raise ValueError(f"adata.obs does not contain a '{ct_col}' column. "
                         f"Available: {list(adata.obs.columns)}")

    cell_mask = adata.obs[ct_col] == args.cell_type
    n_cells = int(cell_mask.sum())
    print(f"  Found {n_cells} cells for cell_type = {args.cell_type}")
    if n_cells == 0:
        raise ValueError(f"No cells found for cell_type '{args.cell_type}'")

    adata_ct = adata[cell_mask].copy()

    # Accept either 'condition' or 'condition_group' column name
    cond_col = None
    for _cand in ["condition", "condition_group"]:
        if _cand in adata_ct.obs.columns:
            cond_col = _cand
            break
    if cond_col is None:
        raise ValueError("adata_ct.obs does not contain a 'condition' or 'condition_group' column. "
                         f"Available: {list(adata_ct.obs.columns)}")

    conditions = adata_ct.obs[cond_col].astype(str).unique()
    print(f"  Conditions found: {conditions}")

    control_label = args.control_condition
    treatment_label = args.treatment_condition

    group_labels = [lab for lab in [control_label, treatment_label] if lab in conditions]
    if len(group_labels) != 2:
        print(
            f"  SKIP: Cell type '{args.cell_type}' has only conditions {list(conditions)} "
            f"(need both {control_label} and {treatment_label}). "
            f"Writing empty outputs and returning."
        )
        # Write minimal outputs so Nextflow output globs are satisfied
        dummy = ad.AnnData()
        dummy.write_h5ad(f"footprints_{safe_label}_SKIPPED.h5ad")
        pd.DataFrame(columns=["mode", "gene", "region"]).to_csv(
            "footprint_summary.csv", index=False
        )
        return

    printer_barcodes = set(map(str, printer.obs_names))

    # FIX-42: Auto-detect barcode normalization mode
    all_bcs_raw = adata_ct.obs_names.astype(str)
    bc_norm_mode = choose_normalization(np.asarray(all_bcs_raw), printer_barcodes)

    grouping = []
    for lab in group_labels:
        mask = adata_ct.obs[cond_col].astype(str) == lab
        bcs_ct_raw = adata_ct.obs_names[mask].astype(str)
        bcs_ct = np.array(
            [apply_normalization(b, bc_norm_mode) for b in bcs_ct_raw], dtype=str
        )
        bcs_common = np.intersect1d(bcs_ct, np.array(list(printer_barcodes)))
        n_common = len(bcs_common)
        print(f"  Condition {lab}: {n_common} cells in printer (out of {int(mask.sum())} in peak matrix)")
        if n_common == 0:
            raise ValueError(f"No cells for condition '{lab}' in printer for cell type {args.cell_type}")
        grouping.append(list(bcs_common))

    uniq_groups = np.array(group_labels)

    # FIX-38: Load TFBS model for differential mode too
    has_tfbs = False
    if not args.skip_tfbs:
        try:
            print("  Loading TF binding score model (pretrained_TFBS_model)...")
            printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
            has_tfbs = True
            print("  [OK] TFBS model loaded")
        except Exception as e:
            print(f"  [WARN] Could not load TFBS model: {e}")

    plots_dir = Path("plots")
    plots_dir.mkdir(exist_ok=True)
    promoter_results = []

    # FIX-39: Pre-initialize uns tracking dicts (same as global mode)
    _init_uns_dict(printer, "footprints")
    _init_uns_dict(printer, "bindingscores")

    # Scale modes
    modes = np.arange(2, 101)

    for gene in target_genes:
        if gene not in gene_regions:
            print(f"Warning: Region not defined for {gene} -- skipping")
            continue

        region = gene_regions[gene]
        tss = _get_tss(region)
        region_str = f"{region['chr']}:{region['start']}-{region['end']}"
        print(f"\n{'='*60}")
        print(f"[Differential] Processing {gene}: {region_str} (TSS={tss})")
        print(f"  Cell type: {label}, Control: {control_label}, Treatment: {treatment_label}")
        print(f"{'='*60}")

        fp_key = f"footprints_{safe_label}_{gene}"
        bs_key = f"tfbs_{safe_label}_{gene}"

        regions_df = pd.DataFrame([{
            "Chromosome": region["chr"],
            "Start": int(region["start"]),
            "End": int(region["end"]),
        }])

        # Compute multiscale footprints
        scp.tl.get_footprint_score(
            printer,
            grouping,
            uniq_groups,
            regions_df,
            modes=modes,
            n_jobs=args.cpus,
            save_key=fp_key,
            backed=True,
            overwrite=True,
        )

        fp_ok = hasattr(printer, "footprintsadata") and fp_key in printer.footprintsadata
        if fp_ok:
            fp_data = printer.footprintsadata[fp_key]
            fp_h5ad_name = f"footprints_{safe_label}_{gene}.h5ad"
            print(f"  Saving promoter footprint data to {fp_h5ad_name}")
            try:
                fp_data = fp_data.copy()
            except Exception:
                pass
            fp_data.write(fp_h5ad_name)

        # Compute TFBS binding scores
        bs_ok = False
        if has_tfbs:
            try:
                scp.tl.get_binding_score(
                    printer,
                    grouping,
                    uniq_groups,
                    regions_df,
                    model_key="TF",
                    n_jobs=args.cpus,
                    contextRadius=100,
                    save_key=bs_key,
                    backed=True,
                    overwrite=True,
                )
                bs_ok = (
                    hasattr(printer, "bindingscoreadata")
                    and bs_key in printer.bindingscoreadata
                )
                if bs_ok:
                    bs_data = printer.bindingscoreadata[bs_key]
                    bs_h5ad_name = f"tfbs_{safe_label}_{gene}.h5ad"
                    print(f"  Saving TFBS data to {bs_h5ad_name}")
                    try:
                        bs_data = bs_data.copy()
                    except Exception:
                        pass
                    bs_data.write(bs_h5ad_name)
            except Exception as e:
                print(f"  [WARN] TFBS binding score failed for {gene}: {e}")

        # =====================================================================
        # PLOT 1: FIX-42 Differential heatmap (treatment - control, RdBu_r)
        # From SDas refactor: computes the difference matrix from
        # fp_data.obsm[region_key] and plots with diverging colormap + TSS line
        # =====================================================================
        if fp_ok:
            try:
                fp_data = printer.footprintsadata[fp_key]

                # Find the region key in obsm — scprinter stores footprint
                # matrices as fp_data.obsm[region_str] where region_str is
                # "chr:start-end"
                region_key = None
                if hasattr(fp_data, 'obsm') and fp_data.obsm is not None:
                    for k in fp_data.obsm_keys():
                        if region['chr'] in str(k):
                            region_key = k
                            break

                if region_key is not None:
                    M = np.array(fp_data.obsm[region_key])
                    # M shape: (n_groups, n_modes, n_positions) or (n_groups, n_modes * n_positions)
                    # obs_names should be [control_label, treatment_label]
                    obs_names = list(fp_data.obs_names.astype(str))
                    i_ctrl = obs_names.index(control_label) if control_label in obs_names else 0
                    i_trt = obs_names.index(treatment_label) if treatment_label in obs_names else 1

                    if M.ndim == 3:
                        diff_matrix = M[i_trt, :, :] - M[i_ctrl, :, :]
                    elif M.ndim == 2:
                        # Reshape: (n_groups, n_modes * n_positions) → (n_groups, n_modes, n_positions)
                        n_modes = len(modes)
                        n_pos = M.shape[1] // n_modes
                        M_3d = M.reshape(M.shape[0], n_modes, n_pos)
                        diff_matrix = M_3d[i_trt, :, :] - M_3d[i_ctrl, :, :]
                    else:
                        raise ValueError(f"Unexpected footprint matrix shape: {M.shape}")

                    # Symmetric color limits centered at 0
                    vmax_diff = np.nanpercentile(np.abs(diff_matrix), 95)
                    if vmax_diff == 0:
                        vmax_diff = 1.0

                    fig_diff, ax_diff = plt.subplots(1, 1, figsize=(12, 6))
                    im_diff = ax_diff.imshow(
                        diff_matrix,
                        aspect="auto",
                        cmap="RdBu_r",
                        vmin=-vmax_diff,
                        vmax=vmax_diff,
                        origin="lower",
                        extent=[
                            int(region["start"]),
                            int(region["end"]),
                            0,
                            diff_matrix.shape[0],
                        ],
                    )
                    plt.colorbar(im_diff, ax=ax_diff, label="Δ Footprint Score",
                                 shrink=0.7)

                    # FIX-42: Vertical line at TSS position
                    ax_diff.axvline(x=tss, color="black", linewidth=1.5,
                                    linestyle="--", alpha=0.7)

                    ax_diff.set_title(
                        f"{gene} — {label} — Differential Footprint "
                        f"({treatment_label} − {control_label})",
                        fontsize=11
                    )
                    _format_msfp_axes(ax_diff, region, gene, modes=modes)
                    plt.tight_layout()
                    diff_path = plots_dir / f"{safe_label}_{gene}_differential_heatmap.png"
                    plt.savefig(diff_path, dpi=200, bbox_inches="tight")
                    plt.close()
                    print(f"  [OK] Differential heatmap (RdBu_r): {diff_path}")
                else:
                    print(f"  [INFO] No obsm region key found for differential heatmap")
            except Exception as e:
                print(f"  [WARN] Differential heatmap failed for {gene}: {e}")
                plt.close("all")

        # =====================================================================
        # PLOT 2: FIX-42 Single-scale line plot (direct footprint comparison)
        # From SDas refactor: extracts a 1D cross-section at --diff-scale
        # to show direct footprint profile overlay for both conditions
        # =====================================================================
        if fp_ok:
            try:
                fp_data = printer.footprintsadata[fp_key]
                region_key = None
                if hasattr(fp_data, 'obsm') and fp_data.obsm is not None:
                    for k in fp_data.obsm_keys():
                        if region['chr'] in str(k):
                            region_key = k
                            break

                if region_key is not None:
                    M = np.array(fp_data.obsm[region_key])
                    obs_names = list(fp_data.obs_names.astype(str))
                    i_ctrl = obs_names.index(control_label) if control_label in obs_names else 0
                    i_trt = obs_names.index(treatment_label) if treatment_label in obs_names else 1

                    # Find the scale index closest to --diff-scale
                    target_scale = args.diff_scale
                    scale_idx = int(np.argmin(np.abs(modes - target_scale)))
                    actual_scale = modes[scale_idx]

                    if M.ndim == 3:
                        profile_ctrl = M[i_ctrl, scale_idx, :]
                        profile_trt = M[i_trt, scale_idx, :]
                    elif M.ndim == 2:
                        n_modes = len(modes)
                        n_pos = M.shape[1] // n_modes
                        M_3d = M.reshape(M.shape[0], n_modes, n_pos)
                        profile_ctrl = M_3d[i_ctrl, scale_idx, :]
                        profile_trt = M_3d[i_trt, scale_idx, :]
                    else:
                        raise ValueError(f"Unexpected shape: {M.shape}")

                    # X-axis: genomic positions relative to TSS
                    x_genomic = np.linspace(
                        int(region["start"]), int(region["end"]),
                        len(profile_ctrl)
                    )
                    x_tss_relative = x_genomic - tss

                    fig_line, ax_line = plt.subplots(1, 1, figsize=(12, 4))
                    ax_line.plot(x_tss_relative, profile_ctrl,
                                color="steelblue", linewidth=1.5, alpha=0.8,
                                label=control_label)
                    ax_line.plot(x_tss_relative, profile_trt,
                                color="indianred", linewidth=1.5, alpha=0.8,
                                label=treatment_label)
                    ax_line.axvline(x=0, color="black", linewidth=1,
                                   linestyle="--", alpha=0.5, label="TSS")
                    ax_line.legend(fontsize=9, loc="upper right")
                    ax_line.set_title(
                        f"{gene} — {label} — Footprint Profile at {actual_scale} bp scale",
                        fontsize=11
                    )
                    ax_line.set_xlabel(
                        f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                        fontsize=9
                    )
                    ax_line.set_ylabel("Footprint Score", fontsize=9)
                    plt.tight_layout()
                    line_path = plots_dir / f"{safe_label}_{gene}_differential_lineplot.png"
                    plt.savefig(line_path, dpi=200, bbox_inches="tight")
                    plt.close()
                    print(f"  [OK] Differential line plot (scale={actual_scale}bp): {line_path}")
            except Exception as e:
                print(f"  [WARN] Differential line plot failed for {gene}: {e}")
                plt.close("all")

        # =====================================================================
        # PLOT 3: Differential TFBS binding score (treatment − control)
        # Shows the difference in TF occupancy prediction between conditions.
        # Uses RdBu_r diverging colormap with TSS vertical line — consistent
        # with the differential footprint heatmap above.
        # =====================================================================
        if bs_ok:
            try:
                bs_adata = printer.bindingscoreadata[bs_key]
                import scipy.sparse as sp
                bs_matrix = bs_adata.X.toarray() if sp.issparse(bs_adata.X) else np.asarray(bs_adata.X)
                if bs_matrix.ndim == 1:
                    bs_matrix = bs_matrix.reshape(1, -1)

                obs_names_bs = list(bs_adata.obs_names.astype(str))
                i_ctrl_bs = obs_names_bs.index(control_label) if control_label in obs_names_bs else 0
                i_trt_bs = obs_names_bs.index(treatment_label) if treatment_label in obs_names_bs else 1

                if bs_matrix.shape[0] >= 2:
                    bs_diff = bs_matrix[i_trt_bs, :] - bs_matrix[i_ctrl_bs, :]
                    vmax_bs = np.nanpercentile(np.abs(bs_diff), 95)
                    if vmax_bs == 0:
                        vmax_bs = 1.0

                    fig_bs, ax_bs = plt.subplots(1, 1, figsize=(12, 3))
                    x_positions = np.linspace(
                        int(region["start"]), int(region["end"]),
                        len(bs_diff)
                    )
                    x_tss_rel = x_positions - tss

                    ax_bs.fill_between(
                        x_tss_rel, bs_diff, 0,
                        where=(bs_diff > 0), color="indianred", alpha=0.6,
                        label=f"↑ {treatment_label}"
                    )
                    ax_bs.fill_between(
                        x_tss_rel, bs_diff, 0,
                        where=(bs_diff < 0), color="steelblue", alpha=0.6,
                        label=f"↑ {control_label}"
                    )
                    ax_bs.axhline(y=0, color="gray", linewidth=0.5, linestyle="-")
                    ax_bs.axvline(x=0, color="black", linewidth=1,
                                  linestyle="--", alpha=0.5, label="TSS")
                    ax_bs.legend(fontsize=8, loc="upper right")
                    ax_bs.set_title(
                        f"{gene} — {label} — Δ TF Binding Score "
                        f"({treatment_label} − {control_label})",
                        fontsize=11
                    )
                    ax_bs.set_xlabel(
                        f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                        fontsize=9
                    )
                    ax_bs.set_ylabel("Δ TF Binding Score", fontsize=9)
                    plt.tight_layout()
                    plot_bs_path = plots_dir / f"{safe_label}_{gene}_tfbs_differential.png"
                    plt.savefig(plot_bs_path, dpi=200, bbox_inches="tight")
                    plt.close()
                    print(f"  [OK] Differential TFBS plot: {plot_bs_path}")
                else:
                    print(f"  [INFO] TFBS matrix has <2 groups; skipping differential TFBS plot")
            except Exception as e:
                print(f"  [WARN] Differential TFBS plot failed for {gene}: {e}")
                plt.close("all")

        promoter_results.append({
            "mode": "differential",
            "cell_type": args.cell_type,
            "gene": gene,
            "region": region_str,
            "tss": tss,
            "n_cells_control": len(grouping[0]),
            "n_cells_treatment": len(grouping[1]),
            "has_tfbs": bs_ok,
        })

    # Save promoter summary
    if promoter_results:
        summary_df = pd.DataFrame(promoter_results)
        summary_df.to_csv("footprint_summary.csv", index=False)
        print(f"\nSaved differential promoter summary ({len(promoter_results)} genes)")
    else:
        pd.DataFrame().to_csv("footprint_summary.csv", index=False)

    # --- Enhancer footprinting (differential only) ---
    # FIX-38 note: Enhancer footprinting stays as a stub for differential mode.
    # The integrated epigenome TF/enhancer machinery footprinter (incorporating
    # ChromVAR, enhancer mapping, Cicero CCANs, and expression data) will be
    # a distinct module.  The stub is preserved for when that module is ready.
    print("\n" + "=" * 60)
    print("Enhancer footprinting (differential mode — stub)")
    print("=" * 60)

    try:
        from run_scprinter_footprinting_enhancers import run_enhancer_analysis
        run_enhancer_analysis(
            printer, grouping, uniq_groups, gene_regions, target_genes,
            control_label, treatment_label, args
        )
    except ImportError:
        print("  Enhancer module not available; skipping enhancer footprinting")
        print("  (Future: integrated epigenome TF/enhancer module with ChromVAR,")
        print("   Cicero CCANs, enhancer mapping, and expression data)")
    except Exception as e:
        print(f"  Enhancer footprinting error: {e}")

    return promoter_results


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()

    # Configure scPrinter cache
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir

    # Initialize scprinter
    global scp
    scp = init_scprinter(args.cache_dir)

    # Resolve genome object
    try:
        genome_obj = getattr(scp.genome, args.genome)
    except AttributeError:
        raise ValueError(f"Genome '{args.genome}' not found in scprinter.genome")

    # Load printer
    printer_path = Path(args.printer_path)
    if not printer_path.exists():
        raise FileNotFoundError(f"Printer file not found: {printer_path}")

    print(f"Loading printer from {printer_path}")
    printer = scp.load_printer(str(printer_path), genome_obj)
    printer.load_disp_model()

    # Parse target genes
    target_genes = [g.strip() for g in args.target_genes.split(",") if g.strip()]
    print(f"Target genes: {target_genes}")

    # Load gene coordinates
    if args.coordinates_file:
        print(f"Loading gene coordinates from {args.coordinates_file}")
        gene_regions = load_gene_coordinates(args.coordinates_file)
        missing = [g for g in target_genes if g not in gene_regions]
        if missing:
            print(f"WARNING: No coordinates for: {missing}")
            target_genes = [g for g in target_genes if g in gene_regions]
            if not target_genes:
                raise ValueError("No target genes have resolved coordinates")
    else:
        raise ValueError("--coordinates-file is required (use gene_coordinate_resolver.py output)")

    # =========================================================================
    # FIX-23/38: Mode detection
    # =========================================================================
    has_control = args.control_condition and args.control_condition.strip()
    has_treatment = args.treatment_condition and args.treatment_condition.strip()

    if has_control and has_treatment:
        run_differential_mode(args, printer, genome_obj, gene_regions, target_genes)
    elif not has_control and not has_treatment:
        run_global_mode(args, printer, genome_obj, gene_regions, target_genes)
    else:
        raise ValueError(
            "Must provide BOTH --control-condition and --treatment-condition (differential mode) "
            "or NEITHER (global mode). Got only one."
        )

    # Clean up
    if hasattr(printer, "close"):
        printer.close()

    print("\nFIX-38: scPRINTER footprinting complete")


if __name__ == "__main__":
    main()
