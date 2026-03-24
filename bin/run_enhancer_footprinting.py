#!/usr/bin/env python3
"""
run_enhancer_footprinting.py — Recipe Steps A9+A10

Run scPRINTER binding scores + multiscale footprints on a single (cell_type, TF)
enhancer region set.  Follows the same pattern as run_scprinter_footprinting.py
(FIX-23b) but operates on enhancer BED regions instead of gene promoters.

Inputs:
    - Per-TF BED file (enhancer regions for one TF, 4-col: chr, start, end, linked_genes)
    - scPrinter printer h5ad
    - Peak matrix h5ad (for barcode groupings)

Outputs:
    - enhancer_footprints_{ct}_{tf}.h5ad
    - enhancer_tfbs_{ct}_{tf}.h5ad
    - enhancer_plots/*.png
    - enhancer_fp_summary.csv
"""

import argparse
import builtins
import os
from pathlib import Path
from typing import Dict

# FIX-75: Ensure HDF5 file locking is disabled before any h5py/anndata import
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# FIX-44: scPRINTER's get_binding_score() uses Python 2's `long` type internally
builtins.long = int

import scprinter as scp
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def init_scprinter(cache_dir):
    """Set scprinter dataset cache path."""
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp


def _get_tss(region):
    """Get TSS position from region dict (copied from run_scprinter_footprinting.py)."""
    if "tss" in region:
        return int(region["tss"])
    return (int(region["start"]) + int(region["end"])) // 2


def _format_tfbs_axes(ax, region: Dict, gene: str, is_enhancer=False):
    """FIX-87: Reformat binding score x-axis from bin indices to genomic coordinates.

    scp.pl.plot_binding_score() uses internal bin indices (0-30).
    This relabels to match the MSFP genomic coordinate system.
    """
    start = int(region["start"])
    end = int(region["end"])
    center = (start + end) // 2 if is_enhancer else _get_tss(region)

    try:
        xlim = ax.get_xlim()
        n_ticks = 7
        tick_positions = np.linspace(xlim[0], xlim[1], n_ticks)
        tick_genomic = start + (tick_positions - xlim[0]) / (xlim[1] - xlim[0]) * (end - start)
        tick_labels = [f"{int(g - center):+d}" for g in tick_genomic]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        if is_enhancer:
            ax.set_xlabel(f"Position in enhancer region (bp) — {region['chr']}:{start}-{end}",
                          fontsize=9)
        else:
            ax.set_xlabel(f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                          fontsize=9)
        # Add center line
        center_pixel = xlim[0] + (center - start) / (end - start) * (xlim[1] - xlim[0])
        if xlim[0] <= center_pixel <= xlim[1]:
            ax.axvline(center_pixel, color='black', linestyle='--', linewidth=0.8, alpha=0.7)
    except Exception as e:
        print(f"  [WARN] _format_tfbs_axes failed: {e}")


def _format_msfp_axes(ax, region: Dict, gene: str, modes=None,
                      xlabel=True, ylabel=True, is_enhancer=False):
    """Post-process multiscale footprint heatmap axes (from run_scprinter_footprinting.py).

    Relabels raw array indices to publication-quality labels:
      - X-axis: distance from region center in bp (enhancer) or TSS (promoter)
      - Y-axis: Footprint scale in bp (2, 10, 20, 50, 100)
      - Center vertical line at midpoint (enhancer) or TSS (promoter)
    """
    if modes is None:
        modes = np.arange(2, 101)
    start = int(region["start"])
    end = int(region["end"])
    center = (start + end) // 2 if is_enhancer else _get_tss(region)

    try:
        xlim = ax.get_xlim()
        n_ticks = 7
        tick_positions = np.linspace(xlim[0], xlim[1], n_ticks)
        tick_genomic = start + (tick_positions - xlim[0]) / (xlim[1] - xlim[0]) * (end - start)
        tick_labels = [f"{int(g - center):+d}" for g in tick_genomic]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        if xlabel:
            if is_enhancer:
                ax.set_xlabel(f"Position in enhancer region (bp) — {region['chr']}:{start}-{end}",
                              fontsize=9)
            else:
                ax.set_xlabel(f"Position relative to TSS (bp) — {gene} ({region['chr']})",
                              fontsize=9)
        center_pixel = xlim[0] + (center - start) / (end - start) * (xlim[1] - xlim[0])
        if xlim[0] <= center_pixel <= xlim[1]:
            label = 'center' if is_enhancer else 'TSS'
            ax.axvline(center_pixel, color='black', linestyle='--', linewidth=0.8, alpha=0.7,
                       label=label)
    except Exception as e:
        print(f"  [WARN] _format_msfp_axes x-axis formatting failed: {e}")

    try:
        y_lo, y_hi = ax.get_ylim()
        scale_ticks = [2, 10, 20, 50, 100]
        n_modes = len(modes)
        tick_positions_y = []
        tick_labels_y = []
        for s in scale_ticks:
            if s in modes:
                idx = np.where(modes == s)[0][0]
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


def normalize_peak_barcode(bc, strategy="strip"):
    """Normalize peak-matrix barcode to match printer.obs_names.

    Strategies:
      'strip'   — 'sample:barcode' → 'barcode' (10x: strip prefix)
      'replace' — 'sample:barcode' → 'sample_barcode' (BD: colon→underscore)
      'identity' — no transformation
    """
    if strategy == "strip" and ":" in bc:
        return bc.split(":", 1)[1]
    elif strategy == "replace" and ":" in bc:
        return bc.replace(":", "_", 1)
    return bc


def detect_barcode_strategy(peak_barcodes, printer_barcodes, n_probe=200):
    """Auto-detect which normalization strategy aligns peak matrix to printer.

    Tests 'strip' (10x), 'replace' (BD colon→underscore), and 'identity'.
    Returns the strategy with the most barcode overlap.
    """
    sample = list(peak_barcodes)[:n_probe]
    printer_set = set(printer_barcodes)
    best_strategy, best_count = "identity", 0
    for strategy in ("strip", "replace", "identity"):
        normed = [normalize_peak_barcode(b, strategy) for b in sample]
        count = sum(1 for b in normed if b in printer_set)
        if count > best_count:
            best_strategy, best_count = strategy, count
    return best_strategy


def parse_args():
    p = argparse.ArgumentParser(description="Enhancer footprinting for a single (cell_type, TF) pair")
    p.add_argument("--region-set", required=True, help="BED file with enhancer regions")
    p.add_argument("--printer-path", required=True, help="scPrinter printer h5ad")
    p.add_argument("--peak-matrix", required=True, help="Filtered peak matrix h5ad")
    p.add_argument("--cell-type", required=True, help="Cell type label")
    p.add_argument("--tf-name", required=True, help="TF name label")
    p.add_argument("--cache-dir", required=True, help="scPrinter cache dir")
    p.add_argument("--genome", default="hg38", help="Genome key")
    p.add_argument("--pfm-path", default="", help="Path to JASPAR PFM file for motif logo rendering")
    p.add_argument("--gtf", default="", help="GTF annotation file for TSS lookup")
    p.add_argument("--cicero-connections", default="", help="Cicero connections TSV (gzipped)")
    p.add_argument("--cpus", type=int, default=4)
    return p.parse_args()


def select_best_region(printer, save_key_bs, regions_df, regions_full,
                       uniq_groups, ct, has_tfbs, target_chrom=None):
    """Select the best region for the insertion plot.

    Strategy: use TF binding scores to find the region with the strongest
    signal in the target cell type.  Falls back to the middle region if
    binding scores are unavailable.

    If target_chrom is provided, restrict selection to regions on that
    chromosome (cis-regulation filter, FIX-82).
    """
    # FIX-82: Filter to same chromosome as target gene for cis-regulation
    if target_chrom:
        cis_mask = regions_df['Chromosome'] == target_chrom
        cis_indices = np.where(cis_mask)[0]
        if len(cis_indices) > 0:
            print(f"  FIX-82: Filtering to {len(cis_indices)}/{len(regions_df)} "
                  f"cis-regulatory regions on {target_chrom}")
        else:
            print(f"  FIX-82: No regions on {target_chrom}, using all regions")
            cis_indices = np.arange(len(regions_df))
    else:
        cis_indices = np.arange(len(regions_df))

    # Default: middle of eligible region list
    best_idx = cis_indices[len(cis_indices) // 2]

    if not has_tfbs:
        return best_idx

    try:
        if not (hasattr(printer, 'bindingscoreadata') and
                save_key_bs in printer.bindingscoreadata):
            return best_idx

        bs_data = printer.bindingscoreadata[save_key_bs]

        # Find the row for our target cell type
        ct_pos = None
        for i, g in enumerate(uniq_groups):
            if g == ct:
                ct_pos = i
                break

        if ct_pos is None or ct_pos >= bs_data.X.shape[0]:
            return best_idx

        # FIX-52: pyanndata 0.4.1 Rust backend does not support slice
        # indexing on .X (panics with "not yet implemented" in slice.rs).
        # Materialize the entire matrix to numpy/scipy first, then index.
        import scipy.sparse as sp
        X_raw = bs_data.X
        if sp.issparse(X_raw):
            X_mat = X_raw.toarray()
        else:
            X_mat = np.asarray(X_raw)
        if X_mat.ndim == 1:
            X_mat = X_mat.reshape(1, -1)
        row = X_mat[ct_pos].flatten()
        scores = row

        # Clamp to actual region count
        scores = scores[:len(regions_df)]

        if len(scores) > 0:
            # Restrict to cis-regulatory regions
            cis_scores = scores[cis_indices]
            valid = ~np.isnan(cis_scores)
            if valid.any():
                best_within_cis = int(np.argmax(np.where(valid, cis_scores, -np.inf)))
                best_idx = cis_indices[best_within_cis]
                print(f"  Selected cis-region {best_idx} by binding score "
                      f"({scores[best_idx]:.3f})")

    except Exception as e:
        print(f"  NOTE: Binding score region selection failed: {e}")

    return best_idx


def format_gene_label(genes_str, max_genes=5):
    """Format a comma-separated gene string into a concise label."""
    if not genes_str or genes_str == 'nan':
        return ''
    gene_list = [g.strip() for g in genes_str.split(',') if g.strip()]
    if not gene_list:
        return ''
    if len(gene_list) <= max_genes:
        return ', '.join(gene_list)
    return ', '.join(gene_list[:max_genes]) + f' (+{len(gene_list) - max_genes} more)'


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
                                p = (counts[nuc][pos] + 0.5) / total
                                row[nuc] = p * np.log2(p / 0.25 + 1e-10)
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


def load_protein_coding_genes(gtf_path):
    """FIX-85: Parse GTF and return set of protein-coding gene names (uppercase).

    Looks for gene_biotype or gene_type attribute == 'protein_coding'.
    """
    import gzip
    pc_genes = set()
    opener = gzip.open if gtf_path.endswith('.gz') else open
    try:
        with opener(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9 or parts[2] != 'gene':
                    continue
                attrs = parts[8]
                gene_name = None
                is_pc = False
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    if attr.startswith(('gene_biotype', 'gene_type')):
                        val = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                        if val == 'protein_coding':
                            is_pc = True
                if gene_name and is_pc:
                    pc_genes.add(gene_name.upper())
    except Exception as e:
        print(f"  [WARN] GTF protein-coding parsing failed: {e}")
    return pc_genes


def parse_tss_from_gtf(gtf_path, tf_name):
    """Look up the TSS for a TF gene from a GTF file.

    Returns (chrom, tss_pos, strand) or None if not found.
    Handles composite TF names like 'FOS::JUND' by trying the first component.
    """
    import gzip

    candidates = [tf_name]
    if '::' in tf_name:
        candidates.extend(tf_name.split('::'))
    candidates_upper = {c.upper() for c in candidates}

    opener = gzip.open if gtf_path.endswith('.gz') else open

    best = None
    try:
        with opener(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9 or parts[2] != 'gene':
                    continue
                attrs = parts[8]
                gene_name = None
                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                        break
                if gene_name and gene_name.upper() in candidates_upper:
                    chrom = parts[0]
                    strand = parts[6]
                    tss = int(parts[3]) if strand == '+' else int(parts[4])
                    best = (chrom, tss, strand)
                    # Prefer exact match on original tf_name
                    if gene_name.upper() == tf_name.upper():
                        return best
    except Exception as e:
        print(f"  [WARN] GTF parsing failed: {e}")

    return best


def load_cicero_connections(conns_path, regions_df):
    """Load Cicero connections and filter to those involving our enhancer regions.

    Returns a DataFrame with columns: Peak1, Peak2, coaccess, plus parsed coordinates.
    """
    import gzip

    if not conns_path or not os.path.isfile(conns_path):
        return None

    try:
        opener = gzip.open if conns_path.endswith('.gz') else open
        conns = pd.read_csv(conns_path, sep='\t')

        # Build a set of region strings for fast lookup
        region_strs = set()
        for _, row in regions_df.iterrows():
            region_strs.add(f"{row['Chromosome']}:{row['Start']}-{row['End']}")

        # Filter to connections where at least one peak is in our enhancer set
        mask = conns['Peak1'].isin(region_strs) | conns['Peak2'].isin(region_strs)
        filtered = conns[mask].copy()

        if filtered.empty:
            return None

        # Parse coordinates for both peaks
        def parse_peak(p):
            chrom, rest = p.split(':')
            start, end = rest.split('-')
            return chrom, int(start), int(end)

        for col, prefix in [('Peak1', 'p1_'), ('Peak2', 'p2_')]:
            parsed = filtered[col].apply(parse_peak)
            filtered[f'{prefix}chrom'] = [p[0] for p in parsed]
            filtered[f'{prefix}start'] = [p[1] for p in parsed]
            filtered[f'{prefix}end'] = [p[2] for p in parsed]

        return filtered.sort_values('coaccess', ascending=False).reset_index(drop=True)

    except Exception as e:
        print(f"  [WARN] Cicero connections loading failed: {e}")
        return None


def plot_coaccessibility_heatmap(cicero_conns, regions_full, ax, tf_name, top_n=20):
    """Plot a co-accessibility heatmap of ALL enhancer regions.

    Each row is one enhancer region from the BED. Color = max Cicero
    coaccessibility score involving that region (0 if no Cicero signal).
    Regions WITH Cicero signal are highlighted; those without remain dark.
    """
    if regions_full is None or regions_full.empty:
        ax.text(0.5, 0.5, 'No enhancer regions',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('Enhancer Co-accessibility', fontsize=10)
        ax.axis('off')
        return

    # Build region key for each enhancer
    region_keys = []
    for _, row in regions_full.iterrows():
        region_keys.append(f"{row['Chromosome']}:{row['Start']}-{row['End']}")

    # Compute max coaccessibility per enhancer region
    region_max_coaccess = {}
    if cicero_conns is not None and not cicero_conns.empty:
        for _, conn in cicero_conns.iterrows():
            p1, p2 = conn['Peak1'], conn['Peak2']
            score = conn['coaccess']
            for pk in [p1, p2]:
                if pk in region_max_coaccess:
                    region_max_coaccess[pk] = max(region_max_coaccess[pk], score)
                else:
                    region_max_coaccess[pk] = score

    # Build score array for all regions (0 if no Cicero signal)
    scores = [region_max_coaccess.get(rk, 0.0) for rk in region_keys]
    n_with_signal = sum(1 for s in scores if s > 0)

    # Build labels (short coords + linked genes)
    labels = []
    for i, rk in enumerate(region_keys):
        chrom, rest = rk.split(':')
        label = f"{chrom}:{rest}"
        if 'linked_genes' in regions_full.columns:
            genes = str(regions_full.iloc[i].get('linked_genes', ''))
            if genes and genes != 'nan':
                gene_list = [g.strip() for g in genes.split(',')][:2]
                label += f" ({', '.join(gene_list)})"
        labels.append(label)

    # Render heatmap
    scores_array = np.array(scores).reshape(-1, 1)
    vmax = max(scores) if max(scores) > 0 else 1
    im = ax.imshow(scores_array, aspect='auto', cmap='YlOrRd', vmin=0, vmax=vmax)

    # Y-axis labels: show all if <=40 regions, otherwise show every Nth
    n_regions = len(labels)
    if n_regions <= 40:
        ax.set_yticks(range(n_regions))
        ax.set_yticklabels(labels, fontsize=5)
    else:
        step = max(1, n_regions // 20)
        tick_pos = list(range(0, n_regions, step))
        ax.set_yticks(tick_pos)
        ax.set_yticklabels([labels[i] for i in tick_pos], fontsize=5)

    ax.set_xticks([])
    ax.set_title(
        f'Enhancer co-accessibility ({n_with_signal}/{n_regions} regions with Cicero signal)',
        fontsize=10)
    plt.colorbar(im, ax=ax, label='Max co-accessibility', shrink=0.8)


def plot_cicero_arcs(cicero_conns, best_region_str, regions_df, ax, window_kb=500):
    """Plot Cicero connection arcs centered on the best enhancer region.

    Shows co-accessibility arcs within ±window_kb of the best region.
    """
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
        ax.text(0.5, 0.5, f'No Cicero connections within ±{window_kb}kb',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'CCAN arcs at {best_region_str}', fontsize=10)
        ax.axis('off')
        return

    # Build set of enhancer midpoints for highlighting
    enh_mids = set()
    for _, row in regions_df.iterrows():
        if row['Chromosome'] == chrom:
            mid = (row['Start'] + row['End']) // 2
            if view_start <= mid <= view_end:
                enh_mids.add(mid)

    # Normalize co-accessibility for color mapping
    max_coaccess = local['coaccess'].max()
    min_coaccess = local['coaccess'].min()

    from matplotlib.patches import Arc as MplArc

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

    # Mark enhancer regions as ticks on the baseline
    for mid in enh_mids:
        ax.plot(mid, 0, '|', color='red', markersize=6, markeredgewidth=1.5)

    # Mark the best region
    ax.axvspan(center_start, center_end, alpha=0.15, color='red', zorder=0)

    ax.set_xlim(view_start, view_end)
    max_arc_height = max(abs(r['p2_start'] - r['p1_start']) * 0.2 for _, r in local.iterrows()) if len(local) > 0 else window * 0.1
    ax.set_ylim(-max_arc_height * 0.1, max_arc_height * 1.2)
    ax.set_xlabel(f'{chrom} position', fontsize=8)
    ax.set_title(f'CCAN arcs at {best_region_str} (±{window_kb}kb, {len(local)} connections)', fontsize=10)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_yticks([])

    # FIX-82: Add colorbar legend for co-accessibility score
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.RdYlBu_r,
        norm=mcolors.Normalize(vmin=min_coaccess, vmax=max_coaccess))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, aspect=15)
    cbar.set_label('Cicero co-accessibility', fontsize=8)
    cbar.ax.tick_params(labelsize=7)


def main():
    args = parse_args()
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir

    global scp
    scp = init_scprinter(args.cache_dir)

    plots_dir = Path("enhancer_plots")
    plots_dir.mkdir(exist_ok=True)

    ct = args.cell_type
    tf = args.tf_name

    print(f"Enhancer footprinting: cell_type={ct}, TF={tf}")

    # Load regions (4-col BED: chr, start, end, linked_genes)
    regions_raw = pd.read_csv(args.region_set, sep='\t', header=None)
    if regions_raw.shape[1] >= 4:
        regions_raw.columns = ['Chromosome', 'Start', 'End', 'linked_genes'] + \
            [f'col{i}' for i in range(4, regions_raw.shape[1])]
        # Collect all unique linked genes across regions
        all_genes = set()
        for g in regions_raw['linked_genes'].dropna().astype(str):
            if g and g != 'nan':
                all_genes.update(gene.strip() for gene in g.split(',') if gene.strip())
        linked_genes_str = ', '.join(sorted(all_genes))
        n_linked_genes = len(all_genes)
    else:
        regions_raw.columns = ['Chromosome', 'Start', 'End']
        linked_genes_str = ''
        n_linked_genes = 0

    # scPRINTER only needs 3-column regions
    regions_df = regions_raw[['Chromosome', 'Start', 'End']].copy()

    # FIX-86: Pad regions by 250bp on each side for ±500bp viewer window
    REGION_PAD = 250
    regions_df['Start'] = (regions_df['Start'] - REGION_PAD).clip(lower=0)
    regions_df['End'] = regions_df['End'] + REGION_PAD
    # Also pad the full dataframe so indices stay aligned
    regions_raw['Start'] = (regions_raw['Start'] - REGION_PAD).clip(lower=0)
    regions_raw['End'] = regions_raw['End'] + REGION_PAD

    # Keep full data for per-region gene lookup
    regions_full = regions_raw.copy()
    print(f"  Loaded {len(regions_df)} enhancer regions (padded ±{REGION_PAD}bp), "
          f"linked to {n_linked_genes} genes")

    if regions_df.empty:
        print("WARNING: Empty region set, writing empty summary")
        pd.DataFrame(columns=['cell_type', 'tf', 'n_regions']).to_csv(
            'enhancer_fp_summary.csv', index=False)
        return

    # Load printer
    print(f"  Loading printer from {args.printer_path}")
    genome_obj = getattr(scp.genome, args.genome)
    printer = scp.load_printer(args.printer_path, genome_obj)
    printer.load_disp_model()

    # Load peak matrix for barcode groupings
    adata = ad.read_h5ad(args.peak_matrix)
    if 'cell_type' not in adata.obs.columns:
        raise ValueError("Peak matrix missing 'cell_type' column")

    # Build barcode -> cell_type table with auto-detected normalization
    printer_bcs = set(map(str, printer.obs_names))
    peak_bcs = adata.obs_names.astype(str).tolist()
    strategy = detect_barcode_strategy(peak_bcs, printer_bcs)
    print(f"  Barcode normalization strategy: {strategy}")

    df = pd.DataFrame({
        'barcode': [normalize_peak_barcode(b, strategy) for b in peak_bcs],
        'cell_type': adata.obs['cell_type'].astype(str).values,
    })

    df = df[df['barcode'].isin(printer_bcs)].copy()

    if df.empty:
        raise ValueError(
            f"No barcodes from peak matrix found in printer (strategy={strategy}). "
            f"Peak examples: {peak_bcs[:3]}, Printer examples: {list(printer_bcs)[:3]}"
        )

    print(f"  Matched {len(df)} barcodes across {df['cell_type'].nunique()} cell types")

    # Create per-cell-type groupings
    barcodeGroups = df[['barcode', 'cell_type']].copy()
    grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcodeGroups)

    order = np.argsort(uniq_groups)
    grouping = [grouping[i] for i in order]
    uniq_groups = uniq_groups[order]

    # Load TFBS binding score model
    has_tfbs = False
    try:
        print("  Loading TF binding score model...")
        printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)
        has_tfbs = True
        print("  [OK] TFBS model loaded")
    except Exception as e:
        print(f"  [WARN] Could not load TFBS model: {e}")

    # Compute binding scores
    save_key_bs = f"enhancer_bs_{ct}_{tf}"
    print(f"  Computing binding scores (key={save_key_bs})...")
    try:
        # FIX-99: Added model_key="TF" and contextRadius=100 to match
        # run_scprinter_footprinting.py (ATAC-only). Without these, binding
        # scores produce 0-dimensional arrays that can't be plotted.
        scp.tl.get_binding_score(
            printer,
            grouping,
            uniq_groups,
            regions_df,
            model_key="TF",
            n_jobs=args.cpus,
            contextRadius=100,
            save_key=save_key_bs,
            backed=False,
            overwrite=True,
        )

        if hasattr(printer, 'bindingscoreadata') and save_key_bs in printer.bindingscoreadata:
            bs_data = printer.bindingscoreadata[save_key_bs]
            bs_h5ad = f"enhancer_tfbs_{ct}_{tf}.h5ad"
            try:
                bs_data = bs_data.copy()
            except Exception:
                pass
            bs_data.write(bs_h5ad)
            print(f"  Saved binding scores to {bs_h5ad}")
    except Exception as e:
        print(f"  WARNING: Binding score computation failed: {e}")

    # Compute multiscale footprints
    save_key_fp = f"enhancer_fp_{ct}_{tf}"
    print(f"  Computing multiscale footprints (key={save_key_fp})...")
    scp.tl.get_footprint_score(
        printer,
        grouping,
        uniq_groups,
        regions_df,
        modes=np.arange(2, 101),
        n_jobs=args.cpus,
        save_key=save_key_fp,
        backed=False,
        overwrite=True,
    )

    fp_saved = False
    if hasattr(printer, 'footprintsadata') and save_key_fp in printer.footprintsadata:
        fp_data = printer.footprintsadata[save_key_fp]
        fp_h5ad = f"enhancer_footprints_{ct}_{tf}.h5ad"
        try:
            fp_data = fp_data.copy()
        except Exception:
            pass
        fp_data.write(fp_h5ad)
        fp_saved = True
        print(f"  Saved footprints to {fp_h5ad}")

        # ---- Plot 1: Full multiscale footprint heatmap (FIX-68a, FIX-72, FIX-82, FIX-85) ----
        # FIX-85: Load protein-coding gene set for target gene filtering
        pc_genes = load_protein_coding_genes(args.gtf) if args.gtf else set()
        if pc_genes:
            print(f"  FIX-85: Loaded {len(pc_genes)} protein-coding genes from GTF")

        # FIX-82/85: Identify target gene chromosome for cis-regulation filter
        target_gene = None
        target_chrom = None
        if n_linked_genes > 0 and 'linked_genes' in regions_full.columns:
            from collections import Counter
            gene_counts = Counter()
            for g in regions_full['linked_genes'].dropna().astype(str):
                if g and g != 'nan':
                    for gene in g.split(','):
                        gene = gene.strip()
                        if gene and gene.upper() != tf.upper():
                            gene_counts[gene] += 1
            # FIX-85: Filter to protein-coding genes only
            if pc_genes:
                pc_counts = {g: c for g, c in gene_counts.items()
                             if g.upper() in pc_genes}
                if pc_counts:
                    print(f"  FIX-85: {len(pc_counts)}/{len(gene_counts)} linked genes "
                          f"are protein-coding")
                    gene_counts = Counter(pc_counts)
                else:
                    print(f"  FIX-85: No protein-coding linked genes found, "
                          f"using all {len(gene_counts)} genes")
            if gene_counts:
                target_gene = gene_counts.most_common(1)[0][0]
                tss_info = parse_tss_from_gtf(args.gtf, target_gene) if args.gtf else None
                if tss_info:
                    target_chrom = tss_info[0]
                    print(f"  Target gene: {target_gene} on {target_chrom} "
                          f"(linked to {gene_counts[target_gene]} enhancer regions)")

        fig1_path = None
        try:
            # FIX-82: Select best cis-regulatory region (same chrom as target gene)
            best_msfp_idx = select_best_region(
                printer, save_key_bs, regions_df, regions_full,
                uniq_groups, ct, has_tfbs, target_chrom=target_chrom)
            ref_row = regions_df.iloc[best_msfp_idx]
            ref_region = f"{ref_row['Chromosome']}:{ref_row['Start']}-{ref_row['End']}"

            # Get linked genes for this specific region
            ref_genes_str = ''
            if 'linked_genes' in regions_full.columns:
                ref_genes_str = str(regions_full.iloc[best_msfp_idx].get('linked_genes', ''))
                if ref_genes_str == 'nan':
                    ref_genes_str = ''

            fig = plt.figure(figsize=(10, 6))
            ax_msfp = fig.add_subplot(111)
            scp.pl.plot_footprints(
                printer,
                save_key=save_key_fp,
                group_names=[ct],
                region=ref_region,
                ax=ax_msfp,
                stack=False,
                cmap="Blues",
                vmin=0.5,
                vmax=2.0,
                add_ticks=True,
            )
            # Invert y-axis so small scales are at bottom
            ymin, ymax = ax_msfp.get_ylim()
            if ymin < ymax:
                ax_msfp.invert_yaxis()
            # FIX-72: Label as enhancer region, not TSS
            region_dict = {
                'chr': ref_row['Chromosome'],
                'start': int(ref_row['Start']),
                'end': int(ref_row['End']),
            }
            _format_msfp_axes(ax_msfp, region_dict, tf, is_enhancer=True)

            title = (f"{ct} — TF {tf} multiscale footprint at {ref_region}"
                     f" ({len(regions_df)} total enhancer regions)")
            if ref_genes_str:
                title += f"\nLinked genes: {format_gene_label(ref_genes_str)}"
            plt.suptitle(title, fontsize=10)
            plt.tight_layout()
            fig1_path = plots_dir / f"{ct}_{tf}_enhancer_footprints_msfp.png"
            plt.savefig(fig1_path, dpi=200)
            plt.close()
            print(f"  Saved plot: {fig1_path}")
        except Exception as e:
            print(f"  WARNING: Footprint plot failed: {e}")
            plt.close('all')

        # ---- Plot 2: Target gene TSS insertion + aggregate enhancer profile ----
        fig2_path = None
        ct_idx = np.where(uniq_groups == ct)[0]
        if len(ct_idx) > 0 and len(regions_df) > 0:
            try:
                best_idx = select_best_region(
                    printer, save_key_bs, regions_df, regions_full,
                    uniq_groups, ct, has_tfbs, target_chrom=target_chrom)
                best_row = regions_df.iloc[best_idx]
                best_region_str = f"{best_row['Chromosome']}:{best_row['Start']}-{best_row['End']}"

                if target_gene:
                    print(f"  Top target gene: {target_gene}")

                # Look up target gene TSS (or fall back to TF TSS if no target)
                tss_gene = target_gene or tf
                tss_info = parse_tss_from_gtf(args.gtf, tss_gene) if args.gtf else None
                has_tss = tss_info is not None

                n_rows = 2 if has_tss else 1
                fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3 * n_rows))
                if n_rows == 1:
                    axes = [axes]

                row_idx = 0

                # Row 1: Target gene TSS-centered insertion profile
                if has_tss:
                    tss_chrom, tss_pos, tss_strand = tss_info
                    tss_region = f"{tss_chrom}:{max(0, tss_pos - 2000)}-{tss_pos + 2000}"
                    if target_gene:
                        tss_label = f"target gene {target_gene}"
                    else:
                        tss_label = f"TF {tf} (own promoter)"
                    print(f"  TSS for {tss_label}: {tss_chrom}:{tss_pos} ({tss_strand})")

                    scp.pl.plot_group_atac(
                        printer,
                        [grouping[ct_idx[0]]],
                        np.arange(1),
                        tss_region,
                        ax=[axes[row_idx]],
                        smooth=5,
                        color="blue",
                    )
                    axes[row_idx].set_title(
                        f"{ct} — Tn5 insertion at {tss_label} TSS (±2kb)\n{tss_region}",
                        fontsize=9)
                    axes[row_idx].axvline(x=2000, color='black', linestyle='--',
                                          linewidth=0.8, alpha=0.5, label='TSS')
                    row_idx += 1

                # Row 2 (or 1 if no TSS): Aggregate enhancer Tn5 profile
                # Sample up to 500 regions for performance (4500+ regions × I/O = slow)
                max_sample = 500
                if len(regions_df) > max_sample:
                    sample_idx = np.random.choice(len(regions_df), max_sample, replace=False)
                    regions_sample = regions_df.iloc[sample_idx]
                    print(f"  Computing aggregate Tn5 insertion (sampled {max_sample}/{len(regions_df)} regions)...")
                else:
                    regions_sample = regions_df
                    print(f"  Computing aggregate Tn5 insertion across {len(regions_df)} enhancer regions...")
                agg_profiles = []
                for i, (_, reg) in enumerate(regions_sample.iterrows()):
                    try:
                        reg_str = f"{reg['Chromosome']}:{reg['Start']}-{reg['End']}"
                        counts = scp.pl.get_group_atac(
                            printer, [grouping[ct_idx[0]]], reg_str
                        )
                        if counts is not None:
                            profile = np.asarray(counts).flatten()
                            if len(profile) > 0:
                                agg_profiles.append(profile)
                    except Exception:
                        continue

                if agg_profiles:
                    # Pad/trim to common length and average
                    min_len = min(len(p) for p in agg_profiles)
                    trimmed = [p[:min_len] for p in agg_profiles]
                    mean_profile = np.mean(trimmed, axis=0)
                    sem_profile = np.std(trimmed, axis=0) / np.sqrt(len(trimmed))

                    from scipy.ndimage import gaussian_filter1d
                    mean_smooth = gaussian_filter1d(mean_profile, sigma=5)
                    sem_smooth = gaussian_filter1d(sem_profile, sigma=5)

                    x = np.arange(len(mean_smooth))
                    axes[row_idx].plot(x, mean_smooth, color='red', linewidth=1.2)
                    axes[row_idx].fill_between(x, mean_smooth - sem_smooth,
                                                mean_smooth + sem_smooth,
                                                color='red', alpha=0.2)
                    n_total_enh = len(regions_df)
                    n_used = len(agg_profiles)
                    sampled_note = f" (sampled from {n_total_enh})" if n_used < n_total_enh else ""
                    axes[row_idx].set_title(
                        f"{ct} — Mean Tn5 insertion across {n_used} "
                        f"{tf} enhancer regions{sampled_note}", fontsize=9)
                    axes[row_idx].set_xlabel('Position (bp)', fontsize=8)
                    axes[row_idx].set_ylabel('Mean insertion', fontsize=8)
                else:
                    # Fallback: single best region
                    scp.pl.plot_group_atac(
                        printer,
                        [grouping[ct_idx[0]]],
                        np.arange(1),
                        best_region_str,
                        ax=[axes[row_idx]],
                        smooth=5,
                        color="red",
                    )
                    axes[row_idx].set_title(
                        f"{ct} — Tn5 insertion at best {tf} enhancer\n{best_region_str}",
                        fontsize=9)

                plt.tight_layout()
                fig2_path = plots_dir / f"{ct}_{tf}_insertion_example.png"
                plt.savefig(fig2_path, dpi=200)
                plt.close()
                print(f"  Saved plot: {fig2_path}")
            except Exception as e:
                print(f"  WARNING: Insertion plot failed: {e}")
                import traceback; traceback.print_exc()
                plt.close('all')

        # ---- Cicero panels ----
        cicero_conns = None
        if args.cicero_connections:
            print(f"  Loading Cicero connections from {args.cicero_connections}")
            cicero_conns = load_cicero_connections(args.cicero_connections, regions_df)
            if cicero_conns is not None:
                print(f"  Found {len(cicero_conns)} Cicero connections involving enhancer regions")
            else:
                print(f"  No Cicero connections overlap enhancer regions")

        # ---- Plot 3: Composite figure (FIX-82 rework) ----
        # Layout: A. TF-Enhancer MSFP (cis-filtered best cCRE)
        #         B. TF-Target Gene MSFP (promoter ±2kb TSS)
        #         C. TF-Enhancer TFBS (binding score at best cCRE)
        #         D. TF-Target Gene TFBS (binding score at promoter)
        #         E. Aggregate enhancer Tn5 insertion profile
        #         F. CCAN arcs (with colorbar)
        #         G. TF motif logo
        #
        # FUTURE DIRECTIONS: Consider plotting the top 5 same-chromosome
        # cis-enhancer regions alongside the single 'best' currently shown
        # in panels A/C. This would provide a more comprehensive view of
        # the TF's cis-regulatory landscape at multiple enhancer loci,
        # potentially as a multi-row grid or overlaid traces.
        try:
            pwm_df = load_jaspar_pwm(args.pfm_path, tf) if args.pfm_path else None

            import matplotlib.image as mpimg
            import matplotlib.gridspec as gridspec

            # Determine which panels are available
            has_enhancer_msfp = fig1_path and fig1_path.exists()
            has_target_msfp = (target_gene is not None and args.gtf and
                               parse_tss_from_gtf(args.gtf, target_gene) is not None)
            has_tfbs_enhancer = (has_tfbs and hasattr(printer, 'bindingscoreadata') and
                                 save_key_bs in printer.bindingscoreadata)
            has_tfbs_target = has_tfbs_enhancer and has_target_msfp
            has_insertion = fig2_path and fig2_path.exists()
            has_cicero_arcs = cicero_conns is not None and len(cicero_conns) > 0
            has_logo = pwm_df is not None

            # Resolve target gene TSS once for reuse across panels
            tgt_tss = None
            tgt_region = None
            if has_target_msfp:
                tgt_tss = parse_tss_from_gtf(args.gtf, target_gene)
                tgt_chrom, tgt_pos, tgt_strand = tgt_tss
                tgt_region = f"{tgt_chrom}:{max(0, tgt_pos - 2000)}-{tgt_pos + 2000}"

            # Resolve best enhancer region string for reuse
            best_row_enh = regions_df.iloc[best_msfp_idx]
            best_enh_region = f"{best_row_enh['Chromosome']}:{best_row_enh['Start']}-{best_row_enh['End']}"

            # Build panel list: (label_suffix, height_ratio, panel_type)
            panels = []
            if has_enhancer_msfp:
                panels.append(('enhancer_msfp', 2, 'raster'))
            if has_target_msfp:
                panels.append(('target_msfp', 2, 'live'))
            if has_tfbs_enhancer:
                panels.append(('tfbs_enhancer', 1.5, 'live'))
            if has_tfbs_target:
                panels.append(('tfbs_target', 1.5, 'live'))
            if has_insertion:
                panels.append(('insertion', 3, 'raster'))
            if has_cicero_arcs:
                panels.append(('ccan_arcs', 2.5, 'live'))
            if has_logo:
                panels.append(('logo', 1, 'live'))

            if not panels:
                raise ValueError("No panels to assemble")

            height_ratios = [p[1] for p in panels]
            fig_comp = plt.figure(figsize=(12, sum(height_ratios) * 2.5),
                                   constrained_layout=True)
            gs = gridspec.GridSpec(len(panels), 1, figure=fig_comp,
                                   height_ratios=height_ratios)

            for panel_idx, (panel_type_name, _, _) in enumerate(panels):
                ax = fig_comp.add_subplot(gs[panel_idx])
                label = chr(ord('A') + panel_idx)

                if panel_type_name == 'enhancer_msfp':
                    img = mpimg.imread(str(fig1_path))
                    ax.imshow(img, aspect='auto')
                    ax.axis('off')
                    ax.set_title(f'{label}. TF-Enhancer MSFP — {tf} at {best_enh_region}',
                                 fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'target_msfp':
                    try:
                        scp.pl.plot_footprints(
                            printer,
                            save_key=save_key_fp,
                            group_names=[ct],
                            region=tgt_region,
                            ax=ax,
                            stack=False,
                            cmap="Blues",
                            vmin=0.5,
                            vmax=2.0,
                            add_ticks=True,
                        )
                        ymin, ymax = ax.get_ylim()
                        if ymin < ymax:
                            ax.invert_yaxis()
                        region_dict = {
                            'chr': tgt_tss[0],
                            'start': max(0, tgt_tss[1] - 2000),
                            'end': tgt_tss[1] + 2000,
                            'tss': tgt_tss[1],
                        }
                        _format_msfp_axes(ax, region_dict, target_gene, is_enhancer=False)
                        ax.set_title(
                            f'{label}. TF-Target Gene MSFP — {tf} at {target_gene} promoter (±2kb TSS)',
                            fontsize=12, fontweight='bold', loc='left', pad=15)
                    except Exception as e:
                        print(f"  [WARN] Target gene MSFP failed: {e}")
                        ax.text(0.5, 0.5, f'Target gene MSFP failed: {e}',
                                ha='center', va='center', transform=ax.transAxes)
                        ax.set_title(f'{label}. TF-Target Gene MSFP',
                                     fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'tfbs_enhancer':
                    try:
                        scp.pl.plot_binding_score(
                            printer,
                            save_key=save_key_bs,
                            group_names=[ct],
                            region=best_enh_region,
                            ax=ax,
                        )
                        # FIX-87: Reformat x-axis from bin indices to genomic coords
                        enh_region_dict = {
                            'chr': best_row_enh['Chromosome'],
                            'start': int(best_row_enh['Start']),
                            'end': int(best_row_enh['End']),
                        }
                        _format_tfbs_axes(ax, enh_region_dict, tf, is_enhancer=True)
                        ax.set_title(
                            f'{label}. TF-Enhancer TFBS — {tf} binding score at {best_enh_region}',
                            fontsize=12, fontweight='bold', loc='left', pad=15)
                    except Exception as e:
                        print(f"  [WARN] Enhancer TFBS plot failed: {e}")
                        ax.text(0.5, 0.5, f'Enhancer TFBS failed: {e}',
                                ha='center', va='center', transform=ax.transAxes)
                        ax.set_title(f'{label}. TF-Enhancer TFBS',
                                     fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'tfbs_target':
                    try:
                        scp.pl.plot_binding_score(
                            printer,
                            save_key=save_key_bs,
                            group_names=[ct],
                            region=tgt_region,
                            ax=ax,
                        )
                        # FIX-87: Reformat x-axis from bin indices to genomic coords
                        tgt_region_dict = {
                            'chr': tgt_tss[0],
                            'start': max(0, tgt_tss[1] - 2000),
                            'end': tgt_tss[1] + 2000,
                            'tss': tgt_tss[1],
                        }
                        _format_tfbs_axes(ax, tgt_region_dict, target_gene, is_enhancer=False)
                        ax.set_title(
                            f'{label}. TF-Target Gene TFBS — {tf} binding score at {target_gene} promoter',
                            fontsize=12, fontweight='bold', loc='left', pad=15)
                    except Exception as e:
                        print(f"  [WARN] Target gene TFBS plot failed: {e}")
                        ax.text(0.5, 0.5, f'Target gene TFBS failed: {e}',
                                ha='center', va='center', transform=ax.transAxes)
                        ax.set_title(f'{label}. TF-Target Gene TFBS',
                                     fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'insertion':
                    img = mpimg.imread(str(fig2_path))
                    ax.imshow(img, aspect='auto')
                    ax.axis('off')
                    ax.set_title(f'{label}. Aggregate Enhancer Tn5 Insertion',
                                 fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'ccan_arcs':
                    plot_cicero_arcs(cicero_conns, best_enh_region, regions_df, ax)
                    ax.set_title(f'{label}. {ax.get_title()}',
                                 fontsize=12, fontweight='bold', loc='left', pad=15)

                elif panel_type_name == 'logo':
                    render_motif_logo(pwm_df, tf, ax)
                    ax.set_title(f'{label}. TF Motif: {tf}',
                                 fontsize=12, fontweight='bold', loc='left', pad=15)

            # Suptitle
            suptitle = f"{ct} — TF {tf} Enhancer Footprinting"
            if target_chrom:
                cis_count = (regions_df['Chromosome'] == target_chrom).sum()
                suptitle += f" | {cis_count} cis-regulatory regions on {target_chrom}"
                suptitle += f" (of {len(regions_df)} total)"
            else:
                suptitle += f" | {len(regions_df)} enhancer regions"
            if target_gene:
                suptitle += f"\nTarget gene: {target_gene}"
                if n_linked_genes > 1:
                    from collections import Counter
                    gene_freq = Counter()
                    for g in regions_full['linked_genes'].dropna().astype(str):
                        if g and g != 'nan':
                            for gene in g.split(','):
                                gene = gene.strip()
                                if gene and gene.upper() != tf.upper():
                                    gene_freq[gene] += 1
                    if len(gene_freq) > 1:
                        other_targets = [g for g, _ in gene_freq.most_common(6)
                                         if g != target_gene][:4]
                        if other_targets:
                            suptitle += f" | also: {', '.join(other_targets)}"
            fig_comp.suptitle(suptitle, fontsize=11, fontweight='bold', y=1.02)

            comp_path = plots_dir / f"{ct}_{tf}_composite.png"
            fig_comp.savefig(comp_path, dpi=200, bbox_inches='tight', facecolor='white')
            plt.close(fig_comp)
            print(f"  Saved composite: {comp_path}")
        except Exception as e:
            print(f"  WARNING: Composite figure failed: {e}")
            import traceback; traceback.print_exc()
            plt.close('all')
    else:
        print(f"  No footprint data generated for key '{save_key_fp}'")

    # Write summary
    summary = pd.DataFrame([{
        'cell_type': ct,
        'tf': tf,
        'n_regions': len(regions_df),
        'n_linked_genes': n_linked_genes,
        'n_cell_types': len(uniq_groups),
        'total_cells': sum(len(g) for g in grouping),
        'footprint_saved': fp_saved,
    }])
    summary.to_csv('enhancer_fp_summary.csv', index=False)
    print(f"\nEnhancer footprinting complete: {ct} / {tf}")

    if hasattr(printer, 'close'):
        printer.close()


if __name__ == '__main__':
    main()
