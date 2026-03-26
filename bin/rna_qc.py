#!/usr/bin/env python
"""RNA QC with demultiplexing and comprehensive metrics matching template workflow"""

import argparse
import pandas as pd
import scanpy as sc
import numpy as np
import zipfile
import os
import re
import scipy.io
from anndata import AnnData
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
from scipy.stats import median_abs_deviation
import json
import shutil
import glob
import gzip
from pandas.api.types import is_bool_dtype

def read_bd_mex(mex_dir):
    """Read BD Rhapsody MEX format manually"""
    # Read the matrix
    matrix = scipy.io.mmread(os.path.join(mex_dir, 'matrix.mtx'))
    
    # Read barcodes (cell indices)
    barcodes = pd.read_csv(os.path.join(mex_dir, 'barcodes.tsv'), header=None, sep='\t')
    
    # Read features
    features = pd.read_csv(os.path.join(mex_dir, 'features.tsv'), header=None, sep='\t')
    
    # Create AnnData object (transpose matrix to be cells x genes)
    adata = AnnData(X=matrix.T.tocsr())
    
    # Add counts layer (copy of X before any normalization)
    adata.layers['counts'] = adata.X.copy()
    
    # Set cell indices as strings for obs index
    adata.obs.index = barcodes[0].astype(str).values
    adata.obs['cell_index'] = barcodes[0].values
    
    # Handle features
    if features.shape[1] >= 2:
        adata.var.index = features[1].values  # Gene symbols
        adata.var['gene_ids'] = features[0].values
    else:
        adata.var.index = features[0].values
    
    if features.shape[1] >= 3:
        adata.var['feature_types'] = features[2].values
    
    # Make var names unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    
    return adata

def is_outlier(adata, metric, nmads):
    """Identify outliers using MAD-based method"""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def save_comprehensive_metrics(metrics, output_path):
    """Save comprehensive QC metrics to CSV file with ALL columns"""
    
    # Define column order (add new filtered columns)
    column_order = [
        'sample',
        'unfiltered_cells', 'unfiltered_genes',
        'pre_mingenes_cells', 'pre_mingenes_genes',
        'post_mingenes_cells', 'frac_filt_mingenes',
        'pre_3cellthresh_cells', 'pre_3cellthresh_genes',
        'post_3cellthresh_cells', 'frac_filt_3cellthresh',
        'Scr_input_cells', 'Scr_input_genes',
        'total_cells_afterScr', 'cells_Scr_out', 'cellfracfiltScr',
        'QC_filt_Input_g',
        'upper_lim_ngbc', 'lower_lim_ngbc', 'upper_lim_umi',
        'post_ngbc_filter_cells', 'post_ngbc_filter_genes', 'ngbc_frac_filt',
        'post_pctcountsmt_filter_cells', 'post_pctcountsmt_filter_genes', 'pctmt_frac_filt',
        'post_totalcounts_filter_cells', 'post_totalcounts_filter_genes', 'tc_frac_filt',
        'final_filtered_total', 'final_filtered_total_genes',
        'basicQC_filtered_cells', 'cellfrac_basicQC',
        'scr_filtered_cells', 'cellfrac_scr',
        'advQC_filtered_cells', 'cellfrac_advQC',
        'median_genes_per_cell', 'median_umi_per_cell', 'median_mt_percent',
        # NEW COLUMNS - Split filtering stats
        'n_valid_donors',
        'filtered_unknown', 'filtered_unassigned', 
        'filtered_doublet', 'filtered_multiplet',
        'filtered_low_cell_count', 'filtered_other',
        'total_filtered_during_split'
    ]
    
    # Create DataFrame with metrics
    df = pd.DataFrame([metrics])
    
    # Ensure all columns exist (fill missing with 0 or empty)
    for col in column_order:
        if col not in df.columns:
            df[col] = 0 if col.startswith(('filtered_', 'n_valid')) else ''
    
    # Reorder columns
    df = df[column_order]
    
    # Save to CSV
    df.to_csv(output_path, index=False)
    
    return df

def append_demux_metrics(adata, lane_metrics, metrics_path):
    """Append per-demux_sample QC rows to the lane-level CSV."""
    
    # Debug output
    print(f"\n=== Appending demux metrics for {lane_metrics['sample']} ===")
    print(f"  Total cells in adata: {adata.n_obs}")
    print(f"  Has demux_sample column: {'demux_sample' in adata.obs.columns}")
    
    if 'demux_sample' not in adata.obs.columns:
        print("  -> No demux_sample column, skipping append")
        return

    demux_vals = adata.obs['demux_sample'].astype(str).unique()
    print(f"  Unique demux values: {demux_vals}")
    
    valid_demux = [s for s in demux_vals if s.lower() not in ('unknown', 'unassigned', 'doublet', 'nan', 'none')]
    print(f"  Valid demux (after filtering): {valid_demux}")
    
    if len(valid_demux) == 0:
        print("  -> No valid demux samples, skipping append")
        return

    # Read existing CSV
    df = pd.read_csv(metrics_path)
    print(f"  Existing CSV has {len(df)} rows")

    # Build per-donor rows
    rows = []
    for ds in valid_demux:
        sub = adata[adata.obs['demux_sample'] == ds].copy()
        if sub.n_obs == 0:
            continue

        m = {
            'sample': f"{lane_metrics['sample']}_{ds}",
            'final_filtered_total': int(sub.n_obs),
            'final_filtered_total_genes': int(sub.n_vars),
            'median_genes_per_cell': float(sub.obs['n_genes_by_counts'].median()),
            'median_umi_per_cell': float(sub.obs['total_counts'].median()),
            'median_mt_percent': float(sub.obs['pct_counts_mt'].median())
        }
        rows.append(m)
        print(f"    Added row for {ds}: {sub.n_obs} cells")

    if rows:
        demux_df = pd.DataFrame(rows)
        df = pd.concat([df, demux_df], ignore_index=True)
        df.to_csv(metrics_path, index=False)
        print(f"  ✓ Wrote {len(rows)} donor rows to CSV (total now: {len(df)} rows)")
    else:
        print("  ⚠ No donor rows generated!")

def add_demux_labels(adata, args):
    """
    Populate adata.obs['demux_sample'] if possible.
    Priority:
      1) args.metadata_file (June): CSV with barcode + label columns
      2) args.souporcell_dir (July): clusters.tsv with singlet assignments
    Fallback: leave missing and downstream code will skip donor-stratified QC.
    """
    import os, re
    import pandas as pd

    # Normalize obs barcodes (strip trailing '-1')
    obs_norm = adata.obs_names.to_series().astype(str).str.replace(r'-1$', '', regex=True)

    # 1) Try a metadata CSV (June)
    meta_path = getattr(args, 'metadata_file', None)
    if meta_path and os.path.exists(meta_path):
        md = pd.read_csv(meta_path)

        # If we have orig.ident and a lane-specific sample name (e.g. "90plus_L1_SMK_june"),
        # restrict metadata to that lane to avoid barcode collisions across lanes.
        if 'orig.ident' in md.columns:
            lane_tag = re.sub(r'_june$', '', getattr(args, 'sample', '') or '')
            if lane_tag:
                md = md[md['orig.ident'] == lane_tag]

        # Heuristic: find a barcode column
        bc_candidates = ['barcode', 'cell', 'rna_barcode', 'CellBarcode', 'cell_barcode', 'Unnamed: 0']
        bc_col = next((c for c in bc_candidates if c in md.columns), None)

        # Heuristic: find a label column
        lbl_candidates = [
            'demux_sample', 'sample', 'donor', 'assignment',
            'DemuxSample', 'Sample', 'Assignment', 'Sample_Name'
        ]
        lbl_col = next((c for c in lbl_candidates if c in md.columns), None)

        if bc_col and lbl_col:
            tmp = md[[bc_col, lbl_col]].dropna()
            tmp[bc_col] = tmp[bc_col].astype(str).str.replace(r'-1$', '', regex=True)
            mapping = dict(zip(tmp[bc_col], tmp[lbl_col].astype(str)))
            adata.obs['demux_sample'] = obs_norm.map(mapping).fillna('unassigned').astype('category')
            return adata

    # 2) Try Souporcell clusters.tsv (July)
    sp_dir = getattr(args, 'souporcell_dir', None)
    if sp_dir:
        clusters_path = os.path.join(sp_dir, 'clusters.tsv')
        if os.path.exists(clusters_path):
            cl = pd.read_csv(clusters_path, sep='\t')
            if {'barcode','assignment'}.issubset(cl.columns):
                if 'status' in cl.columns:
                    cl = cl.loc[cl['status'] == 'singlet', ['barcode','assignment']].copy()
                else:
                    cl = cl[['barcode','assignment']].copy()
                cl['barcode_norm'] = cl['barcode'].astype(str).str.replace(r'-1$', '', regex=True)
                # Optional lane-aware label (e.g., L2_Donor_3)
                lane_id = None
                m = re.search(r'(L\d+)', getattr(args, 'sample', '') or '')
                if m:
                    lane_id = m.group(1)
                def mk_label(x):
                    s = str(x)
                    return f"{lane_id}_Donor_{s}" if lane_id and s.isdigit() else s
                cl['label'] = cl['assignment'].apply(mk_label)
                mapping = dict(zip(cl['barcode_norm'], cl['label']))
                adata.obs['demux_sample'] = obs_norm.map(mapping).fillna('unassigned').astype('category')
                return adata

    # Nothing added — leave as-is
    return adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--min_genes', type=int, default=100)
    parser.add_argument('--min_cells', type=int, default=3)
    parser.add_argument('--mt_threshold', type=float, default=10)
    parser.add_argument('--output', required=True)
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--plots_dir', required=True)
    parser.add_argument('--metadata_file', help='Seurat metadata CSV for June samples')
    parser.add_argument('--souporcell_dir', help='Souporcell output directory for July samples')
    parser.add_argument('--species', default='human', choices=['human', 'mouse'])
    
    args = parser.parse_args()
    args.sample = re.sub(r'^_\d+_', '', args.sample)
    sample = args.sample
    plots_dir = args.plots_dir
    
    # Create plots directory
    os.makedirs(args.plots_dir, exist_ok=True)
    
    # Decide how to read input based on extension
    input_ext = os.path.splitext(args.input)[1].lower()

    if input_ext == ".zip":
        # ------------------------------------------------------------------
        # Original path: BD Rhapsody MEX archive 
        # ------------------------------------------------------------------
        print(f"Input appears to be a ZIP (MEX) file: {args.input}")
        
        # Extract and read MEX files
        with zipfile.ZipFile(args.input, 'r') as zip_ref:
            zip_ref.extractall('temp_mex')
        
        # Find if files are gzipped
        matrix_files = glob.glob('temp_mex/**/matrix.mtx*', recursive=True)
        
        if matrix_files and matrix_files[0].endswith('.gz'):
            # Decompress files
            for gz_file in glob.glob('temp_mex/**/*.gz', recursive=True):
                with gzip.open(gz_file, 'rb') as f_in:
                    with open(gz_file[:-3], 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        
        # Read the data in BD MEX format
        adata = read_bd_mex('temp_mex')

    else:
        # ------------------------------------------------------------------
        # New path: HDF5 input (e.g., CellBender-filtered output)
        # ------------------------------------------------------------------
        print(f"Input is not a ZIP file; assuming HDF5/AnnData: {args.input}")
        
        # Try 10x-style HDF5 first (common for CellBender outputs),
        # fall back to h5ad if that fails.
        try:
            adata = sc.read_10x_h5(args.input)
            print("Loaded input with scanpy.read_10x_h5")
        except Exception as e1:
            print(f"  read_10x_h5 failed ({e1}); trying scanpy.read_h5ad...")
            adata = sc.read_h5ad(args.input)
            print("Loaded input with scanpy.read_h5ad")
        
        # Ensure counts layer and unique names
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
    
    # Store initial counts BEFORE any filtering
    unfiltered_cells = adata.n_obs
    unfiltered_genes = adata.n_vars
    
    # Initialize comprehensive metrics dictionary
    metrics = {
        'sample': args.sample,
        'unfiltered_cells': unfiltered_cells,
        'unfiltered_genes': unfiltered_genes
    }
    
    # Add sample name
    adata.obs["sample"] = args.sample
    
    # Populate demux labels if possible (June metadata or July Souporcell)
    adata = add_demux_labels(adata, args)

    # Calculate comprehensive QC metrics
    # Define gene categories based on species
    if args.species == 'human':
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains(r"^HB[^P]")
    else:  # mouse
        adata.var["mt"] = adata.var_names.str.startswith("mt-")
        adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
        adata.var["hb"] = adata.var_names.str.contains(r"(?i)^HB(?!P)")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    
    # Fix NaN in pct_counts_mt
    adata.obs['pct_counts_mt'] = adata.obs['pct_counts_mt'].fillna(0)
    
    # Track cells and genes at each step
    # Step 1: Pre min_genes
    metrics['pre_mingenes_cells'] = adata.n_obs
    metrics['pre_mingenes_genes'] = adata.n_vars
    
    # Filter cells by min_genes
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    
    metrics['post_mingenes_cells'] = adata.n_obs
    metrics['frac_filt_mingenes'] = (metrics['pre_mingenes_cells'] - metrics['post_mingenes_cells']) / metrics['unfiltered_cells'] * 100
    
    # Step 2: Pre 3-cell threshold
    metrics['pre_3cellthresh_cells'] = adata.n_obs
    metrics['pre_3cellthresh_genes'] = adata.n_vars
    
    # Filter genes by min_cells
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    
    metrics['post_3cellthresh_cells'] = adata.n_obs
    metrics['post_3cellthresh_genes'] = adata.n_vars
    metrics['frac_filt_3cellthresh'] = (metrics['pre_3cellthresh_cells'] - metrics['post_3cellthresh_cells']) / metrics['unfiltered_cells'] * 100
    
    # Step 3: Scrublet (per demultiplexed sample where possible)
    metrics['Scr_input_cells'] = adata.n_obs
    metrics['Scr_input_genes'] = adata.n_vars

    print(f"Running Scrublet (per demux_sample) for {args.sample}...")

    # Initialize Scrublet output columns
    adata.obs['doublet_scores'] = np.nan
    adata.obs['predicted_doublets'] = False

    # Minimum cells to keep in a demux group *after* Scrublet filtering
    MIN_CELLS_AFTER_SCRUBLET = 100

    if 'demux_sample' in adata.obs.columns:
        demux_vals = adata.obs['demux_sample'].astype(str).unique()
        # Only run Scrublet on "real" demuxed samples, not junk labels
        valid_demux = [
            s for s in demux_vals
            if s.lower() not in ('unknown', 'unassigned', 'doublet', 'nan', 'none')
        ]

        if len(valid_demux) == 0:
            print("  No valid demux_sample values found; skipping Scrublet.")
        else:
            for ds in valid_demux:
                mask = adata.obs['demux_sample'].astype(str) == ds
                sub = adata[mask].copy()
                n_cells_before = sub.n_obs
                n_genes = sub.n_vars

                print(
                    f"  Scrublet on demux_sample={ds} "
                    f"with {n_cells_before} cells and {n_genes} genes"
                )

                # If there are no genes or no cells, Scrublet cannot run
                if n_cells_before == 0 or n_genes == 0:
                    print(
                        f"    [Skipping Scrublet] demux_sample={ds}: "
                        f"n_cells={n_cells_before}, n_genes={n_genes}; "
                        "keeping all cells (no doublet filtering)."
                    )
                    continue

                # Try running Scrublet
                try:
                    scrub = scr.Scrublet(sub.X)
                    ds_scores, ds_pred = scrub.scrub_doublets()
                except ValueError as e:
                    # Typical scrublet failure when there are 0 features
                    if "0 feature(s)" in str(e):
                        print(
                            f"    [Scrublet failed] demux_sample={ds}: {e}. "
                            "Keeping all cells (no doublet filtering)."
                        )
                        continue
                    else:
                        print(
                            f"    [Scrublet failed] demux_sample={ds}: {e}. "
                            "Keeping all cells (no doublet filtering)."
                        )
                        continue
                except Exception as e:
                    print(
                        f"    [Scrublet failed] demux_sample={ds}: {e}. "
                        "Keeping all cells (no doublet filtering)."
                    )
                    continue

                # Assign Scrublet scores back to the full adata
                adata.obs.loc[sub.obs_names, 'doublet_scores'] = ds_scores

                # Decide whether to filter based on how many cells would remain
                keep_mask_sub = ~ds_pred
                n_cells_after = int(keep_mask_sub.sum())

                if n_cells_after < MIN_CELLS_AFTER_SCRUBLET:
                    # Do NOT filter this demux group; keep all cells
                    print(
                        f"    [Not filtering] demux_sample={ds}: "
                        f"Scrublet would leave {n_cells_after} cells "
                        f"(< {MIN_CELLS_AFTER_SCRUBLET}); keeping all cells."
                    )
                    # We still ran Scrublet and stored scores, but we do not mark doublets
                    adata.obs.loc[sub.obs_names, 'predicted_doublets'] = False
                else:
                    print(
                        f"    [Filtering] demux_sample={ds}: "
                        f"{n_cells_before} -> {n_cells_after} cells after doublet removal"
                    )
                    # Mark predicted doublets for filtering
                    adata.obs.loc[sub.obs_names, 'predicted_doublets'] = ds_pred

                # Try to save a per-demux Scrublet histogram; never let it crash
                ds_safe = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(ds))
                try:
                    try:
                        fig, axes = scrub.plot_histogram()
                    except AttributeError:
                        # Scrublet did not set threshold_; plot without threshold line
                        fig, axes = plt.subplots(1, 1, figsize=(6, 4))
                        axes.hist(
                            scrub.doublet_scores_obs_,
                            bins=50,
                            alpha=0.7,
                            color="gray"
                        )
                        axes.set_xlabel("Scrublet doublet score")
                        axes.set_ylabel("Number of cells")
                        axes.set_title(f"{sample}_{ds_safe}: Scrublet scores (no threshold)")
                    fig.savefig(
                        os.path.join(
                            plots_dir,
                            f"{sample}_{ds_safe}_scrublet_histogram.png"
                        )
                    )
                    plt.close(fig)
                except Exception as e:
                    print(
                        f"    [Warning] Could not generate Scrublet histogram for "
                        f"demux_sample={ds}: {e}"
                    )
    else:
        print("  No demux_sample column; skipping Scrublet entirely (no lane-level Scrublet).")

    # Filter out predicted doublets (only where we explicitly set True)
    pred = adata.obs['predicted_doublets'].fillna(False).astype(bool)
    adata = adata[~pred].copy()

    metrics['total_cells_afterScr'] = adata.n_obs
    metrics['cells_Scr_out'] = metrics['Scr_input_cells'] - metrics['total_cells_afterScr']
    metrics['cellfracfiltScr'] = metrics['cells_Scr_out'] / metrics['unfiltered_cells'] * 100
    
    # Step 4: Advanced QC
    metrics['QC_filt_Input_g'] = adata.n_vars
    
    # Calculate quantiles
    metrics['upper_lim_ngbc'] = np.quantile(adata.obs.n_genes_by_counts.values, 0.975)
    metrics['lower_lim_ngbc'] = np.quantile(adata.obs.n_genes_by_counts.values, 0.025)
    metrics['upper_lim_umi'] = np.quantile(adata.obs.total_counts.values, 0.95)
    
    # Filter by n_genes_by_counts quantiles
    pre_ngbc = adata.n_obs
    adata = adata[
        (adata.obs.n_genes_by_counts >= metrics['lower_lim_ngbc']) & 
        (adata.obs.n_genes_by_counts <= metrics['upper_lim_ngbc'])
    ].copy()
    
    metrics['post_ngbc_filter_cells'] = adata.n_obs
    metrics['post_ngbc_filter_genes'] = adata.n_vars
    metrics['ngbc_frac_filt'] = (pre_ngbc - adata.n_obs) / metrics['unfiltered_cells'] * 100
    
    # FIX-E4: mt_threshold <= 0 disables MT filtering (previously 0 filtered all cells)
    pre_mt = adata.n_obs
    if args.mt_threshold > 0:
        adata = adata[adata.obs["pct_counts_mt"] < args.mt_threshold].copy()
    else:
        print(f"  MT filtering disabled (mt_threshold={args.mt_threshold})")

    metrics['post_pctcountsmt_filter_cells'] = adata.n_obs
    metrics['post_pctcountsmt_filter_genes'] = adata.n_vars
    metrics['pctmt_frac_filt'] = (pre_mt - adata.n_obs) / metrics['unfiltered_cells'] * 100
    
    # Filter by total_counts upper quantile
    pre_tc = adata.n_obs
    adata = adata[adata.obs.total_counts < metrics['upper_lim_umi']].copy()
    
    metrics['post_totalcounts_filter_cells'] = adata.n_obs
    metrics['post_totalcounts_filter_genes'] = adata.n_vars
    metrics['tc_frac_filt'] = (pre_tc - adata.n_obs) / metrics['unfiltered_cells'] * 100
    
    # Final metrics
    metrics['final_filtered_total'] = adata.n_obs
    metrics['final_filtered_total_genes'] = adata.n_vars
    
    # Calculate summary statistics
    metrics['basicQC_filtered_cells'] = metrics['unfiltered_cells'] - metrics['post_3cellthresh_cells']
    metrics['cellfrac_basicQC'] = metrics['basicQC_filtered_cells'] / metrics['unfiltered_cells'] * 100
    
    metrics['scr_filtered_cells'] = metrics['cells_Scr_out']
    metrics['cellfrac_scr'] = metrics['cellfracfiltScr']
    
    metrics['advQC_filtered_cells'] = metrics['total_cells_afterScr'] - metrics['final_filtered_total']
    metrics['cellfrac_advQC'] = metrics['advQC_filtered_cells'] / metrics['unfiltered_cells'] * 100
    
    # Add median statistics
    if adata.n_obs > 0:
        metrics['median_genes_per_cell'] = adata.obs['n_genes_by_counts'].median()
        metrics['median_umi_per_cell'] = adata.obs['total_counts'].median()
        metrics['median_mt_percent'] = adata.obs['pct_counts_mt'].median()
    
    # Step 5: Identify MAD-based outliers (but don't filter)
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5) |
        is_outlier(adata, "log1p_n_genes_by_counts", 5) |
        is_outlier(adata, "pct_counts_mt", 3)
    )
    
    outlier_count = adata.obs["outlier"].sum()
    metrics['mad_outliers'] = outlier_count
    
    # Store lane-level metrics in adata.uns for downstream use
    adata.uns['qc_stats'] = metrics

    # Save lane-level metrics
    lane_df = save_comprehensive_metrics(metrics, args.metrics)

    # Append per-demux metrics (extra rows) 
    append_demux_metrics(adata, metrics, args.metrics)
    
    # ============================================
    # SPLIT AND SAVE PER DEMULTIPLEXED SAMPLE
    # ============================================    
    fallback_label = 'All'    

    # Define junk sample patterns to skip
    SKIP_PATTERNS = ['multiplet', 'undetermined', '_mm', 'sampletag']
    MIN_CELLS_THRESHOLD = 50    

    # Track filtered cells for reporting
    filtered_stats = {
        'unknown': 0,
        'unassigned': 0,
        'doublet': 0,
        'multiplet': 0,
        'low_cell_count': 0,
        'other': 0
    }    

    if 'demux_sample' in adata.obs.columns:
        print(f"\nSplitting into individual demultiplexed samples...")
        demux_samples = adata.obs['demux_sample'].unique()
        print(f"Found {len(demux_samples)} demultiplexed samples")    

        # Filter out junk and low-cell samples
        valid_samples = []
        for s in demux_samples:
            s_lower = str(s).lower()
            sample_adata = adata[adata.obs['demux_sample'] == s]
            n_cells = sample_adata.n_obs
            
            # Track why samples are being filtered
            filtered = False
            reason = None
            
            # Skip known junk patterns
            if any(pattern in s_lower for pattern in SKIP_PATTERNS):
                filtered = True
                reason = 'junk_pattern'
                if 'multiplet' in s_lower:
                    filtered_stats['multiplet'] += n_cells
                else:
                    filtered_stats['other'] += n_cells
                print(f"  Skipping {s} ({n_cells} cells) - matches exclusion pattern")
            
            # Skip 'unknown' and 'unassigned'
            elif s_lower in ('unknown', 'unassigned'):
                filtered = True
                reason = s_lower
                filtered_stats[s_lower] += n_cells
                print(f"  Skipping {s} ({n_cells} cells) - {s_lower}")
            
            # Skip 'doublet'
            elif s_lower in ('doublet',):
                filtered = True
                reason = 'doublet'
                filtered_stats['doublet'] += n_cells
                print(f"  Skipping {s} ({n_cells} cells) - doublet")
            
            # Check cell count threshold
            elif n_cells < MIN_CELLS_THRESHOLD:
                filtered = True
                reason = 'low_cell_count'
                filtered_stats['low_cell_count'] += n_cells
                print(f"  Skipping {s} ({n_cells} cells < {MIN_CELLS_THRESHOLD} threshold)")
            
            if not filtered:
                valid_samples.append(s)
        
        print(f"\nValid samples after filtering: {len(valid_samples)}")
        
        # Print filtered cell summary
        total_filtered = sum(filtered_stats.values())
        if total_filtered > 0:
            print(f"\n=== Filtered Cell Summary ===")
            print(f"Total filtered: {total_filtered} cells")
            for category, count in filtered_stats.items():
                if count > 0:
                    pct = (count / metrics['final_filtered_total']) * 100
                    print(f"  {category}: {count} cells ({pct:.1f}%)")
            print(f"==========================")    

        # Add filtered stats to metrics for CSV output
        metrics['filtered_unknown'] = filtered_stats['unknown']
        metrics['filtered_unassigned'] = filtered_stats['unassigned']
        metrics['filtered_doublet'] = filtered_stats['doublet']
        metrics['filtered_multiplet'] = filtered_stats['multiplet']
        metrics['filtered_low_cell_count'] = filtered_stats['low_cell_count']
        metrics['filtered_other'] = filtered_stats['other']
        metrics['total_filtered_during_split'] = total_filtered
        metrics['n_valid_donors'] = len(valid_samples)
        
        # ============================================
        # UPDATE LANE-LEVEL ROW WITH FILTERED STATS
        # ============================================
        # Read existing CSV (which has lane + donor rows)
        df = pd.read_csv(args.metrics)
        
        # Update only the first row (lane-level) with new filtered stats
        for key in ['filtered_unknown', 'filtered_unassigned', 'filtered_doublet', 
                    'filtered_multiplet', 'filtered_low_cell_count', 'filtered_other',
                    'total_filtered_during_split', 'n_valid_donors']:
            if key in metrics:
                df.loc[0, key] = metrics[key]
        
        # Save updated CSV (preserves donor rows!)
        df.to_csv(args.metrics, index=False)
        print(f"\n✓ Updated lane-level metrics with filtered stats")

        if len(valid_samples) == 0:
            # Fallback: save a single aggregate file with a suffix to satisfy Nextflow pattern
            sample_output = args.output.replace('.h5ad', f'_{fallback_label}.h5ad')
            adata.layers["counts"] = adata.X.copy()
            adata.write(sample_output)
            print(f"  No valid demux groups; saved aggregate dataset to {sample_output}")
        else:
            for sample_name in valid_samples:
                unique_sample_id = f"{args.sample}_{sample_name}"
                sample_adata = adata[adata.obs['demux_sample'] == sample_name].copy()
                sample_adata.layers["counts"] = sample_adata.X.copy()
                sample_output = args.output.replace('.h5ad', f'_{sample_name}.h5ad')
                sample_adata.write(sample_output)
                print(f"  Saved {unique_sample_id}: {sample_adata.n_obs} cells")
    else:
        # If no demultiplexing labels, save a single file with a suffix to match Nextflow's expected pattern
        adata.layers["counts"] = adata.X.copy()
        sample_output = args.output.replace('.h5ad', f'_{fallback_label}.h5ad')
        adata.write(sample_output)
        print(f"No demux_sample column found, saved full dataset to {sample_output}")

    print(f"\nQC Complete for {args.sample}:")
    print(f"  Initial: {metrics['unfiltered_cells']} cells")
    print(f"  Final: {metrics['final_filtered_total']} cells ({metrics['final_filtered_total']/metrics['unfiltered_cells']*100:.1f}%)")
    print(f"  Metrics saved to: {args.metrics}")

    # Clean up
    if os.path.exists('temp_mex'):
        shutil.rmtree('temp_mex')

if __name__ == "__main__":
    main()
