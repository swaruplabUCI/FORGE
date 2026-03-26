#!/usr/bin/env python3
# v2: auto-detect celltypist/scanvi prediction columns (cache invalidation)
"""
Multi-Sample Multiome MuData Construction Pipeline (CLI-enabled)
Combines RNA + ATAC samples into unified MuData object
"""

import numpy as np
import pandas as pd
import scanpy as sc
import mudata as md
import anndata as ad
from pathlib import Path
import warnings
import json
import gc
import shutil
from datetime import datetime
import psutil
import os
import argparse
import sys
import re

warnings.filterwarnings('ignore')

def log_memory(log_file):
    process = psutil.Process(os.getpid())
    mem_gb = process.memory_info().rss / 1024**3
    log_message(f"  Memory usage: {mem_gb:.2f} GB", log_file)

def log_message(message: str, log_file, level: str = "INFO"):
    """Log messages to file and console"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_entry = f"[{timestamp}] {level}: {message}"
    print(log_entry)
    with open(log_file, 'a') as f:
        f.write(log_entry + '\n')

def parse_args():
    parser = argparse.ArgumentParser(
        description='Build unified MuData from RNA and ATAC h5ad files'
    )
    
    parser.add_argument('--rna_files', nargs='+', type=str, required=True,
                   help='List of RNA h5ad files (space-separated)')
    parser.add_argument('--atac_files', nargs='+', type=str, required=True,
                   help='List of ATAC h5ad files (space-separated)')
    parser.add_argument('--scanvi_file', type=str, default=None,
                       help='Path to SCANVI predictions h5ad (optional)')
    parser.add_argument('--metadata_file', type=str, required=True,
                       help='Path to sample metadata CSV')
    
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for MuData')
    parser.add_argument('--output_name', type=str, default='integrated.h5mu',  # CHANGED
                       help='Output filename (default: integrated.h5mu)')
    
    parser.add_argument('--batch_size', type=int, default=10,
                       help='Batch size for concatenation (default: 10)')
    parser.add_argument('--min_retention', type=float, default=0.0,
                       help='Minimum RNA/ATAC overlap retention rate (default: 0.0)')
    
    return parser.parse_args()

# ============================================================================
# NEW FUNCTION: Extract sample ID from filename
# ============================================================================
def extract_sample_id(filename):
    """
    Extract standardized sample ID from RNA or ATAC filename
    
    Examples:
        90plus_L1_SMK_june_filtered_S16.h5ad  → S16_june
        S16_june.h5ad                         → S16_june
        L1_Donor_0_july.h5ad                  → L1_Donor_0_july
        90plus_WTA_10plex_L1_nov_filtered_L1_Donor_0.h5ad → L1_Donor_0_nov
    """
    stem = Path(filename).stem
    
    # Pattern 1: RNA files with _filtered_SXX pattern (June batch)
    match = re.search(r'_filtered_(S\d+)$', stem)
    if match:
        sample_id = match.group(1)
        # Extract batch from filename
        if 'june' in stem.lower():
            return f"{sample_id}_june"
        elif 'july' in stem.lower():
            return f"{sample_id}_july"
        elif 'nov' in stem.lower():
            return f"{sample_id}_nov"
    
    # Pattern 2: RNA files with _filtered_LX_Donor_Y pattern (Nov batch)
    match = re.search(r'_filtered_(L\d+_Donor_\d+)$', stem)
    if match:
        sample_id = match.group(1)
        if 'nov' in stem.lower():
            return f"{sample_id}_nov"
        elif 'july' in stem.lower():
            return f"{sample_id}_july"
    
    # Pattern 3: ATAC files (already simplified)
    # S16_june.h5ad, L1_Donor_0_july.h5ad, etc.
    if re.match(r'(S\d+|L\d+_Donor_\d+)_(june|july|nov)', stem):
        return stem
    
    # Fallback: return as-is
    return stem

# Add bin directory to path if utils_barcode_matching exists
bin_dir = Path(__file__).parent
if (bin_dir / "utils_barcode_matching.py").exists():
    sys.path.insert(0, str(bin_dir))
    from utils_barcode_matching import (
        standardize_barcode,
        extract_sample_id_from_path,
        match_rna_atac_barcodes,
        validate_multiome_overlap
    )

# Then used during sample pairing (lines ~200-250):
def process_sample_pair(rna_file, atac_file, sample_id):
    # Load data
    rna_adata = sc.read_h5ad(rna_file)
    atac_adata = sc.read_h5ad(atac_file)
    
    # Use barcode matching utility
    common_barcodes, rna_idx, atac_idx = match_rna_atac_barcodes(
        rna_barcodes=rna_adata.obs_names.tolist(),
        atac_barcodes=atac_adata.obs_names.tolist(),
        sample_id=sample_id,
        method='positional'  # or 'direct' for 10x data
    )
    
    # Validate overlap
    if not validate_multiome_overlap(
        rna_n=rna_adata.n_obs,
        atac_n=atac_adata.n_obs,
        overlap_n=len(common_barcodes),
        min_retention=0.5
    ):
        log_message(f"⚠ Low RNA/ATAC overlap for {sample_id}", LOG_FILE, "WARN")


def main():
    args = parse_args()
    
    # Setup directories
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_DIR = OUTPUT_DIR / "temp_batches"
    TEMP_DIR.mkdir(exist_ok=True)
    
    LOG_FILE = OUTPUT_DIR / f"mudata_construction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    log_message("="*80, LOG_FILE)
    log_message("MUDATA CONSTRUCTION PIPELINE", LOG_FILE)
    log_message("="*80, LOG_FILE)
    log_message(f"RNA files: {len(args.rna_files)} provided", LOG_FILE)
    log_message(f"ATAC files: {len(args.atac_files)} provided", LOG_FILE)
    log_message(f"Metadata: {args.metadata_file}", LOG_FILE)
    log_message(f"Output: {OUTPUT_DIR / args.output_name}", LOG_FILE)
    
    # Load metadata
    log_message("\nLoading metadata...", LOG_FILE)
    metadata_df = pd.read_csv(args.metadata_file)
    log_message(f"Loaded {len(metadata_df)} entries", LOG_FILE)
    
    # ========================================================================
    # NEW: Build sample-matched file pairs
    # ========================================================================
    log_message("\n" + "="*80, LOG_FILE)
    log_message("BUILDING SAMPLE-MATCHED FILE PAIRS", LOG_FILE)
    log_message("="*80, LOG_FILE)
    
    rna_files = [Path(f) for f in args.rna_files]
    atac_files = [Path(f) for f in args.atac_files]
    
    # Filter out aggregated ATAC files
    atac_files = [f for f in atac_files if f.name not in 
                  ['peak_matrix.h5ad', 'atac_complete.h5ad', 'atac_initial_qc.h5ad']]
    
    log_message(f"RNA files: {len(rna_files)}", LOG_FILE)
    log_message(f"ATAC files (after filtering): {len(atac_files)}", LOG_FILE)
    
    # Build lookup dictionaries
    rna_lookup = {}
    for f in rna_files:
        sample_id = extract_sample_id(f.name)
        rna_lookup[sample_id] = f
    
    atac_lookup = {}
    for f in atac_files:
        sample_id = extract_sample_id(f.name)
        atac_lookup[sample_id] = f
    
    log_message(f"\nRNA samples detected: {len(rna_lookup)}", LOG_FILE)
    log_message(f"ATAC samples detected: {len(atac_lookup)}", LOG_FILE)
    
    # Find matching samples
    common_samples = sorted(set(rna_lookup.keys()) & set(atac_lookup.keys()))
    log_message(f"\nMatching samples (RNA ∩ ATAC): {len(common_samples)}", LOG_FILE)
    
    if len(common_samples) == 0:
        log_message("\n❌ ERROR: No matching samples found!", LOG_FILE, "ERROR")
        log_message("\nRNA samples:", LOG_FILE)
        for sid in sorted(rna_lookup.keys())[:10]:
            log_message(f"  - {sid}", LOG_FILE)
        log_message("\nATAC samples:", LOG_FILE)
        for sid in sorted(atac_lookup.keys())[:10]:
            log_message(f"  - {sid}", LOG_FILE)
        sys.exit(1)
    
    log_message("\nFirst 10 matched samples:", LOG_FILE)
    for sid in common_samples[:10]:
        log_message(f"  {sid}: RNA={rna_lookup[sid].name}, ATAC={atac_lookup[sid].name}", LOG_FILE)
    
    # ========================================================================
    # Load SCANVI predictions if provided
    # ========================================================================
    scanvi_predictions = {}
    if args.scanvi_file and Path(args.scanvi_file).exists():
        log_message(f"\nLoading SCANVI predictions from {args.scanvi_file}", LOG_FILE)
        try:
            scanvi_adata = sc.read_h5ad(args.scanvi_file)
            pred_col = None
            for candidate in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_prediction']:
                if candidate in scanvi_adata.obs.columns:
                    pred_col = candidate
                    break
            if pred_col:
                scanvi_predictions = scanvi_adata.obs[pred_col].to_dict()
                log_message(f"✓ Loaded predictions from '{pred_col}' for {len(scanvi_predictions)} cells", LOG_FILE)
        except Exception as e:
            log_message(f"⚠ Could not load SCANVI file: {e}", LOG_FILE, "WARN")
    
    # ========================================================================
    # NEW: Process MATCHED sample pairs
    # ========================================================================
    log_message("\n" + "="*80, LOG_FILE)
    log_message("PROCESSING MATCHED SAMPLE PAIRS", LOG_FILE)
    log_message("="*80, LOG_FILE)
    
    mudata_list = []
    sample_stats = []
    
    for idx, sample_id in enumerate(common_samples, 1):
        rna_file = rna_lookup[sample_id]
        atac_file = atac_lookup[sample_id]
        
        log_message(f"\n[{idx}/{len(common_samples)}] Processing {sample_id}...", LOG_FILE)
        log_message(f"  RNA:  {rna_file.name}", LOG_FILE)
        log_message(f"  ATAC: {atac_file.name}", LOG_FILE)
        
        try:
            # Load data
            rna_adata = sc.read_h5ad(rna_file)
            atac_adata = sc.read_h5ad(atac_file)
            
            log_message(f"  RNA cells: {rna_adata.n_obs:,}", LOG_FILE)
            log_message(f"  ATAC cells: {atac_adata.n_obs:,}", LOG_FILE)
            
            # ============================================================
            # FIX-P0-5 (C2): Smart barcode matching with normalization
            # ============================================================
            rna_sample_bc = str(rna_adata.obs_names[0])
            atac_sample_bc = str(atac_adata.obs_names[0])

            log_message(f"  RNA barcode example: {rna_sample_bc}", LOG_FILE)
            log_message(f"  ATAC barcode example: {atac_sample_bc}", LOG_FILE)

            # FIX-P0-5: Normalize barcodes by stripping -N suffix before matching
            # This handles 10x multiome where RNA has '-1' suffix but ATAC does not
            def _strip_suffix(bc):
                return re.sub(r'-\d+$', '', str(bc))

            rna_raw = list(rna_adata.obs_names)
            atac_raw = list(atac_adata.obs_names)
            rna_normalized = [_strip_suffix(bc) for bc in rna_raw]
            atac_normalized = [_strip_suffix(bc) for bc in atac_raw]

            # Build lookup: normalized → original index
            rna_norm_to_idx = {bc: i for i, bc in enumerate(rna_normalized)}
            atac_norm_to_idx = {bc: i for i, bc in enumerate(atac_normalized)}

            common_normalized = sorted(set(rna_normalized) & set(atac_normalized))

            if len(common_normalized) > 0:
                # Use normalized matching — subset to common cells via original indices
                rna_idx = [rna_norm_to_idx[bc] for bc in common_normalized]
                atac_idx = [atac_norm_to_idx[bc] for bc in common_normalized]
                rna_adata = rna_adata[rna_idx, :].copy()
                atac_adata = atac_adata[atac_idx, :].copy()
                # Assign consistent barcode names (use normalized form)
                rna_adata.obs_names = common_normalized
                atac_adata.obs_names = common_normalized
                common_barcodes = common_normalized

                # FIX-P0-5: Warn if join yield is suspiciously low
                join_yield = len(common_barcodes) / min(len(rna_raw), len(atac_raw))
                if join_yield < 0.10:
                    log_message(
                        f"  ⚠ Very low RNA-ATAC barcode overlap ({join_yield*100:.1f}%) for {sample_id}. "
                        f"Check that RNA and ATAC data come from the same multiome experiment.",
                        LOG_FILE, "WARN"
                    )
                log_message(f"  ✓ {len(common_barcodes):,} cells matched via normalized barcodes", LOG_FILE)
            else:
                # No overlap even after normalization — fall back to positional matching
                log_message("  ⚠ No barcode overlap after normalization - attempting positional matching", LOG_FILE, "WARN")
                log_message("  ⚠ WARNING: Positional matching assumes cells are in the same order!", LOG_FILE, "WARN")

                n_cells = min(rna_adata.n_obs, atac_adata.n_obs)
                synthetic_barcodes = [f"{sample_id}:cell_{i:06d}" for i in range(n_cells)]

                rna_adata = rna_adata[:n_cells, :].copy()
                atac_adata = atac_adata[:n_cells, :].copy()
                rna_adata.obs_names = synthetic_barcodes
                atac_adata.obs_names = synthetic_barcodes
                common_barcodes = synthetic_barcodes

                log_message(f"  ✓ Created {len(common_barcodes)} positionally matched cell pairs", LOG_FILE)

            if len(common_barcodes) == 0:
                log_message(f"  ⚠ Still no overlapping cells, skipping", LOG_FILE, "WARN")
                continue
            
            # Subset to common cells
            rna_sub = rna_adata[common_barcodes, :].copy()
            atac_sub = atac_adata[common_barcodes, :].copy()
            
            # Add SCANVI predictions
            if scanvi_predictions:
                cell_types = [scanvi_predictions.get(bc, 'Unknown') for bc in common_barcodes]
                rna_sub.obs['scanvi_prediction'] = cell_types
                atac_sub.obs['scanvi_prediction'] = cell_types
            
            # Make barcodes globally unique across samples
            new_barcodes = [f"{sample_id}:{bc}" for bc in common_barcodes]
            rna_sub.obs_names = new_barcodes
            atac_sub.obs_names = new_barcodes
            
            # Create MuData
            mdata = md.MuData({'rna': rna_sub, 'atac': atac_sub})
            mdata.obs['sample_id'] = sample_id
            mdata.obs['sample_idx'] = idx
            
            mudata_list.append(mdata)
            
            retention = len(common_barcodes) / min(rna_adata.n_obs, atac_adata.n_obs)
            sample_stats.append({
                'sample_id': sample_id,
                'n_rna': rna_adata.n_obs,
                'n_atac': atac_adata.n_obs,
                'n_multiome': len(common_barcodes),
                'retention_rate': retention
            })
            
            log_message(f"  ✓ {len(common_barcodes):,} cells ({retention*100:.1f}% retention)", LOG_FILE)
            
            if idx % 10 == 0:
                log_memory(LOG_FILE)
                gc.collect()
                
        except Exception as e:
            log_message(f"  ✗ Error: {e}", LOG_FILE, "ERROR")
            import traceback
            log_message(traceback.format_exc(), LOG_FILE, "ERROR")
            continue
    
    log_message(f"\n✓ Successfully processed {len(mudata_list)} sample pairs", LOG_FILE)
    
    if len(mudata_list) == 0:
        log_message("\n❌ ERROR: No valid MuData objects created!", LOG_FILE, "ERROR")
        sys.exit(1)
    
    # ========================================================================
    # Concatenate samples
    # ========================================================================
    log_message("\n" + "="*80, LOG_FILE)
    log_message("CONCATENATING SAMPLES", LOG_FILE)
    log_message("="*80, LOG_FILE)
    log_memory(LOG_FILE)
    
    if len(mudata_list) < 50:
        log_message("Using direct concatenation", LOG_FILE)
        combined_mdata = md.concat(
            mudata_list,
            join='inner',
            merge='same',
            label=None,
            keys=None,
            index_unique=None
        )
    else:
        # Batched concatenation for large datasets
        log_message(f"Using batched concatenation (batch_size={args.batch_size})", LOG_FILE)
        
        n_batches = (len(mudata_list) + args.batch_size - 1) // args.batch_size
        batch_files = []
        
        for batch_idx in range(n_batches):
            start_idx = batch_idx * args.batch_size
            end_idx = min((batch_idx + 1) * args.batch_size, len(mudata_list))
            batch_mdatas = mudata_list[start_idx:end_idx]
            
            batch_combined = md.concat(
                batch_mdatas,
                join='inner',
                merge='same',
                label=None,
                keys=None,
                index_unique=None
            )
            
            batch_file = TEMP_DIR / f"batch_{batch_idx:03d}.h5mu"
            batch_combined.write_h5mu(batch_file)
            batch_files.append(batch_file)
            
            log_message(f"  Batch {batch_idx+1}/{n_batches}: {batch_combined.n_obs:,} cells", LOG_FILE)
            
            del batch_combined, batch_mdatas
            gc.collect()
        
        # Hierarchical merge
        current_files = batch_files
        merge_round = 1
        
        while len(current_files) > 1:
            next_files = []
            
            for i in range(0, len(current_files), 2):
                if i + 1 < len(current_files):
                    mdata1 = md.read_h5mu(current_files[i])
                    mdata2 = md.read_h5mu(current_files[i+1])
                    
                    merged = md.concat([mdata1, mdata2], join='inner', merge='same',
                                      label=None, keys=None, index_unique=None)
                    
                    merged_file = TEMP_DIR / f"merged_round{merge_round}_{i//2:03d}.h5mu"
                    merged.write_h5mu(merged_file)
                    next_files.append(merged_file)
                    
                    current_files[i].unlink()
                    current_files[i+1].unlink()
                    
                    del mdata1, mdata2, merged
                    gc.collect()
                else:
                    next_files.append(current_files[i])
            
            current_files = next_files
            merge_round += 1
            log_message(f"  Merge round {merge_round-1} complete: {len(next_files)} files", LOG_FILE)
        
        combined_mdata = md.read_h5mu(current_files[0])
        current_files[0].unlink()
    
    log_message(f"\n✓ Final MuData: {combined_mdata.n_obs:,} cells", LOG_FILE)
    log_memory(LOG_FILE)

    # ========================================================================
    # CRITICAL FIX: Propagate cell type annotations to global obs
    # ========================================================================
    log_message("\nPropagating cell type annotations to global obs...", LOG_FILE)
    
    # Check for cell type prediction column in RNA modality and propagate to global
    pred_col_found = None
    for modality_name in ['rna', 'atac']:
        if modality_name not in combined_mdata.mod:
            continue
        mod_obs = combined_mdata.mod[modality_name].obs
        for candidate in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_prediction']:
            if candidate in mod_obs.columns:
                pred_col_found = (modality_name, candidate)
                break
        if pred_col_found:
            break

    if pred_col_found:
        mod_name, col_name = pred_col_found
        values = combined_mdata.mod[mod_name].obs[col_name].values
        combined_mdata.obs['cell_type'] = values
        combined_mdata.obs['scanvi_prediction'] = values  # backward compat
        log_message(f" Propagated '{col_name}' from {mod_name} modality to global obs", LOG_FILE)
        log_message(f"  Unique cell types: {combined_mdata.obs['cell_type'].nunique()}", LOG_FILE)
    else:
        log_message(" No cell type prediction column found in modalities to propagate", LOG_FILE, "WARN")
    
    log_message(f"\n Final MuData: {combined_mdata.n_obs:,} cells", LOG_FILE)
    log_memory(LOG_FILE)
    
    # ========================================================================
    # Save outputs with Nextflow-compatible names
    # ========================================================================
    log_message("\n" + "="*80, LOG_FILE)
    log_message("SAVING OUTPUTS", LOG_FILE)
    log_message("="*80, LOG_FILE)
    
    # Save main MuData file
    output_file = OUTPUT_DIR / args.output_name

    # --------------------------------------------------------------------
    # Ensure 'sample_id' exists in modality obs (needed by MultiVI)
    # --------------------------------------------------------------------
    if "sample_id" in combined_mdata.obs.columns:
        # Propagate to RNA modality
        rna_obs = combined_mdata.mod["rna"].obs
        if "sample_id" not in rna_obs.columns:
            rna_obs["sample_id"] = combined_mdata.obs.loc[rna_obs.index, "sample_id"].values
        
        # Propagate to ATAC modality (not strictly required for batch, but useful)
        atac_obs = combined_mdata.mod["atac"].obs
        if "sample_id" not in atac_obs.columns:
            atac_obs["sample_id"] = combined_mdata.obs.loc[atac_obs.index, "sample_id"].values
    else:
        raise RuntimeError("Expected 'sample_id' in combined_mdata.obs when building MuData")
    
    # Save main MuData file
    output_file = OUTPUT_DIR / args.output_name
    combined_mdata.write_h5mu(output_file)
    log_message(f"✓ Saved: {output_file}", LOG_FILE)
    
    # Save statistics as CSV
    stats_df = pd.DataFrame(sample_stats)
    stats_csv = OUTPUT_DIR / "sample_statistics.csv"
    stats_df.to_csv(stats_csv, index=False)
    log_message(f"✓ Saved statistics CSV: {stats_csv}", LOG_FILE)
    
    # NEW: Save statistics as JSON (Nextflow expects this)
    stats_json = {
        'n_samples': len(sample_stats),
        'total_cells': int(combined_mdata.n_obs),
        'total_rna_genes': int(combined_mdata.mod['rna'].n_vars),
        'total_atac_peaks': int(combined_mdata.mod['atac'].n_vars),
        'mean_cells_per_sample': float(stats_df['n_multiome'].mean()),
        'mean_retention_rate': float(stats_df['retention_rate'].mean()),
        'samples': sample_stats
    }
    
    stats_json_file = OUTPUT_DIR / "mudata_stats.json"
    with open(stats_json_file, 'w') as f:
        json.dump(stats_json, f, indent=2)
    log_message(f" Saved statistics JSON: {stats_json_file}", LOG_FILE)
    
    # Save integration summary
    summary = {
        'timestamp': datetime.now().isoformat(),
        'output_file': str(output_file),
        'total_samples': len(sample_stats),
        'total_cells': int(combined_mdata.n_obs),
        'total_rna_genes': int(combined_mdata.mod['rna'].n_vars),
        'total_atac_peaks': int(combined_mdata.mod['atac'].n_vars),
        'mean_cells_per_sample': float(stats_df['n_multiome'].mean()),
        'mean_retention_rate': float(stats_df['retention_rate'].mean()),
    }
    
    summary_file = OUTPUT_DIR / "integration_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    log_message(f"✓ Saved summary: {summary_file}", LOG_FILE)
    
    # Cleanup temp directory
    if TEMP_DIR.exists():
        shutil.rmtree(TEMP_DIR)
    
    log_message("\n" + "="*80, LOG_FILE)
    log_message(" MUDATA CONSTRUCTION COMPLETE!", LOG_FILE)
    log_message("="*80, LOG_FILE)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
