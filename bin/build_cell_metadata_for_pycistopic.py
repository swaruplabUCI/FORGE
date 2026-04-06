#!/usr/bin/env python3
import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc

def sample_key_from_obs_sample(s: str) -> str:
    """Derive a sample key from obs['sample'] by stripping file extensions."""
    s = os.path.basename(str(s))
    # Strip common h5ad suffixes
    for suffix in ['_filtered_All.h5ad', '_filtered.h5ad', '.h5ad']:
        if s.endswith(suffix):
            return s[:-len(suffix)]
    return s

def update_metadata_with_batch_suffix(meta_df: pd.DataFrame) -> pd.DataFrame:
    """Pass through sample_id as sample_id_unique (no batch suffix needed)."""
    meta = meta_df.copy()
    meta["sample_id_unique"] = meta["sample_id"]
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-h5ad", required=True)
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--cell-type-key", default=None,
                    help="Explicit cell type column from Nextflow (e.g. celltypist_prediction)")
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    print(f"[fix] Loading RNA h5ad: {args.rna_h5ad}", file=sys.stderr)
    # Fix anndata version compat: newer anndata (>=0.10) writes uns/log1p/base
    # with null encoding that older anndata (in scenicplus container) can't read.
    # Copy the file first to avoid modifying the original (may be a Nextflow symlink).
    import shutil, h5py
    rna_path = args.rna_h5ad
    try:
        with h5py.File(rna_path, 'r') as f:
            if 'uns/log1p' in f and 'base' in f['uns/log1p']:
                local_copy = "rna_input_compat.h5ad"
                shutil.copy2(rna_path, local_copy)
                with h5py.File(local_copy, 'a') as fc:
                    del fc['uns/log1p/base']
                rna_path = local_copy
                print("[fix] Stripped uns/log1p/base from local copy (anndata compat)", file=sys.stderr)
    except Exception as e:
        print(f"[fix] h5py compat check: {e}", file=sys.stderr)
    adata = sc.read_h5ad(rna_path)

    if "sample" not in adata.obs.columns:
        raise ValueError("obs['sample'] not found in h5ad")

    # Resolve cell type column: prefer explicit key from Nextflow, then smart fallback
    cell_type_col = None
    candidates = ['celltypist_prediction', 'scanvi_prediction', 'cell_type_prediction']

    # 1) If Nextflow passed an explicit key and it exists with real labels, use it
    if args.cell_type_key and args.cell_type_key in adata.obs.columns:
        vals = adata.obs[args.cell_type_key].dropna().unique()
        if len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown'):
            cell_type_col = args.cell_type_key
            print(f"[fix] Using explicit cell_type_key='{cell_type_col}'", file=sys.stderr)

    # 2) Fallback: pick first candidate with real labels (skip all-Unknown)
    if cell_type_col is None:
        for candidate in candidates:
            if candidate in adata.obs.columns:
                vals = adata.obs[candidate].dropna().unique()
                if len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown'):
                    cell_type_col = candidate
                    print(f"[fix] Auto-detected cell type column: '{cell_type_col}'", file=sys.stderr)
                    break
                else:
                    print(f"[fix] Skipping '{candidate}' (all Unknown)", file=sys.stderr)

    if cell_type_col is None:
        raise ValueError("No cell type column with real labels found in h5ad obs "
                         f"(tried explicit='{args.cell_type_key}', fallback={candidates})")

    # Alias to 'scanvi_prediction' for downstream consistency
    if cell_type_col != 'scanvi_prediction':
        adata.obs['scanvi_prediction'] = adata.obs[cell_type_col]
        print(f"[fix] Aliased '{cell_type_col}' -> 'scanvi_prediction'", file=sys.stderr)

    obs = adata.obs.copy()
    barcodes = adata.obs_names.astype(str)

    print(f"[fix] Loading sample metadata: {args.sample_metadata}", file=sys.stderr)
    meta = pd.read_csv(args.sample_metadata)

    for col in ("sample_id", "batch", "sample_type"):
        if col not in meta.columns:
            raise ValueError(f"sample_metadata missing '{col}' column")

    meta = update_metadata_with_batch_suffix(meta)

    demux_meta = meta[meta["sample_type"] == "demux"].copy()
    manifest_ids = set(demux_meta["sample_id_unique"].unique())

    # FIX-70: Use obs['sample_id'] directly if it already matches the manifest
    # (avoids broken regex derivation from obs['sample'] library names)
    if "sample_id" in obs.columns:
        obs_ids = set(obs["sample_id"].dropna().unique())
        overlap = obs_ids & manifest_ids
        if len(overlap) > 0:
            print(f"[fix] obs['sample_id'] already matches manifest "
                  f"({len(overlap)}/{len(manifest_ids)} IDs overlap) — using directly",
                  file=sys.stderr)
            obs["sample_key"] = obs["sample_id"]
        else:
            print("[fix] obs['sample_id'] does not match manifest, deriving from obs['sample']",
                  file=sys.stderr)
            obs["sample_key"] = obs["sample"].map(sample_key_from_obs_sample)
    else:
        print("[fix] Deriving sample_key from obs['sample']", file=sys.stderr)
        obs["sample_key"] = obs["sample"].map(sample_key_from_obs_sample)

    # If derived sample_keys don't match manifest, try stripping h5ad filename suffixes
    derived_keys = set(obs["sample_key"].unique())
    if not (derived_keys & manifest_ids):
        print("[fix] No sample_key matches manifest — trying filename suffix stripping", file=sys.stderr)
        suffix_patterns = ['_filtered_All.h5ad', '_filtered.h5ad', '.h5ad']
        def strip_suffixes(s):
            for suffix in suffix_patterns:
                if s.endswith(suffix):
                    return s[:-len(suffix)]
            return s
        obs["sample_key"] = obs["sample_key"].map(strip_suffixes)
        stripped_keys = set(obs["sample_key"].unique())
        if stripped_keys & manifest_ids:
            print(f"[fix] Suffix stripping matched: {stripped_keys & manifest_ids}", file=sys.stderr)
        else:
            print(f"[fix] Still no match after stripping. Keys: {stripped_keys}, Manifest: {manifest_ids}",
                  file=sys.stderr)

    print("[fix] Example sample_key values:", file=sys.stderr)
    for sk in obs["sample_key"].unique().tolist()[:10]:
        print(f"  {sk}", file=sys.stderr)

    print("[fix] Unique demux sample_id_unique:", file=sys.stderr)
    print(
        demux_meta[["sample_id", "batch", "sample_id_unique"]]
        .drop_duplicates()
        .head(20)
        .to_string(index=False),
        file=sys.stderr,
    )

    # Merge on sample_key -> sample_id_unique
    # Drop sample_id from obs to avoid _x/_y suffix conflict during merge
    obs_for_merge = obs.drop(columns=["sample_id"], errors="ignore")
    merged = obs_for_merge.merge(
        demux_meta[["sample_id", "sample_id_unique"]],
        left_on="sample_key",
        right_on="sample_id_unique",
        how="left",
        validate="many_to_one",
    )

    n_missing = merged["sample_id"].isna().sum()
    if n_missing > 0:
        print(f"[fix] WARNING: {n_missing} cells not matched.", file=sys.stderr)
        print("[fix] Example unmatched sample_keys:", file=sys.stderr)
        bad = merged.loc[merged["sample_id"].isna(), "sample_key"].unique()
        for s in bad[:10]:
            print(f"  {s}", file=sys.stderr)
        raise SystemExit(1)

    # Use sample_id_unique (e.g. L1_Donor_0_july, L1_Donor_0_nov) to match mudata_stats.json
    cell_meta = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": merged["sample_id_unique"].astype(str).values,
            "cell_type": merged["scanvi_prediction"].astype(str).values,
        }
    ).set_index("barcode")

    print(
        f"[fix] Final cell metadata: {cell_meta.shape[0]} cells, "
        f"{cell_meta['sample_id'].nunique()} samples, "
        f"{cell_meta['cell_type'].nunique()} cell types.",
        file=sys.stderr,
    )

    cell_meta.to_csv(args.out_tsv, sep="\t")
    print(f"[fix] Wrote {args.out_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()
