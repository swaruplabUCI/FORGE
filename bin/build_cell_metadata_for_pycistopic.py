#!/usr/bin/env python3
"""
Build cell-level metadata TSV for pycisTopic from a CellTypist-annotated RNA h5ad
and a sample manifest CSV.

Outputs: barcode → sample_id, cell_type
"""
import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-h5ad", required=True)
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--cell-type-key", default=None,
                    help="Explicit cell type column from Nextflow (e.g. celltypist_prediction)")
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    print(f"[build_meta] Loading RNA h5ad: {args.rna_h5ad}", file=sys.stderr)
    adata = sc.read_h5ad(args.rna_h5ad)

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
            print(f"[build_meta] Using explicit cell_type_key='{cell_type_col}'", file=sys.stderr)

    # 2) Fallback: pick first candidate with real labels (skip all-Unknown)
    if cell_type_col is None:
        for candidate in candidates:
            if candidate in adata.obs.columns:
                vals = adata.obs[candidate].dropna().unique()
                if len(vals) > 1 or (len(vals) == 1 and vals[0] != 'Unknown'):
                    cell_type_col = candidate
                    print(f"[build_meta] Auto-detected cell type column: '{cell_type_col}'", file=sys.stderr)
                    break
                else:
                    print(f"[build_meta] Skipping '{candidate}' (all Unknown)", file=sys.stderr)

    if cell_type_col is None:
        raise ValueError("No cell type column with real labels found in h5ad obs "
                         f"(tried explicit='{args.cell_type_key}', fallback={candidates})")

    # Alias to 'scanvi_prediction' for downstream consistency
    if cell_type_col != 'scanvi_prediction':
        adata.obs['scanvi_prediction'] = adata.obs[cell_type_col]
        print(f"[build_meta] Aliased '{cell_type_col}' -> 'scanvi_prediction'", file=sys.stderr)

    obs = adata.obs.copy()
    barcodes = adata.obs_names.astype(str)

    # Defensive: strip h5ad filename suffixes from barcodes if present
    # (upstream concat_batches.py should already produce clean names, but guard against it)
    h5ad_suffixes = ['_filtered_All.h5ad', '_filtered.h5ad', '.h5ad']
    def _strip_h5ad_suffix(bc):
        for sfx in h5ad_suffixes:
            if bc.endswith(sfx):
                return bc[:-len(sfx)]
        return bc
    n_before = sum(1 for b in barcodes if any(b.endswith(s) for s in h5ad_suffixes))
    if n_before > 0:
        barcodes = barcodes.map(_strip_h5ad_suffix)
        print(f"[build_meta] Stripped h5ad filename suffix from {n_before} barcodes "
              f"(e.g. {adata.obs_names[0]} -> {barcodes[0]})", file=sys.stderr)

    # Load manifest and resolve sample_id for each cell
    print(f"[build_meta] Loading sample metadata: {args.sample_metadata}", file=sys.stderr)
    meta = pd.read_csv(args.sample_metadata)

    for col in ("sample_id",):
        if col not in meta.columns:
            raise ValueError(f"sample_metadata missing '{col}' column")

    manifest_ids = set(meta["sample_id"].drop_duplicates().unique())
    print(f"[build_meta] Manifest sample_ids: {sorted(manifest_ids)}", file=sys.stderr)

    # Resolve sample_key: try obs['sample_id'] first, then obs['sample']
    if "sample_id" in obs.columns:
        obs_ids = set(obs["sample_id"].dropna().unique())
        if obs_ids & manifest_ids:
            obs["sample_key"] = obs["sample_id"]
            print(f"[build_meta] Using obs['sample_id'] directly "
                  f"({len(obs_ids & manifest_ids)}/{len(manifest_ids)} match)", file=sys.stderr)
        else:
            obs["sample_key"] = obs["sample"].astype(str)
            print("[build_meta] obs['sample_id'] doesn't match manifest, using obs['sample']",
                  file=sys.stderr)
    else:
        obs["sample_key"] = obs["sample"].astype(str)
        print("[build_meta] No obs['sample_id'], using obs['sample']", file=sys.stderr)

    # If sample_keys don't match manifest, try stripping h5ad filename suffixes
    derived_keys = set(obs["sample_key"].unique())
    if not (derived_keys & manifest_ids):
        print("[build_meta] No sample_key matches manifest — trying suffix stripping", file=sys.stderr)
        def strip_suffixes(s):
            # Same regex as concat_batches.py: strips .h5ad and _filtered/_qc_filtered variants
            return re.sub(r'(_qc_filtered|_filtered(_\w+)?)?\.h5ad$', '', s)
        obs["sample_key"] = obs["sample_key"].map(strip_suffixes)
        stripped_keys = set(obs["sample_key"].unique())
        if stripped_keys & manifest_ids:
            print(f"[build_meta] Suffix stripping matched: {stripped_keys & manifest_ids}", file=sys.stderr)
        else:
            print(f"[build_meta] Still no match. Keys: {stripped_keys}, Manifest: {manifest_ids}",
                  file=sys.stderr)

    # Verify all cells matched
    matched_keys = set(obs["sample_key"].unique()) & manifest_ids
    n_unmatched = obs[~obs["sample_key"].isin(manifest_ids)].shape[0]
    if n_unmatched > 0:
        bad_keys = obs.loc[~obs["sample_key"].isin(manifest_ids), "sample_key"].unique()
        print(f"[build_meta] WARNING: {n_unmatched} cells not matched to manifest.", file=sys.stderr)
        print(f"[build_meta] Unmatched keys: {list(bad_keys[:10])}", file=sys.stderr)
        raise SystemExit(1)

    cell_meta = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": obs["sample_key"].values,
            "cell_type": obs["scanvi_prediction"].astype(str).values,
        }
    ).set_index("barcode")

    print(
        f"[build_meta] Final cell metadata: {cell_meta.shape[0]} cells, "
        f"{cell_meta['sample_id'].nunique()} samples, "
        f"{cell_meta['cell_type'].nunique()} cell types.",
        file=sys.stderr,
    )

    cell_meta.to_csv(args.out_tsv, sep="\t")
    print(f"[build_meta] Wrote {args.out_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()
