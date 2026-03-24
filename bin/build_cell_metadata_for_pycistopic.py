#!/usr/bin/env python3
import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc

def sample_key_from_obs_sample(s: str) -> str:
    """
    Derive a sample key from obs['sample'] that includes batch info.
    Handles:
      1) June SMK:   '90plus_L1_SMK_june_filtered_S16.h5ad' -> 'S16'
      2) July multiome: '90plus_WTA_8plex_L1_july_filtered_L1_Donor_0.h5ad' -> 'L1_Donor_0_july'
      3) Nov multiome:  '90plus_WTA_10plex_L4_nov_filtered_L4_Donor_0.h5ad' -> 'L4_Donor_0_nov'
      4) Special L10 nov: '...L10_nov_filtered_L1_Donor_0.h5ad' -> 'L10_Donor_0_nov'
    """
    s = os.path.basename(str(s))

    # 1) June SMK: scan entire string for Sxx
    if "_SMK_june_" in s:
        m = re.search(r"S(\d+)", s)
        if m:
            return f"S{m.group(1)}_june"

    # 2) July multiome
    if "_july_" in s:
        m = re.search(r"L\d+_Donor_\d+", s)
        if m:
            return f"{m.group(0)}_july"

    # 3) Nov multiome
    if "_nov_" in s:
        # special case: L10_nov mislabeled donor prefix
        if "L10_nov" in s:
            dm = re.search(r"Donor_(\d+)", s)
            if dm:
                return f"L10_Donor_{dm.group(1)}_nov"

        m = re.search(r"L\d+_Donor_\d+", s)
        if m:
            return f"{m.group(0)}_nov"

    return s

def update_metadata_with_batch_suffix(meta_df: pd.DataFrame) -> pd.DataFrame:
    meta = meta_df.copy()
    meta["sample_id_unique"] = meta["sample_id"]

    mask_demux = meta["sample_type"] == "demux"

    mask_june = (meta["batch"] == "june") & mask_demux & meta["sample_id"].str.match(r"^S\d+$")
    meta.loc[mask_june, "sample_id_unique"] = meta.loc[mask_june, "sample_id"] + "_june"

    mask_july = (meta["batch"] == "july") & mask_demux & meta["sample_id"].str.match(r"L\d+_Donor_\d+$")
    meta.loc[mask_july, "sample_id_unique"] = meta.loc[mask_july, "sample_id"] + "_july"

    mask_nov = (meta["batch"] == "nov") & mask_demux & meta["sample_id"].str.match(r"L\d+_Donor_\d+$")
    meta.loc[mask_nov, "sample_id_unique"] = meta.loc[mask_nov, "sample_id"] + "_nov"

    # June SMK: 'S16' etc. stay as-is
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-h5ad", required=True)
    ap.add_argument("--sample-metadata", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    print(f"[fix] Loading RNA h5ad: {args.rna_h5ad}", file=sys.stderr)
    adata = sc.read_h5ad(args.rna_h5ad)

    if "sample" not in adata.obs.columns:
        raise ValueError("obs['sample'] not found in h5ad")

    # Auto-detect cell type prediction column
    cell_type_col = None
    for candidate in ['scanvi_prediction', 'celltypist_prediction', 'cell_type_prediction']:
        if candidate in adata.obs.columns:
            cell_type_col = candidate
            break
    if cell_type_col is None:
        raise ValueError("No cell type prediction column found in h5ad obs "
                         "(tried scanvi_prediction, celltypist_prediction, cell_type_prediction)")
    if cell_type_col != 'scanvi_prediction':
        # Alias so downstream code can use 'scanvi_prediction' consistently
        adata.obs['scanvi_prediction'] = adata.obs[cell_type_col]
        print(f"[fix] Using '{cell_type_col}' as cell type column (aliased to scanvi_prediction)", file=sys.stderr)

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
