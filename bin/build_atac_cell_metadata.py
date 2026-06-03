#!/usr/bin/env python3
"""
build_atac_cell_metadata.py

ATAC-only replacement for build_cell_metadata_for_pycistopic.py.

Reads cell type labels and sample IDs directly from the peak_matrix h5ad
obs (no RNA h5ad required). Used by PYCISTOPIC_ATAC_PREPARE to build the
cell metadata TSV that run_pycistopic_prepare.py consumes.

Expected obs columns in the h5ad:
  - obs_names : barcodes (may be prefixed 'sample_id:barcode')
  - obs['sample'] or obs['sample_id'] : sample identifier
  - obs[cell_type_key] : cell type label (e.g. 'taxonomy', 'cell_type')

Output TSV schema (matches RNA version output by build_cell_metadata_for_pycistopic.py):
  barcode \t sample_id \t cell_type
"""

import argparse
import sys

import anndata as ad
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Build pyCisTopic cell metadata from ATAC peak matrix h5ad (no RNA)."
    )
    ap.add_argument("--atac-h5ad",      required=True,
                    help="peak_matrix_annotated.h5ad (or umap-cleaned variant)")
    ap.add_argument("--sample-metadata", required=True,
                    help="Sample manifest CSV — must have 'sample_id' column")
    ap.add_argument("--cell-type-key",   default="cell_type",
                    help="obs column with cell type labels (default: cell_type)")
    ap.add_argument("--condition-key",   default="condition",
                    help="obs column with condition labels (default: condition). "
                         "Pass 'none' to skip.")
    ap.add_argument("--out-tsv",         required=True,
                    help="Output TSV path")
    args = ap.parse_args()

    print(f"[build_atac_meta] Loading ATAC h5ad: {args.atac_h5ad}", file=sys.stderr)
    adata = ad.read_h5ad(args.atac_h5ad)
    print(f"[build_atac_meta] {adata.n_obs} cells, obs columns: {list(adata.obs.columns)}",
          file=sys.stderr)

    # Resolve cell type column
    if args.cell_type_key not in adata.obs.columns:
        raise ValueError(
            f"--cell-type-key '{args.cell_type_key}' not found in obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Resolve sample column — try manifest-matched names first.
    # 'demux_sample' (e.g. 'sample1') matches manifest sample_id directly;
    # 'sample' in this dataset is 'sample1_batch1' (includes batch suffix).
    for _col in ("sample_id", "demux_sample", "sample"):
        if _col in adata.obs.columns:
            sample_col = _col
            break
    else:
        raise ValueError(
            "obs must contain 'sample_id', 'demux_sample', or 'sample' column. "
            f"Found: {list(adata.obs.columns)}"
        )

    # Load manifest to validate sample overlap
    print(f"[build_atac_meta] Loading manifest: {args.sample_metadata}", file=sys.stderr)
    manifest = pd.read_csv(args.sample_metadata)
    if "sample_id" not in manifest.columns:
        raise ValueError(
            f"manifest missing 'sample_id' column. Found: {list(manifest.columns)}"
        )
    manifest_ids = set(manifest["sample_id"].astype(str))

    # Build output dataframe
    barcodes = adata.obs_names.tolist()
    cell_types = adata.obs[args.cell_type_key].astype(str).tolist()
    sample_ids = adata.obs[sample_col].astype(str).tolist()

    # Strip sample prefix from barcodes if present (e.g. 'sample1:AACG...' → 'AACG...')
    # The fragment files contain raw barcodes; pycisTopic needs matching format.
    stripped_barcodes = []
    n_stripped = 0
    for bc in barcodes:
        if ":" in bc:
            stripped = bc.split(":", 1)[1]
            stripped_barcodes.append(stripped)
            n_stripped += 1
        else:
            stripped_barcodes.append(bc)
    if n_stripped > 0:
        print(
            f"[build_atac_meta] Stripped sample prefix from {n_stripped}/{len(barcodes)} "
            f"barcodes (e.g. '{barcodes[0]}' → '{stripped_barcodes[0]}')",
            file=sys.stderr,
        )

    # Validate sample overlap with manifest
    obs_ids = set(sample_ids)
    overlap = obs_ids & manifest_ids
    if not overlap:
        print(
            f"[build_atac_meta] WARNING: No overlap between obs sample IDs and manifest. "
            f"obs: {sorted(obs_ids)[:5]}, manifest: {sorted(manifest_ids)[:5]}",
            file=sys.stderr,
        )
    else:
        print(
            f"[build_atac_meta] sample_id overlap: {len(overlap)}/{len(manifest_ids)} manifest samples",
            file=sys.stderr,
        )

    row = {"barcode": stripped_barcodes, "sample_id": sample_ids, "cell_type": cell_types}

    if args.condition_key and args.condition_key.lower() != "none":
        if args.condition_key in adata.obs.columns:
            row["condition"] = adata.obs[args.condition_key].astype(str).tolist()
            print(
                f"[build_atac_meta] condition col '{args.condition_key}': "
                f"{sorted(set(row['condition']))}",
                file=sys.stderr,
            )
        else:
            print(
                f"[build_atac_meta] WARNING: --condition-key '{args.condition_key}' not found "
                f"in obs. Condition column will be absent from output.",
                file=sys.stderr,
            )

    out_df = pd.DataFrame(row).set_index("barcode")

    print(
        f"[build_atac_meta] Output: {len(out_df)} cells, "
        f"{out_df['sample_id'].nunique()} samples, "
        f"{out_df['cell_type'].nunique()} cell types",
        file=sys.stderr,
    )

    out_df.to_csv(args.out_tsv, sep="\t")
    print(f"[build_atac_meta] Wrote {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
