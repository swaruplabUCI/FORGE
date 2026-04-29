#!/usr/bin/env python3
"""
export_atac_bigwigs.py — per-cell-type pseudobulk bigWigs from ATAC fragments.

Uses snapatac2's snap.ex.export_coverage on the AnnDataSet (which carries
per-cell insertion data) to produce normalized bigWig coverage tracks, one
file per group. Cell-type labels are joined from the annotated peak matrix
(scATAnno / celltypist / marker — column name passed via --cell-type-col),
so we stay consistent with the rest of the pipeline's grouping.

Output layout:
  bigwigs/{sanitized_cell_type}.bw
  bigwigs/manifest.json   # { cell_type: path, ... }

The ATAC-only equivalent of FORGE's PYCISTOPIC_PREPARE bigwig stage — no
RNA dependency. Container: snapatac_extended.sif (ships snapatac2 ≥ 2.8
and the BigWig writer backend).
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import snapatac2 as snap


def sanitize(s):
    """Filesystem-safe group name; keeps letters, digits, _.-"""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s))


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--anndataset", required=True,
                   help="atac_complete.h5ads from ATAC_FINAL_PIPELINE")
    p.add_argument("--annotated-peaks", required=True,
                   help="peak_matrix_annotated.h5ad providing cell-type labels")
    p.add_argument("--cell-type-col", default="cell_type_prediction",
                   help="obs column on annotated_peaks to group by")
    p.add_argument("--min-cells", type=int, default=100,
                   help="Skip groups with fewer than this many cells")
    p.add_argument("--bin-size", type=int, default=10,
                   help="snap.ex.export_coverage bin_size (default 10 bp)")
    p.add_argument("--normalization", default="RPKM",
                   choices=["RPKM", "CPM", "BPM", "None"],
                   help="Coverage normalization mode (None for raw counts)")
    p.add_argument("--blacklist", default=None,
                   help="Optional BED of blacklisted regions to exclude")
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--outdir", default="bigwigs")
    return p.parse_args()


def main():
    args = parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[export_bw] Loading AnnDataSet: {args.anndataset}")
    adata = snap.read_dataset(args.anndataset, mode="r")
    print(f"[export_bw]   n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    print(f"[export_bw] Loading labels from: {args.annotated_peaks}")
    ann = ad.read_h5ad(args.annotated_peaks)
    if args.cell_type_col not in ann.obs.columns:
        raise KeyError(
            f"--cell-type-col='{args.cell_type_col}' not in annotated peaks obs "
            f"(available: {list(ann.obs.columns)})"
        )

    # Align labels to AnnDataSet's obs order. Barcodes should match exactly
    # (both come from the same fragment pipeline) — any AnnDataSet cell not
    # in the annotated peak matrix is labeled 'Unknown' and dropped.
    labels = ann.obs[args.cell_type_col].astype(str)
    aligned = labels.reindex(list(adata.obs_names)).fillna("Unknown")
    aligned_list = aligned.tolist()

    counts = pd.Series(aligned_list).value_counts()
    print(f"[export_bw] Cell-type distribution (top 10):")
    for ct, n in counts.head(10).items():
        print(f"           {ct}: {n}")
    if len(counts) > 10:
        print(f"           ... and {len(counts) - 10} more")

    keep_groups = counts[(counts >= args.min_cells) & (counts.index != "Unknown")].index.tolist()
    if not keep_groups:
        raise RuntimeError(
            f"No cell-type groups with >= {args.min_cells} cells — nothing to export."
        )
    print(f"[export_bw] Exporting bigWigs for {len(keep_groups)} groups "
          f"(>= {args.min_cells} cells)")

    # Sanitize labels BEFORE handing them to snapatac2: the Rust backend uses
    # the raw groupname to form both the intermediate .zst and final .bw
    # filenames, so any '/' or ' ' in Allen fine labels (e.g. "L2/3 IT ENT
    # Glut") crashes with RuntimeError("invalid filename: .../<label>.zst").
    raw_to_safe = {ct: sanitize(ct) for ct in keep_groups}
    collisions = {}
    for ct, safe in raw_to_safe.items():
        collisions.setdefault(safe, []).append(ct)
    for safe, cts in collisions.items():
        if len(cts) > 1:
            for i, ct in enumerate(cts[1:], start=1):
                raw_to_safe[ct] = f"{safe}__{i}"
    safe_aligned = [raw_to_safe.get(ct, sanitize(ct)) for ct in aligned_list]
    safe_selections = [raw_to_safe[ct] for ct in keep_groups]

    norm = None if args.normalization == "None" else args.normalization
    blacklist = args.blacklist if args.blacklist and os.path.exists(args.blacklist) else None

    written = snap.ex.export_coverage(
        adata,
        groupby=safe_aligned,
        selections=safe_selections,
        bin_size=args.bin_size,
        blacklist=blacklist,
        normalization=norm,
        output_format="bigwig",
        out_dir=str(out_dir),
        prefix="",
        suffix=".bw",
        n_jobs=args.n_jobs,
    )

    # `written` is keyed by sanitized group name. Map back to original labels
    # for the manifest so downstream consumers can still look up Allen names.
    safe_to_raw = {safe: raw for raw, safe in raw_to_safe.items()}
    manifest = {}
    for safe, path in written.items():
        raw = safe_to_raw.get(safe, safe)
        manifest[raw] = Path(path).name
        print(f"[export_bw]   {raw} -> {Path(path).name}")

    manifest_path = out_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump({
            "bigwigs": manifest,
            "cell_type_col": args.cell_type_col,
            "bin_size": args.bin_size,
            "normalization": args.normalization,
            "n_groups": len(manifest),
        }, f, indent=2)
    print(f"[export_bw] Wrote manifest to {manifest_path}")
    print(f"[export_bw] Done — {len(manifest)} bigWigs in {out_dir}/")


if __name__ == "__main__":
    main()
