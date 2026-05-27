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
  bigwigs/{sanitized_cell_type}__{sanitized_condition}.bw   # when --condition-col given
  bigwigs/manifest.json

Manifest structure:
  {
    "bigwigs":       { cell_type: filename, ... },        # per-CT tracks (always)
    "by_condition":  { cell_type: { condition: filename } }, # when --condition-col given
    "cell_type_col": "...",
    "condition_col": "...",   # null when --condition-col not given
    "bin_size":      10,
    "normalization": "RPKM",
    "n_groups":      N
  }

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


def _deduplicate_safe_map(raw_to_safe):
    """Resolve any sanitize() collisions by appending __1, __2, …"""
    collisions = {}
    for ct, safe in raw_to_safe.items():
        collisions.setdefault(safe, []).append(ct)
    for safe, cts in collisions.items():
        if len(cts) > 1:
            for i, ct in enumerate(cts[1:], start=1):
                raw_to_safe[ct] = f"{safe}__{i}"
    return raw_to_safe


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--anndataset", required=True,
                   help="atac_complete.h5ads from ATAC_FINAL_PIPELINE")
    p.add_argument("--annotated-peaks", required=True,
                   help="peak_matrix_annotated.h5ad providing cell-type labels")
    p.add_argument("--cell-type-col", default="cell_type_prediction",
                   help="obs column on annotated_peaks to group by (cell type)")
    p.add_argument("--condition-col", default=None,
                   help="obs column on annotated_peaks carrying condition labels "
                        "(e.g. 'condition'). When provided, also exports per-CT "
                        "per-condition BigWigs; adds 'by_condition' to manifest.")
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


def _export_groups(adata, out_dir, aligned_list, safe_aligned,
                   raw_to_safe, keep_groups,
                   bin_size, blacklist, norm, n_jobs, prefix=""):
    """Call snap.ex.export_coverage for one groupby layer; return manifest dict."""
    safe_to_raw    = {safe: raw for raw, safe in raw_to_safe.items()}
    safe_selections = [raw_to_safe[ct] for ct in keep_groups]

    written = snap.ex.export_coverage(
        adata,
        groupby=safe_aligned,
        selections=safe_selections,
        bin_size=bin_size,
        blacklist=blacklist,
        normalization=norm,
        output_format="bigwig",
        out_dir=str(out_dir),
        prefix=prefix,
        suffix=".bw",
        n_jobs=n_jobs,
    )

    manifest = {}
    for safe, path in written.items():
        raw = safe_to_raw.get(safe, safe)
        manifest[raw] = Path(path).name
        print(f"[export_bw]   {raw} -> {Path(path).name}")
    return manifest


def main():
    args = parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    norm      = None if args.normalization == "None" else args.normalization
    blacklist = args.blacklist if args.blacklist and os.path.exists(args.blacklist) else None

    print(f"[export_bw] Loading AnnDataSet: {args.anndataset}")
    adata = snap.read_dataset(args.anndataset, mode="r")
    print(f"[export_bw]   n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    print(f"[export_bw] Loading labels from: {args.annotated_peaks}")
    ann = ad.read_h5ad(args.annotated_peaks)

    for col in [args.cell_type_col] + ([args.condition_col] if args.condition_col else []):
        if col not in ann.obs.columns:
            raise KeyError(
                f"Column '{col}' not in annotated peaks obs "
                f"(available: {list(ann.obs.columns)})"
            )

    # ── Per-CT BigWigs (always) ────────────────────────────────────────────────
    ct_labels   = ann.obs[args.cell_type_col].astype(str)
    ct_aligned  = ct_labels.reindex(list(adata.obs_names)).fillna("Unknown")
    ct_list     = ct_aligned.tolist()

    ct_counts = pd.Series(ct_list).value_counts()
    print(f"[export_bw] Cell-type distribution (top 10):")
    for ct, n in ct_counts.head(10).items():
        print(f"           {ct}: {n}")
    if len(ct_counts) > 10:
        print(f"           ... and {len(ct_counts) - 10} more")

    ct_keep = ct_counts[(ct_counts >= args.min_cells) & (ct_counts.index != "Unknown")].index.tolist()
    if not ct_keep:
        raise RuntimeError(
            f"No cell-type groups with >= {args.min_cells} cells — nothing to export.")
    print(f"[export_bw] Exporting {len(ct_keep)} per-CT BigWigs")

    ct_raw_to_safe = _deduplicate_safe_map({ct: sanitize(ct) for ct in ct_keep})
    ct_safe_list   = [ct_raw_to_safe.get(ct, sanitize(ct)) for ct in ct_list]

    ct_manifest = _export_groups(
        adata, out_dir, ct_list, ct_safe_list,
        ct_raw_to_safe, ct_keep,
        args.bin_size, blacklist, norm, args.n_jobs,
    )

    # ── Per-CT × condition BigWigs (when --condition-col given) ───────────────
    by_condition_manifest = {}

    if args.condition_col:
        cond_labels   = ann.obs[args.condition_col].astype(str)
        cond_aligned  = cond_labels.reindex(list(adata.obs_names)).fillna("Unknown")
        cond_list     = cond_aligned.tolist()

        # Composite label: CT___condition (triple underscore to avoid clashes)
        composite_list = [
            f"{ct}___{cond}" if ct != "Unknown" and cond != "Unknown" else "Unknown"
            for ct, cond in zip(ct_list, cond_list)
        ]
        comp_counts = pd.Series(composite_list).value_counts()
        comp_keep   = comp_counts[
            (comp_counts >= args.min_cells) & (comp_counts.index != "Unknown")
        ].index.tolist()

        if not comp_keep:
            print(f"[export_bw][WARN] --condition-col given but no (CT, condition) "
                  f"stratum has >= {args.min_cells} cells; skipping per-condition export.")
        else:
            print(f"[export_bw] Exporting {len(comp_keep)} per-CT×condition BigWigs")
            # Sanitize composite labels — the separator ___ survives sanitize()
            comp_raw_to_safe = _deduplicate_safe_map(
                {comp: sanitize(comp) for comp in comp_keep})
            comp_safe_list = [
                comp_raw_to_safe.get(comp, sanitize(comp)) for comp in composite_list]

            comp_manifest = _export_groups(
                adata, out_dir, composite_list, comp_safe_list,
                comp_raw_to_safe, comp_keep,
                args.bin_size, blacklist, norm, args.n_jobs,
            )

            # Re-nest as { ct: { condition: filename } }
            for raw_comp, filename in comp_manifest.items():
                ct_part, cond_part = raw_comp.split("___", 1)
                by_condition_manifest.setdefault(ct_part, {})[cond_part] = filename

    # ── Manifest ──────────────────────────────────────────────────────────────
    manifest = {
        "bigwigs":      ct_manifest,
        "by_condition": by_condition_manifest,   # {} when --condition-col not given
        "cell_type_col": args.cell_type_col,
        "condition_col": args.condition_col,
        "bin_size":      args.bin_size,
        "normalization": args.normalization,
        "n_groups":      len(ct_manifest),
    }
    manifest_path = out_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[export_bw] Wrote manifest to {manifest_path}")
    print(f"[export_bw] Done — {len(ct_manifest)} per-CT BigWigs"
          + (f", {sum(len(v) for v in by_condition_manifest.values())} per-condition BigWigs"
             if by_condition_manifest else "")
          + f" in {out_dir}/")


if __name__ == "__main__":
    main()
