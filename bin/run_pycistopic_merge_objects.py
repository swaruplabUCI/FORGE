#!/usr/bin/env python
"""
run_pycistopic_merge_objects.py — Phase 3a of the pyCisTopic fan-out pipeline.

Merges all per-group CistopicObject pkl files produced by Phase 2 into a
single atlas object and saves merged_cistopic.pkl for downstream parallel
LDA fan-out jobs (one SLURM job per topic count).
"""

import argparse
import glob
import logging
import os
import pickle
import sys

import pandas as pd
from pycisTopic.cistopic_class import CistopicObject


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycistopic_merge_objects")


def parse_args():
    p = argparse.ArgumentParser(
        description="pyCisTopic Phase 3a: merge per-group objects, save merged pkl."
    )
    p.add_argument("--pkl-dir",       required=True,
                   help="Directory containing *_cistopic.pkl files from Phase 2.")
    p.add_argument("--cell-metadata", required=True,
                   help="cell_metadata_for_pycistopic.safe.tsv from Phase 1.")
    p.add_argument("--outdir",        required=True)
    return p.parse_args()


def main():
    args = parse_args()
    logger = setup_logger()
    os.makedirs(args.outdir, exist_ok=True)

    pkl_files = sorted(glob.glob(os.path.join(args.pkl_dir, "*_cistopic.pkl")))
    if not pkl_files:
        raise FileNotFoundError(f"No *_cistopic.pkl files in {args.pkl_dir}")
    logger.info("Found %d pkl files to merge.", len(pkl_files))

    cistopic_objects = []
    for pf in pkl_files:
        logger.info("  Loading %s", os.path.basename(pf))
        with open(pf, "rb") as f:
            obj = pickle.load(f)
        logger.info("    → %d cells, %d regions", len(obj.cell_names), len(obj.region_names))
        cistopic_objects.append(obj)

    logger.info("Merging %d CistopicObjects...", len(cistopic_objects))
    if len(cistopic_objects) == 1:
        cistopic_obj = cistopic_objects[0]
    else:
        # pycisTopic's CistopicObject.merge() appends a numeric suffix (_1, _2, …)
        # to every object's `project` to disambiguate cell_names across the merged
        # objects. When all per-group objects come from the same sample (the common
        # single-sample BD/10x case, where each object is one cell type of one sample),
        # this mangles cell_names from "<bc>-<sample>" to "<bc>-<sample>_N". The RNA
        # side keeps "<sample>:<bc>", so SCENIC+'s prepare_GEX_ACC later finds zero
        # shared cells and aborts. Restore the canonical pre-merge names: capture each
        # source object's cell_names in merge order, then re-apply after the merge.
        # This is exact (these ARE the original names) and self-verifying — a barcode
        # belongs to exactly one group, so the concatenated names are globally unique.
        original_cell_names = [n for obj in cistopic_objects for n in obj.cell_names]
        original_sample_ids = [
            s for obj in cistopic_objects
            for s in (obj.cell_data["sample_id"].astype(str).tolist()
                      if "sample_id" in obj.cell_data.columns
                      else [None] * len(obj.cell_names))
        ]
        cistopic_obj = cistopic_objects[0].merge(
            cistopic_objects[1:],
            is_acc=1, project="merged", copy=True, split_pattern="-",
        )
        merged_barcodes = cistopic_obj.cell_data["barcode"].astype(str).tolist() \
            if "barcode" in cistopic_obj.cell_data.columns else None
        aligned = merged_barcodes is not None and all(
            name.split("-", 1)[0] == bc
            for name, bc in zip(original_cell_names, merged_barcodes))
        if len(original_cell_names) != len(cistopic_obj.cell_names):
            logger.warning(
                "Cell count mismatch (%d originals vs %d merged) — leaving merge "
                "suffixes in place.", len(original_cell_names), len(cistopic_obj.cell_names))
        elif not aligned:
            logger.warning(
                "Merge reordered cells (restored barcode token != merged barcode "
                "column) — leaving merge suffixes in place to avoid mis-labeling.")
        else:
            cistopic_obj.cell_names = list(original_cell_names)
            cistopic_obj.cell_data.index = pd.Index(
                original_cell_names, name=cistopic_obj.cell_data.index.name)
            # Also restore the sample_id column — merge rewrites it to match the
            # disambiguated project (sample_N), which would otherwise desync from the
            # restored cell_names and break the barcode+sample_id keys that merge_lda /
            # finalize_lda reconstruct downstream.
            if "sample_id" in cistopic_obj.cell_data.columns and None not in original_sample_ids:
                cistopic_obj.cell_data["sample_id"] = original_sample_ids
            logger.info(
                "Restored %d canonical cell_names after merge (stripped pycisTopic "
                "disambiguation suffix).", len(original_cell_names))
    logger.info("Merged object: %d cells, %d regions.",
                len(cistopic_obj.cell_names), len(cistopic_obj.region_names))

    # Annotate cells from metadata
    logger.info("Loading cell metadata from %s", args.cell_metadata)
    cell_data = pd.read_table(args.cell_metadata, index_col=0, sep="\t")
    if "barcode" in cell_data.columns:
        cell_data.index = cell_data["barcode"].astype(str)
    raw_idx = cell_data.index.astype(str)
    if raw_idx.str.contains(":").any():
        cell_data.index = raw_idx.str.split(":", n=1).str[-1]

    cols_to_add = [c for c in cell_data.columns if c not in cistopic_obj.cell_data.columns]
    shared = cistopic_obj.cell_data.index.astype(str).intersection(cell_data.index.astype(str))
    logger.info("Barcode overlap: %d / %d cistopic barcodes match metadata.",
                len(shared), len(cistopic_obj.cell_data))
    if cols_to_add and len(shared) > 0:
        cistopic_obj.cell_data = cistopic_obj.cell_data.join(cell_data[cols_to_add], how="left")
        logger.info("Added annotation columns: %s", cols_to_add)

    out_pkl = os.path.join(args.outdir, "merged_cistopic.pkl")
    logger.info("Saving merged_cistopic.pkl → %s", out_pkl)
    with open(out_pkl, "wb") as f:
        pickle.dump(cistopic_obj, f)
    logger.info("Phase 3a complete.")


if __name__ == "__main__":
    main()
