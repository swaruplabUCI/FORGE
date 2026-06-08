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
        cistopic_obj = cistopic_objects[0].merge(
            cistopic_objects[1:],
            is_acc=1, project="merged", copy=True, split_pattern="-",
        )
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
