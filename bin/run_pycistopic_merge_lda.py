#!/usr/bin/env python
"""
run_pycistopic_merge_lda.py — Phase 3 of the pyCisTopic 3-phase pipeline.

Loads all per-group CistopicObject pkl files produced by Phase 2,
merges them into a single atlas object, runs MALLET LDA topic modeling
on the merged sparse count matrix, binarizes topics, computes DARs,
writes region_sets/, and optionally writes gene_activity.h5ad.

LDA is intentionally run on the MERGED object (not per-group) to capture
cross-group chromatin programs with maximum statistical power.
"""

import argparse
import glob
import logging
import os
import pickle
import sys
import tempfile

import numpy as np
import pandas as pd
import scanpy as sc
import pyranges as pr
from pycisTopic.cistopic_class import CistopicObject
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features,
)
from pycisTopic.gene_activity import get_gene_activity
import subprocess


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycistopic_merge_lda")


def parse_args():
    p = argparse.ArgumentParser(
        description="pyCisTopic Phase 3: merge per-group objects, run LDA, DARs, region_sets."
    )
    p.add_argument("--pkl-dir", required=True,
                   help="Directory containing *_cistopic.pkl files from Phase 2.")
    p.add_argument("--cell-metadata", required=True,
                   help="cell_metadata_for_pycistopic.safe.tsv from Phase 1.")
    p.add_argument("--tss-bed", required=True,
                   help="TSS BED (qc/tss.bed) from Phase 1.")
    p.add_argument("--blacklist", required=True,
                   help="Decompressed blacklist.bed from Phase 1.")
    p.add_argument("--cell-type-col", default="cell_type_safe")
    p.add_argument("--condition-col",  default="condition")
    p.add_argument("--mallet-path",    default="Mallet-202108/bin/mallet")
    p.add_argument("--topics",         default="10,20,30,40")
    p.add_argument("--selected-topics", type=int, default=40)
    p.add_argument("--n-cpu",          type=int, default=8)
    p.add_argument("--do-gene-activity", action="store_true")
    p.add_argument("--species",        default="mmusculus")
    p.add_argument("--outdir",         required=True)
    return p.parse_args()


def ensure_mallet(mallet_path, logger):
    if os.path.exists(mallet_path):
        return mallet_path
    logger.info("MALLET not found at %s — downloading...", mallet_path)
    subprocess.check_call([
        "wget", "https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz"
    ])
    subprocess.check_call(["tar", "-xf", "Mallet-202108-bin.tar.gz"])
    new_path = os.path.join("Mallet-202108", "bin", "mallet")
    if not os.path.exists(new_path):
        raise FileNotFoundError(f"MALLET binary not found at {new_path} after extraction.")
    return new_path


def main():
    args = parse_args()
    logger = setup_logger()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load all per-group pkl files ────────────────────────────────────────
    pkl_files = sorted(glob.glob(os.path.join(args.pkl_dir, "*_cistopic.pkl")))
    if not pkl_files:
        raise FileNotFoundError(f"No *_cistopic.pkl files found in {args.pkl_dir}")
    logger.info("Found %d pkl files to merge.", len(pkl_files))

    cistopic_objects = []
    for pf in pkl_files:
        logger.info("  Loading %s", os.path.basename(pf))
        with open(pf, "rb") as f:
            obj = pickle.load(f)
        logger.info("    → %d cells, %d regions", len(obj.cell_names), len(obj.region_names))
        cistopic_objects.append(obj)

    # ── Merge ───────────────────────────────────────────────────────────────
    logger.info("Merging %d per-group CistopicObjects...", len(cistopic_objects))
    if len(cistopic_objects) == 1:
        cistopic_obj = cistopic_objects[0]
    else:
        cistopic_obj = cistopic_objects[0].merge(
            cistopic_objects[1:],
            is_acc=1,
            project="merged",
            copy=True,
            split_pattern="___",
        )
    logger.info("Merged object: %d cells, %d regions.",
                len(cistopic_obj.cell_names), len(cistopic_obj.region_names))

    # ── Merge cell-type / condition annotations ─────────────────────────────
    logger.info("Loading cell metadata from %s", args.cell_metadata)
    cell_data = pd.read_table(args.cell_metadata, index_col=0, sep="\t")
    if "barcode" in cell_data.columns:
        # Prefer the compound key barcode___sample_id to match pycisTopic's split_pattern.
        # This avoids ambiguous reverse-mapping when the same raw barcode appears in
        # multiple samples (pd.Series.get() returns a Series for duplicate keys).
        if "sample_id" in cell_data.columns:
            cell_data.index = (cell_data["barcode"].astype(str) + "___"
                               + cell_data["sample_id"].astype(str))
        else:
            cell_data.index = cell_data["barcode"].astype(str)
    raw_idx = cell_data.index.astype(str)
    if raw_idx.str.contains(":").any():
        cell_data.index = raw_idx.str.split(":", n=1).str[-1]

    cto_barcodes = cistopic_obj.cell_data.index.astype(str)
    meta_barcodes = cell_data.index.astype(str)
    shared = cto_barcodes.intersection(meta_barcodes)
    logger.info("Barcode overlap (direct): %d / %d cistopic barcodes match metadata.",
                len(shared), len(cto_barcodes))

    # Fallback: strip the split_pattern='___' suffix for single-sample objects whose
    # metadata only has raw barcodes (no sample_id column).
    if len(shared) == 0 and cto_barcodes.str.contains("___").any():
        stripped = cto_barcodes.str.split("___").str[0]
        shared_stripped = stripped.intersection(meta_barcodes)
        logger.info("After '___' strip: %d / %d match.", len(shared_stripped), len(cto_barcodes))
        if len(shared_stripped) > 0:
            # Build a dict (not Series) to avoid duplicate-key ambiguity.
            bc_map = dict(zip(stripped.values, cto_barcodes.values))
            cell_data = cell_data.copy()
            cell_data.index = cell_data.index.map(lambda x: bc_map.get(x, x))
            meta_barcodes = cell_data.index.astype(str)
            shared = cto_barcodes.intersection(meta_barcodes)
            logger.info("After reindex: shared=%d", len(shared))

    cols_to_add = [c for c in cell_data.columns if c not in cistopic_obj.cell_data.columns]
    if cols_to_add and len(shared) > 0:
        cistopic_obj.cell_data = cistopic_obj.cell_data.join(cell_data[cols_to_add], how="left")
        logger.info("Added annotation columns: %s", cols_to_add)

    # ── LDA (MALLET) ────────────────────────────────────────────────────────
    logger.info("Running MALLET LDA topic models...")
    mallet_path = ensure_mallet(args.mallet_path, logger)
    os.environ["MALLET_MEMORY"] = "200G"

    n_topics = [int(x) for x in args.topics.split(",") if x.strip()]
    mallet_tmp  = os.path.join(args.outdir, "mallet_tmp")
    mallet_save = os.path.join(args.outdir, "mallet_models")
    os.makedirs(mallet_tmp, exist_ok=True)
    os.makedirs(mallet_save, exist_ok=True)

    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=n_topics,
        n_cpu=args.n_cpu,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=mallet_tmp,
        save_path=mallet_save,
        mallet_path=mallet_path,
    )
    model = evaluate_models(models, select_model=args.selected_topics, return_model=True)
    cistopic_obj.add_LDA_model(model)

    # ── Topic binarization ──────────────────────────────────────────────────
    logger.info("Binarizing topics...")
    region_bin_topics_top_3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000, plot=False)
    region_bin_topics_otsu   = binarize_topics(cistopic_obj, method="otsu", plot=False)

    # ── Impute + DARs ───────────────────────────────────────────────────────
    logger.info("Imputing accessibility...")
    imputed_acc_obj = impute_accessibility(cistopic_obj, scale_factor=10**6)
    norm_acc_obj    = normalize_scores(imputed_acc_obj, scale_factor=10**4)

    logger.info("Finding highly variable regions...")
    variable_regions = find_highly_variable_features(
        norm_acc_obj,
        min_disp=0.05, min_mean=0.0125, max_mean=3,
        max_disp=np.inf, n_bins=20, n_top_features=None, plot=False,
    )

    logger.info("Computing DARs by '%s'...", args.cell_type_col)
    if args.cell_type_col in cistopic_obj.cell_data.columns:
        markers_dict = find_diff_features(
            cistopic_obj, imputed_acc_obj,
            variable=args.cell_type_col,
            var_features=variable_regions,
            contrasts=None,
            adjpval_thr=0.05,
            log2fc_thr=np.log2(1.5),
            n_cpu=args.n_cpu,
            _temp_dir=tempfile.gettempdir(),
            split_pattern="___",
        )
    else:
        logger.warning("cell_type_col '%s' not in cell_data — skipping DARs.", args.cell_type_col)
        markers_dict = {}

    # ── Region sets ─────────────────────────────────────────────────────────
    logger.info("Writing region sets...")
    region_sets_dir  = os.path.join(args.outdir, "region_sets")
    topics_otsu_dir  = os.path.join(region_sets_dir, "Topics_otsu")
    topics_top_dir   = os.path.join(region_sets_dir, "Topics_top_3k")
    dars_dir         = os.path.join(region_sets_dir, "DARs_cell_type")
    for d in [topics_otsu_dir, topics_top_dir, dars_dir]:
        os.makedirs(d, exist_ok=True)

    for topic, df in region_bin_topics_otsu.items():
        region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            os.path.join(topics_otsu_dir, f"{topic}.bed"), sep="\t", header=False, index=False
        )
    for topic, df in region_bin_topics_top_3k.items():
        region_names_to_coordinates(df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            os.path.join(topics_top_dir, f"{topic}.bed"), sep="\t", header=False, index=False
        )
    for cell_type, dar_df in markers_dict.items():
        if len(dar_df) == 0:
            continue
        region_names_to_coordinates(dar_df.index).sort_values(["Chromosome", "Start", "End"]).to_csv(
            os.path.join(dars_dir, f"{cell_type}.bed"), sep="\t", header=False, index=False
        )

    # ── Optional gene activity ──────────────────────────────────────────────
    if args.do_gene_activity:
        logger.info("Computing gene activity...")
        if args.species in ("mmusculus", "mouse"):
            chromsizes_url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
        else:
            chromsizes_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
        chromsizes = pd.read_table(chromsizes_url, header=None, names=["Chromosome", "End"])
        chromsizes.insert(1, "Start", 0)
        chromsizes_pr = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

        pr_annotation = pd.read_table(args.tss_bed).rename(
            {"Name": "Gene", "# Chromosome": "Chromosome"}, axis=1
        )
        pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
        pr_annotation = pr.PyRanges(pr_annotation)

        gene_act, _ = get_gene_activity(
            imputed_acc_obj, pr_annotation, chromsizes_pr,
            use_gene_boundaries=True,
            upstream=[1000, 100000], downstream=[1000, 100000],
            distance_weight=True, decay_rate=1,
            extend_gene_body_upstream=10000, extend_gene_body_downstream=500,
            gene_size_weight=False, gene_size_scale_factor="median",
            remove_promoters=False, average_scores=True,
            scale_factor=1, extend_tss=[10, 10],
            gini_weight=True, return_weights=True,
            project="Gene_activity",
        )
        adata = sc.AnnData(
            X=gene_act.mtx.T,
            obs=pd.DataFrame(index=gene_act.cell_names),
            var=pd.DataFrame(index=gene_act.feature_names),
        )
        adata.write_h5ad(os.path.join(args.outdir, "pycistopic_gene_activity.h5ad"))
        logger.info("Gene activity written.")

    # ── Save CistopicObject ─────────────────────────────────────────────────
    logger.info("Saving cistopic_obj.pkl...")
    with open(os.path.join(args.outdir, "cistopic_obj.pkl"), "wb") as f:
        pickle.dump(cistopic_obj, f)

    logger.info("Phase 3 complete.")


if __name__ == "__main__":
    main()
