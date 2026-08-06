#!/usr/bin/env python
"""
run_pycistopic_finalize_lda.py — Phase 3c of the pyCisTopic fan-out pipeline.

Loads the merged CistopicObject and all pre-computed Topic*.pkl models,
selects the best model (or the user-specified one), adds it to the object,
binarizes topics, computes DARs, writes region_sets/, and optionally
computes gene activity.
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
import pyranges as pr
import scanpy as sc

from pycisTopic.lda_models import evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features,
)
from pycisTopic.gene_activity import get_gene_activity


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycistopic_finalize_lda")


def parse_args():
    p = argparse.ArgumentParser(
        description="pyCisTopic Phase 3c: evaluate models, binarize, DARs, region_sets."
    )
    p.add_argument("--merged-pkl",      required=True,
                   help="merged_cistopic.pkl from Phase 3a.")
    p.add_argument("--topic-pkl-dir",   required=True,
                   help="Directory containing Topic*.pkl files from Phase 3b.")
    p.add_argument("--cell-metadata",   required=True)
    p.add_argument("--tss-bed",         required=True)
    p.add_argument("--blacklist",       required=True)
    p.add_argument("--cell-type-col",   default="cell_type_safe")
    p.add_argument("--condition-col",   default="condition")
    p.add_argument("--selected-topics", type=int, default=None)
    p.add_argument("--n-cpu",           type=int, default=8)
    p.add_argument("--do-gene-activity", action="store_true")
    p.add_argument("--species",         default="mmusculus")
    p.add_argument("--outdir",          required=True)
    return p.parse_args()


def main():
    args = parse_args()
    logger = setup_logger()
    os.makedirs(args.outdir, exist_ok=True)

    # ── Load merged object ────────────────────────────────────────────────
    logger.info("Loading merged_cistopic.pkl from %s", args.merged_pkl)
    with open(args.merged_pkl, "rb") as f:
        cistopic_obj = pickle.load(f)
    logger.info("Loaded: %d cells, %d regions",
                len(cistopic_obj.cell_names), len(cistopic_obj.region_names))

    # ── Load all Topic pkl files ──────────────────────────────────────────
    topic_pkls = sorted(glob.glob(os.path.join(args.topic_pkl_dir, "Topic*.pkl")))
    if not topic_pkls:
        raise FileNotFoundError(f"No Topic*.pkl files in {args.topic_pkl_dir}")
    logger.info("Loading %d topic models...", len(topic_pkls))
    models = []
    for tp in topic_pkls:
        with open(tp, "rb") as f:
            m = pickle.load(f)
        logger.info("  Loaded %s (n_topic=%d)", os.path.basename(tp), m.n_topic)
        models.append(m)

    # ── Select best model ─────────────────────────────────────────────────
    logger.info("Evaluating models, selecting n_topics=%s...",
                args.selected_topics if args.selected_topics is not None
                else "auto (pycisTopic chooses)")
    model = evaluate_models(
        models,
        select_model=args.selected_topics,
        return_model=True,
        plot=False,
    )
    cistopic_obj.add_LDA_model(model)
    logger.info("Model with %d topics added to CistopicObject.", model.n_topic)

    # ── Topic binarization ────────────────────────────────────────────────
    logger.info("Binarizing topics...")
    region_bin_topics_top_3k = binarize_topics(
        cistopic_obj, method="ntop", ntop=3000, plot=False
    )
    region_bin_topics_otsu = binarize_topics(
        cistopic_obj, method="otsu", plot=False
    )

    # ── Impute + normalize + DARs ─────────────────────────────────────────
    # Restrict imputation to the union of top-3k regions across all topics.
    # The full atlas (938K regions × 117K cells × int32) is ~411 GB; the
    # binarized union is typically 50-80K regions → ~35 GB, within budget.
    selected_regions = sorted({
        r
        for df in region_bin_topics_top_3k.values()
        for r in df.index.tolist()
    })
    logger.info(
        "Imputing accessibility for %d HVR candidates "
        "(union of top-3k per topic; full atlas = %d regions).",
        len(selected_regions), len(cistopic_obj.region_names),
    )
    imputed_acc_obj = impute_accessibility(
        cistopic_obj, selected_regions=selected_regions, scale_factor=10**6
    )
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
            split_pattern="-",
        )
    else:
        logger.warning("'%s' not in cell_data — skipping DARs.", args.cell_type_col)
        markers_dict = {}

    # ── Region sets ───────────────────────────────────────────────────────
    logger.info("Writing region sets...")
    region_sets_dir = os.path.join(args.outdir, "region_sets")
    for subdir in ("Topics_otsu", "Topics_top_3k", "DARs_cell_type"):
        os.makedirs(os.path.join(region_sets_dir, subdir), exist_ok=True)

    for topic, df in region_bin_topics_otsu.items():
        region_names_to_coordinates(df.index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(region_sets_dir, "Topics_otsu", f"{topic}.bed"),
            sep="\t", header=False, index=False,
        )
    for topic, df in region_bin_topics_top_3k.items():
        region_names_to_coordinates(df.index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(region_sets_dir, "Topics_top_3k", f"{topic}.bed"),
            sep="\t", header=False, index=False,
        )
    for cell_type, dar_df in markers_dict.items():
        if len(dar_df) == 0:
            continue
        region_names_to_coordinates(dar_df.index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(region_sets_dir, "DARs_cell_type", f"{cell_type}.bed"),
            sep="\t", header=False, index=False,
        )

    # ── Optional gene activity ────────────────────────────────────────────
    if args.do_gene_activity:
        logger.info("Computing gene activity...")
        if args.species in ("mmusculus", "mouse"):
            chromsizes_url = (
                "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
            )
        else:
            chromsizes_url = (
                "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
            )
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

    # ── Save CistopicObject ───────────────────────────────────────────────
    logger.info("Saving cistopic_obj.pkl...")
    with open(os.path.join(args.outdir, "cistopic_obj.pkl"), "wb") as f:
        pickle.dump(cistopic_obj, f)

    logger.info("Phase 3c (finalize) complete.")


if __name__ == "__main__":
    main()
