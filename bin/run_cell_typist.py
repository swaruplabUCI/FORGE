#!/usr/bin/env python3
"""
CellTypist-based automated cell type annotation for RNA data.
Used as the default annotator when no reference atlas is provided for scANVI.

CellTypist uses pre-trained logistic regression models on curated reference
atlases. Models are downloaded automatically on first use.

Usage:
    python run_celltypist.py --input concat.h5ad --output annotated.h5ad \
        --species human --model Immune_All_Low.pkl
"""

import argparse
import scanpy as sc
from h5ad_compat import sanitize_adata
import logging
import sys

logger = logging.getLogger(__name__)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Run CellTypist cell type annotation on RNA h5ad"
    )
    ap.add_argument("--input", required=True, help="Input h5ad file")
    ap.add_argument("--output", required=True, help="Output annotated h5ad file")
    ap.add_argument("--species", default="human", choices=["human", "mouse"])
    ap.add_argument(
        "--model",
        default=None,
        help="CellTypist model name (e.g., Immune_All_Low.pkl). "
        "If not specified, auto-selects based on species.",
    )
    ap.add_argument(
        "--majority_voting",
        action="store_true",
        default=True,
        help="Use majority voting for over-clustering refinement (default: True)",
    )
    return ap.parse_args()


# Default model selection by species
DEFAULT_MODELS = {
    "human": "Immune_All_Low.pkl",
    "mouse": "Immune_All_Low.pkl",
}


def run_celltypist_annotation(adata, species, model_name=None, majority_voting=True):
    """
    Run CellTypist annotation on an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input data (raw counts or log-normalized).
    species : str
        'human' or 'mouse'.
    model_name : str or None
        CellTypist model to use. If None, selects default for species.
    majority_voting : bool
        Whether to use majority voting refinement.

    Returns
    -------
    adata : AnnData
        Annotated data with 'celltypist_prediction' and
        'cell_type_prediction' columns in .obs.
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        logger.error(
            "celltypist is not installed. Install with: pip install celltypist"
        )
        sys.exit(1)

    # Select model
    if model_name is None:
        model_name = DEFAULT_MODELS.get(species, "Immune_All_Low.pkl")
    logger.info(f"Using CellTypist model: {model_name}")

    # Load model — skip the download entirely when an absolute path is provided
    import os, scipy.sparse as sp
    if os.path.isabs(model_name) or os.path.exists(model_name):
        logger.info(f"Loading model from path: {model_name}")
    else:
        # Fetch ONLY the model we are about to load. download_models() defaults to
        # model=None, and None means "download every available model" — 61 of them.
        # Measured on the tutorial dataset:
        #     model=None (all)          61 files, 32m08s — longer than CellBender
        #     model=Immune_All_Low.pkl   1 file, 2.7 MB, 5.9 s
        # Both load an identical model (98 cell types); only transfer time differs.
        logger.info(f"Downloading CellTypist model: {model_name}")
        models.download_models(force_update=False, model=model_name)
    model = models.Model.load(model=model_name)

    # CellTypist expects log-normalized data
    # Check if data looks like raw counts (large integers)
    adata_ct = adata.copy()
    if "counts" in adata_ct.layers:
        logger.info("Using 'counts' layer, normalizing for CellTypist...")
        adata_ct.X = adata_ct.layers["counts"].copy()
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)
    elif adata_ct.X.max() > 50:
        logger.info("Data looks like raw counts, normalizing...")
        sc.pp.normalize_total(adata_ct, target_sum=1e4)
        sc.pp.log1p(adata_ct)
    else:
        logger.info("Data appears already log-normalized")

    # CellTypist's gene-matching step can produce a COO matrix internally,
    # which then crashes on the `indata[indata > 10] = 10` clip. Converting
    # to dense sidesteps all sparse-format issues; 2-3k cells is fine for RAM.
    if sp.issparse(adata_ct.X):
        logger.info(f"Converting {type(adata_ct.X).__name__} to dense array for CellTypist")
        adata_ct.X = adata_ct.X.toarray()

    # Run annotation
    logger.info("Running CellTypist annotation...")
    predictions = celltypist.annotate(
        adata_ct, model=model, majority_voting=majority_voting
    )

    # Transfer predictions back to original adata
    result = predictions.to_adata()
    adata.obs["celltypist_prediction"] = result.obs["predicted_labels"]

    if majority_voting and "majority_voting" in result.obs.columns:
        adata.obs["celltypist_majority_voting"] = result.obs["majority_voting"]
        # Use majority voting as the primary prediction
        adata.obs["cell_type_prediction"] = result.obs["majority_voting"]
    else:
        adata.obs["cell_type_prediction"] = result.obs["predicted_labels"]

    # Transfer confidence scores if available
    if "conf_score" in result.obs.columns:
        adata.obs["celltypist_confidence"] = result.obs["conf_score"]

    # Summary
    pred_col = "cell_type_prediction"
    n_types = adata.obs[pred_col].nunique()
    logger.info(f"Annotated {adata.n_obs} cells into {n_types} cell types")
    logger.info(f"\nCell type distribution:\n{adata.obs[pred_col].value_counts()}")

    return adata


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logger.info(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    adata = run_celltypist_annotation(
        adata,
        species=args.species,
        model_name=args.model,
        majority_voting=args.majority_voting,
    )

    logger.info(f"Saving annotated data to {args.output}")
    sanitize_adata(adata, args.output)
    adata.write_h5ad(args.output)

    # Write report
    report_path = args.output.replace(".h5ad", "_celltypist_report.txt")
    with open(report_path, "w") as f:
        f.write("CellTypist Annotation Report\n")
        f.write(f"{'=' * 40}\n")
        f.write(f"Species: {args.species}\n")
        f.write(f"Model: {args.model or DEFAULT_MODELS.get(args.species)}\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Cell types found: {adata.obs['cell_type_prediction'].nunique()}\n\n")
        f.write("Distribution:\n")
        f.write(str(adata.obs["cell_type_prediction"].value_counts()))
    logger.info(f"Report saved to {report_path}")


if __name__ == "__main__":
    main()
