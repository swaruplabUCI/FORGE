#!/usr/bin/env python

import argparse
import os
import sys
import logging
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

# scprinter.datasets eagerly downloads CisBP motif files at import time via pooch.
# On HPC compute nodes there's no internet. The pooch registry caches to
# $XDG_CACHE_HOME/scprinter (/tmp/cache/scprinter). We symlink the pre-populated
# cache there BEFORE importing scprinter so pooch finds every file locally.
_xdg = os.environ.get("XDG_CACHE_HOME", "/tmp/cache")
_pooch_target = os.path.join(_xdg, "scprinter")
_precache = os.environ.get("SCPRINTER_CACHE_DIR", "/dfs7/swaruplab/lesolano/ref/scprinter")
if not os.path.exists(_pooch_target) and os.path.isdir(_precache):
    os.makedirs(_xdg, exist_ok=True)
    os.symlink(_precache, _pooch_target)

import scprinter as scp

def normalize_barcodes(atac, rna, logger):
    """
    Normalize ATAC and/or RNA barcodes so that they share a common key.

    Updated logic (robust for multi-sample runs):
      - If ATAC obs_names look like 'sample:barcode', we PRESERVE that full ID
        to keep them unique.
      - We then create matching 'sample:barcode' IDs in RNA using its sample column.
      - Alignment is done on these composite, unique IDs.
    """
    atac_names = list(atac.obs_names)

    # Quick pattern check on ATAC
    has_colon = sum(':' in x for x in atac_names)
    if has_colon > 0.5 * len(atac_names):  # >50% have colon
        logger.info(
            "Detected 'sample:barcode' pattern in ATAC obs_names; "
            "preserving full 'sample:barcode' IDs for alignment."
        )
    else:
        logger.info(
            "No 'sample:barcode' pattern detected in ATAC obs_names; "
            "using raw obs_names for alignment."
        )

    # We do NOT modify atac.obs_names here anymore
    return atac, rna

def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("scprinter_dorc")

def init_scprinter(cache_dir: str):
    """
    Initialize scprinter as in footprinting script:
    - set scprinter.datasets.datasets().path to cache_dir
    """
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp

def parse_args():
    p = argparse.ArgumentParser(description="Run scPrinter DORC analysis")
    p.add_argument("--atac", required=True, help="ATAC peak matrix h5ad")
    p.add_argument("--rna", required=True, help="RNA h5ad")
    p.add_argument("--genome", required=True, help="Genome key for scPrinter (e.g. 'hg38')")
    p.add_argument("--covariate", default="none", help="Covariate key for regression (or 'none')")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--cache-dir", required=True, help="Path to scPrinter cache directory")
    p.add_argument(
        "--window-pad-size",
        type=int,
        default=250000,
        help="Window padding size in bp for peak-gene correlation (default: 250000)",
    )
    p.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of parallel jobs for scPrinter.fast_gene_peak_corr (default: 8)",
    )
    p.add_argument(
        "--n-bg",
        type=int,
        default=500,
        help="Number of background peaks for Z-score estimation (default: 500)",
    )
    p.add_argument(
        "--pval-cutoff",
        type=float,
        default=0.05,
        help="pvalZ cutoff for significant peak-gene pairs (default: 0.05)",
    )
    return p.parse_args()

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def main():
    args = parse_args()
    logger = setup_logger()

    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # Configure scPrinter cache
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir

    # Initialize scprinter & genome
    global scp
    scp = init_scprinter(args.cache_dir)
    try:
        genome_obj = getattr(scp.genome, args.genome)
    except AttributeError:
        raise ValueError(f"Genome '{args.genome}' not found in scprinter.genome")

    logger.info(f"Loading ATAC AnnData from {args.atac}")
    adata_atac = ad.read_h5ad(args.atac)

    logger.info(f"Loading RNA AnnData from {args.rna}")
    adata_rna = ad.read_h5ad(args.rna)

    # Normalize barcodes if needed (e.g., ATAC 'sample:barcode' vs RNA 'barcode')
    adata_atac, adata_rna = normalize_barcodes(adata_atac, adata_rna, logger)

    # 1. Align barcodes / obs indices using sample:barcode IDs
    logger.info("Aligning barcodes between ATAC and RNA...")

    # ATAC: obs_names already 'sample:barcode' in this pipeline
    atac_ids = pd.Index(adata_atac.obs_names, name="sample:barcode")

    # RNA: either obs_names are already 'sample:barcode', or we build them
    rna_obs_names = pd.Index(adata_rna.obs_names)
    if (rna_obs_names.str.contains(":").sum() > 0.5 * len(rna_obs_names)):
        logger.info(
            "Detected 'sample:barcode' pattern in RNA obs_names; "
            "using RNA obs_names directly for alignment."
        )
        rna_ids = rna_obs_names
    else:
        sample_col_candidates = ["sample", "sample_id", "sample_key"]
        sample_col = next((c for c in sample_col_candidates if c in adata_rna.obs.columns), None)
        if sample_col is None:
            raise ValueError(
                f"RNA AnnData is missing a sample column; tried {sample_col_candidates} "
                f"and found only {list(adata_rna.obs.columns)}"
            )

        rna_ids = pd.Index(
            adata_rna.obs[sample_col].astype(str) + ":" + adata_rna.obs_names.astype(str),
            name="sample:barcode",
        )
    # Drop duplicates conservatively (keep first occurrence)
    atac_ids_unique = atac_ids[~atac_ids.duplicated()]
    rna_ids_unique = rna_ids[~rna_ids.duplicated()]

    shared_ids = atac_ids_unique.intersection(rna_ids_unique)
    logger.info(f"Found {len(shared_ids)} shared sample:barcode IDs.")

    if len(shared_ids) == 0:
        raise ValueError("No shared sample:barcode IDs between ATAC and RNA AnnData.")

    # Subset ATAC and RNA to shared IDs
    # For ATAC, obs_names already equal atac_ids
    adata_atac = adata_atac[shared_ids].copy()

    # For RNA, build a mapping from 'sample:barcode' -> position, then subset
    rna_id_to_pos = pd.Series(index=rna_ids_unique, data=np.arange(len(rna_ids_unique)))
    rna_pos = rna_id_to_pos.loc[shared_ids].values
    adata_rna = adata_rna[rna_pos].copy()

    # Sanity check: obs_names are now aligned in length and order
    logger.info(
        f"ATAC cells after alignment: {adata_atac.n_obs}, "
        f"RNA cells after alignment: {adata_rna.n_obs}"
    )

    # If covariate given, try to propagate (optional)
    if args.covariate != "none":
        if args.covariate in adata_atac.obs.columns:
            adata_rna.obs[args.covariate] = adata_atac.obs[args.covariate]
        elif args.covariate in adata_rna.obs.columns:
            adata_atac.obs[args.covariate] = adata_rna.obs[args.covariate]
        else:
            logger.warning(f"Covariate '{args.covariate}' not found in obs of ATAC or RNA; skipping.")

    # 2. Normalize RNA
    logger.info("Normalizing RNA counts (total-counts + log1p)...")
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)

    # 3. Compute peak-gene correlations
    # Cap n_jobs to numba's thread limit to avoid ValueError from pynndescent
    import numba
    max_threads = numba.config.NUMBA_NUM_THREADS
    n_jobs = min(args.n_jobs, max_threads)
    logger.info(f"Computing peak-gene correlations (n_jobs={n_jobs}, numba_max={max_threads})...")
    start_cols = set(adata_atac.obs.columns)
    dorc_all = scp.dorc.fast_gene_peak_corr(
        adata_atac,
        adata_rna,
        genome=getattr(scp.genome, args.genome),
        tss_df=scp.datasets.FigR_hg38TSSRanges if args.genome == "hg38" else None,
        gene_list=None,
        window_pad_size=args.window_pad_size,
        n_jobs=n_jobs,
        n_bg=args.n_bg,
        pval_cut=None,
        pos_only=True,
        multimapping=False,
    )
    logger.info(
        f"Found {len(dorc_all)} peak-gene pairs; "
        f"{dorc_all['Gene'].nunique()} genes, {dorc_all['PeakRanges'].nunique()} peaks."
    )

    dorc_all_path = os.path.join(outdir, "dorc_all.tsv.gz")
    dorc_all.to_csv(dorc_all_path, sep="\t", index=False)

    # 4. Filter by pvalZ
    logger.info(f"Filtering peak-gene pairs by pvalZ <= {args.pval_cutoff}...")
    dorc_sig = dorc_all[dorc_all["pvalZ"] <= args.pval_cutoff].copy()
    logger.info(f"{len(dorc_sig)} significant peak-gene pairs after pvalZ filtering.")

    dorc_sig_path = os.path.join(outdir, "dorc_significant.tsv.gz")
    dorc_sig.to_csv(dorc_sig_path, sep="\t", index=False)

    # 5. J-plot to get DORC gene list
    logger.info("Generating J-plot and selecting DORC genes...")
    plt.figure(figsize=(6, 5))
    gene_list = scp.dorc.dorc_j_plot(
        dorc_sig,
        cutoff=7,
        label_top=25,
        return_gene_list=True,
    )
    jplot_pdf = os.path.join(outdir, "dorc_jplot.pdf")
    jplot_png = os.path.join(outdir, "dorc_jplot.png")
    plt.tight_layout()
    plt.savefig(jplot_pdf)
    plt.savefig(jplot_png, dpi=300)
    plt.close()

    genes_path = os.path.join(outdir, "dorc_jplot_genes.txt")
    with open(genes_path, "w") as f:
        for g in gene_list:
            f.write(f"{g}\n")
    logger.info(
        "DORC J-plot saved to %s and %s; %d DORC genes selected.",
        jplot_pdf, jplot_png, len(gene_list)
    )

    # 6. Compute DORC scores per cell (global, across all cells)
    logger.info("Computing DORC scores per cell...")
    adata_dorc = scp.dorc.get_dorc_score(
        adata_atac,
        dorc_sig,
        normalize_atac=True,
        gene_list=gene_list,
    )

    # Preserve UMAP if available in ATAC
    if "X_umap" in adata_atac.obsm:
        adata_dorc.obsm["X_umap"] = adata_atac.obsm["X_umap"]

    dorc_scores_path = os.path.join(outdir, "dorc_scores.h5ad")
    adata_dorc.write_h5ad(dorc_scores_path)

    # 6b. UMAP plots of DORC scores for top genes
    if "X_umap" in adata_dorc.obsm:
        logger.info("Generating UMAP plots for top DORC genes...")
        umap_dir = os.path.join(outdir, "dorc_umap_plots")
        ensure_dir(umap_dir)

        n_top = min(10, len(gene_list))
        top_genes = [g for g in gene_list if g in adata_dorc.var_names][:n_top]

        if len(top_genes) == 0:
            logger.warning(
                "No DORC genes from J-plot found in adata_dorc.var_names; skipping UMAP plots."
            )
        else:
            logger.info(
                "Plotting UMAP for %d DORC genes: %s",
                len(top_genes),
                ", ".join(top_genes),
            )
            for g in top_genes:
                plt.figure(figsize=(5, 4))
                sc.pl.umap(
                    adata_dorc,
                    color=g,
                    show=False,
                    title=f"DORC score: {g}",
                )
                out_pdf = os.path.join(umap_dir, f"umap_dorc_{g}.pdf")
                out_png = os.path.join(umap_dir, f"umap_dorc_{g}.png")
                plt.tight_layout()
                plt.savefig(out_pdf)
                plt.savefig(out_png, dpi=300)
                plt.close()
    else:
        logger.warning(
            "X_umap not present in adata_dorc.obsm; skipping UMAP plots."
        )

    # 7. Write log summary
    log_txt = os.path.join(outdir, "dorc_logs.txt")
    with open(log_txt, "w") as f:
        f.write(f"Total peak-gene pairs: {len(dorc_all)}\n")
        f.write(f"Significant (pvalZ<={args.pval_cutoff}): {len(dorc_sig)}\n")
        f.write(f"Number of DORC genes (J-plot): {len(gene_list)}\n")

    logger.info("scPrinter DORC analysis complete.")

if __name__ == "__main__":
    main()
