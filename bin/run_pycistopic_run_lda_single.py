#!/usr/bin/env python
"""
run_pycistopic_run_lda_single.py — Phase 3b of the pyCisTopic fan-out pipeline.

Loads the merged CistopicObject and runs MALLET LDA for exactly one topic
count, saving the resulting CistopicLDAModel as Topic{n}.pkl in --outdir.

Intended to be invoked as a parallel fan-out: one SLURM job per topic count.
Each job is independent and can be retried individually without repeating
the merge or other models.
"""

import argparse
import logging
import os
import pickle
import subprocess
import sys

from pycisTopic.lda_models import run_cgs_models_mallet


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycistopic_run_lda_single")


def parse_args():
    p = argparse.ArgumentParser(
        description="pyCisTopic Phase 3b: run single MALLET LDA model."
    )
    p.add_argument("--merged-pkl",    required=True,
                   help="Path to merged_cistopic.pkl from Phase 3a.")
    p.add_argument("--n-topics",      type=int, required=True,
                   help="Number of topics for this model.")
    p.add_argument("--n-cpu",         type=int, default=8,
                   help="MALLET --num-threads (per-model parallelism).")
    p.add_argument("--n-iter",        type=int, default=500,
                   help="MALLET --num-iterations.")
    p.add_argument("--random-state",  type=int, default=555)
    p.add_argument("--alpha",         type=float, default=50.0)
    p.add_argument("--eta",           type=float, default=0.1)
    p.add_argument("--mallet-path",   default="Mallet-202108/bin/mallet")
    p.add_argument("--mallet-memory", default="100G",
                   help="Java heap for MALLET (MALLET_MEMORY env var).")
    p.add_argument("--outdir",        default=".",
                   help="Output directory; Topic{n}.pkl written here.")
    return p.parse_args()


def ensure_mallet(mallet_path, logger):
    if os.path.exists(mallet_path):
        return mallet_path
    logger.info("MALLET not found at %s — downloading...", mallet_path)
    subprocess.check_call([
        "wget", "-q",
        "https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz",
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

    out_pkl = os.path.join(args.outdir, f"Topic{args.n_topics}.pkl")

    logger.info("Loading merged_cistopic.pkl from %s", args.merged_pkl)
    with open(args.merged_pkl, "rb") as f:
        cistopic_obj = pickle.load(f)
    logger.info("Loaded: %d cells, %d regions",
                len(cistopic_obj.cell_names), len(cistopic_obj.region_names))

    mallet_path = ensure_mallet(args.mallet_path, logger)
    os.environ["MALLET_MEMORY"] = args.mallet_memory

    mallet_tmp = os.path.join(args.outdir, "mallet_tmp")
    os.makedirs(mallet_tmp, exist_ok=True)

    logger.info("Running MALLET LDA: %d topics, %d iter, %d threads",
                args.n_topics, args.n_iter, args.n_cpu)
    run_cgs_models_mallet(
        cistopic_obj,
        n_topics=[args.n_topics],
        n_cpu=args.n_cpu,
        n_iter=args.n_iter,
        random_state=args.random_state,
        alpha=args.alpha,
        alpha_by_topic=True,
        eta=args.eta,
        eta_by_topic=False,
        tmp_path=mallet_tmp,
        save_path=args.outdir,
        mallet_path=mallet_path,
    )

    if not os.path.exists(out_pkl):
        raise RuntimeError(f"Expected output {out_pkl} not found after MALLET run.")
    logger.info("Topic%d model saved to %s", args.n_topics, out_pkl)


if __name__ == "__main__":
    main()
