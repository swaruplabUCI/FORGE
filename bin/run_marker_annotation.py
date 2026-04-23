#!/usr/bin/env python3
"""
Marker-gene-based cell type annotation for RNA data.

Used when no CellTypist model exists for the tissue (e.g., mouse kidney).
Scores each cell against curated per-cell-type marker gene lists with
sc.tl.score_genes, then assigns labels by argmax under min-score and
top-vs-second-best margin gates. Cells failing both gates are labeled
'unknown'. Near-duplicate subtypes (e.g., PT subtypes) share a
`collapse_group` in the marker CSV; within-group margin failures fall
back to the parent label rather than 'unknown'.

Input CSV schema (long format):
    cell_type,gene,rank,source_cluster,collapse_group
    Proximal Tubule,Slc34a1,1,PT,PT
    ...

Output obs columns on the annotated h5ad:
    cell_type_marker          -- primary string label (consumed by
                                 build_mudata_batched.py precedence chain)
    cell_type_marker_score    -- top score (float)
    score_<cell_type>         -- per-set score (one column per marker set)

Usage:
    python run_marker_annotation.py --input concat.h5ad --markers markers.csv \\
        --output marker_annotated.h5ad --min_score 0.0 --min_margin 0.1
"""

import argparse
import logging
import sys

import numpy as np
import pandas as pd
import scanpy as sc

from h5ad_compat import sanitize_adata

logger = logging.getLogger(__name__)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Marker-gene cell type annotation on RNA h5ad"
    )
    ap.add_argument("--input", required=True, help="Input h5ad file")
    ap.add_argument("--markers", required=True, help="Marker CSV (see schema in docstring)")
    ap.add_argument("--output", required=True, help="Output annotated h5ad file")
    ap.add_argument("--min_score", type=float, default=0.0,
                    help="Minimum winner score_genes value; below → 'unknown' (default 0.0)")
    ap.add_argument("--min_margin", type=float, default=0.1,
                    help="Minimum winner − runner-up gap; below → collapse or 'unknown' (default 0.1)")
    ap.add_argument("--ctrl_size", type=int, default=50,
                    help="score_genes control pool size (default 50)")
    return ap.parse_args()


def load_markers(path):
    """Parse marker CSV into ordered dict {cell_type: [genes]} and collapse map."""
    df = pd.read_csv(path)
    required = {"cell_type", "gene"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Marker CSV missing required columns: {missing}")

    if "collapse_group" not in df.columns:
        df["collapse_group"] = ""
    df["collapse_group"] = df["collapse_group"].fillna("").astype(str)

    marker_sets = {}
    collapse_map = {}  # cell_type -> parent label (or "" if none)
    for ct, sub in df.groupby("cell_type", sort=False):
        if "rank" in sub.columns:
            sub = sub.sort_values("rank", kind="stable")
        genes = sub["gene"].dropna().astype(str).tolist()
        marker_sets[ct] = genes
        group_vals = [g for g in sub["collapse_group"].unique() if g]
        collapse_map[ct] = group_vals[0] if group_vals else ""
    return marker_sets, collapse_map


def prepare_expression(adata):
    """Return a working copy that is log-normalized (score_genes assumption)."""
    adata_sc = adata.copy()
    if "counts" in adata_sc.layers:
        logger.info("Using 'counts' layer; renormalizing for scoring")
        adata_sc.X = adata_sc.layers["counts"].copy()
        sc.pp.normalize_total(adata_sc, target_sum=1e4)
        sc.pp.log1p(adata_sc)
    elif adata_sc.X.max() > 50:
        logger.info("X looks like raw counts; normalizing")
        sc.pp.normalize_total(adata_sc, target_sum=1e4)
        sc.pp.log1p(adata_sc)
    else:
        logger.info("X appears already log-normalized")
    return adata_sc


def score_all_sets(adata_sc, marker_sets, ctrl_size):
    """Score each marker set; return (scores_df cells×celltypes, dropped_summary)."""
    var_set = set(adata_sc.var_names)
    scores = {}
    dropped = {}
    for ct, genes in marker_sets.items():
        present = [g for g in genes if g in var_set]
        dropped_genes = [g for g in genes if g not in var_set]
        if len(present) < 3:
            logger.warning(
                f"Skipping '{ct}': only {len(present)} of {len(genes)} markers present (<3)"
            )
            dropped[ct] = {"kept": len(present), "total": len(genes), "skipped": True}
            continue
        score_key = f"score_{ct}"
        sc.tl.score_genes(
            adata_sc,
            gene_list=present,
            score_name=score_key,
            ctrl_size=min(ctrl_size, max(len(present) * 10, ctrl_size)),
            random_state=0,
        )
        scores[ct] = adata_sc.obs[score_key].to_numpy()
        dropped[ct] = {
            "kept": len(present),
            "total": len(genes),
            "skipped": False,
            "dropped_genes": dropped_genes,
        }
    if not scores:
        raise RuntimeError("No marker sets had enough genes to score.")
    scores_df = pd.DataFrame(scores, index=adata_sc.obs_names)
    return scores_df, dropped


def assign_labels(scores_df, collapse_map, min_score, min_margin):
    """Argmax with min-score and margin gates; collapse-group fallback."""
    cell_types = scores_df.columns.to_numpy()
    vals = scores_df.to_numpy()
    order = np.argsort(-vals, axis=1)  # descending
    top_idx = order[:, 0]
    top_score = vals[np.arange(len(vals)), top_idx]
    if vals.shape[1] >= 2:
        second_idx = order[:, 1]
        second_score = vals[np.arange(len(vals)), second_idx]
    else:
        second_idx = np.full(len(vals), -1)
        second_score = np.full(len(vals), -np.inf)

    labels = []
    for i in range(len(vals)):
        winner = cell_types[top_idx[i]]
        if top_score[i] < min_score:
            labels.append("unknown")
            continue
        if vals.shape[1] < 2 or (top_score[i] - second_score[i]) >= min_margin:
            labels.append(winner)
            continue
        # Margin failed. If runner-up is in the same collapse group, fall back
        # to the group parent label; otherwise 'unknown'.
        runner = cell_types[second_idx[i]] if second_idx[i] >= 0 else None
        w_group = collapse_map.get(winner, "")
        r_group = collapse_map.get(runner, "") if runner else ""
        if w_group and w_group == r_group:
            labels.append(w_group)
        else:
            labels.append("unknown")
    return np.array(labels), top_score


def write_report(path, adata, marker_sets, dropped, labels, top_score,
                 min_score, min_margin, markers_path):
    counts = pd.Series(labels).value_counts()
    with open(path, "w") as f:
        f.write("Marker-Gene Annotation Report\n")
        f.write("=" * 40 + "\n")
        f.write(f"Markers: {markers_path}\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Cell types scored: {sum(1 for d in dropped.values() if not d['skipped'])}\n")
        f.write(f"Cell types skipped (<3 markers present): "
                f"{sum(1 for d in dropped.values() if d['skipped'])}\n")
        f.write(f"min_score: {min_score}  min_margin: {min_margin}\n\n")
        f.write("Label distribution:\n")
        f.write(counts.to_string())
        f.write("\n\nTop-score summary (winner score):\n")
        ts = pd.Series(top_score)
        f.write(ts.describe().to_string())
        f.write("\n\nMarker set coverage (kept/total):\n")
        for ct, d in dropped.items():
            mark = " [SKIPPED]" if d["skipped"] else ""
            f.write(f"  {ct}: {d['kept']}/{d['total']}{mark}\n")
            dg = d.get("dropped_genes", [])
            if dg:
                f.write(f"    dropped: {', '.join(dg)}\n")


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    logger.info(f"Loading h5ad: {args.input}")
    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded {adata.n_obs} cells × {adata.n_vars} genes")

    logger.info(f"Loading markers: {args.markers}")
    marker_sets, collapse_map = load_markers(args.markers)
    logger.info(f"Loaded {len(marker_sets)} marker sets: {list(marker_sets)}")

    adata_sc = prepare_expression(adata)

    scores_df, dropped = score_all_sets(adata_sc, marker_sets, args.ctrl_size)

    labels, top_score = assign_labels(
        scores_df, collapse_map, args.min_score, args.min_margin
    )

    adata.obs["cell_type_marker"] = pd.Categorical(labels)
    adata.obs["cell_type_marker_score"] = top_score
    for ct in scores_df.columns:
        adata.obs[f"score_{ct}"] = scores_df[ct].to_numpy()

    n_unknown = int((labels == "unknown").sum())
    logger.info(
        f"Assigned {adata.n_obs - n_unknown}/{adata.n_obs} cells a label "
        f"({n_unknown} unknown, {n_unknown / max(adata.n_obs, 1):.1%})"
    )
    logger.info(f"Label counts:\n{pd.Series(labels).value_counts().to_string()}")

    report_path = args.output.replace(".h5ad", "_marker_report.txt")
    write_report(
        report_path, adata, marker_sets, dropped, labels, top_score,
        args.min_score, args.min_margin, args.markers,
    )
    logger.info(f"Report: {report_path}")

    logger.info(f"Writing annotated h5ad: {args.output}")
    sanitize_adata(adata, args.output)
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
