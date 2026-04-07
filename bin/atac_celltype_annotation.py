#!/usr/bin/env python3
"""
Automated cell type calling based on marker gene activity scores.
Supports tissue-specific marker sets (PBMC, brain, general) loaded
from an external JSON config, with optional user-provided overrides.
"""

import argparse
import pandas as pd
import numpy as np
import json
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster_scores", required=True,
                    help="CSV of cluster-averaged gene scores")
    ap.add_argument("--species", default="human", choices=["mouse", "human"])
    ap.add_argument("--tissue_type", default="pbmc",
                    help="Tissue type for marker selection (pbmc, brain, general)")
    ap.add_argument("--marker_file", default=None,
                    help="Path to marker_genes.json config (default: bundled markers)")
    ap.add_argument("--resolution", default="0_5",
                    help="Clustering resolution key (e.g., 0_5 for leiden_0_5)")
    ap.add_argument("--out", required=True,
                    help="Output JSON with cluster -> cell type mapping")
    return ap.parse_args()



def get_marker_genes(species, tissue_type, marker_file=None):
    """
    Load marker genes for the given species and tissue type.
    Requires either a user-provided marker_file or a bundled JSON in configs/.
    """
    # Option 1: User-provided marker file (full override)
    if marker_file and Path(marker_file).exists():
        print(f"Loading markers from user-provided file: {marker_file}")
        with open(marker_file) as f:
            all_markers = json.load(f)
        # User file can be flat (just cell_type: [genes]) or nested (species/tissue)
        if species in all_markers and isinstance(all_markers[species], dict):
            if tissue_type in all_markers[species]:
                return all_markers[species][tissue_type]
            else:
                print(f"  WARNING: tissue_type '{tissue_type}' not in file, using first available")
                first_key = list(all_markers[species].keys())[0]
                return all_markers[species][first_key]
        else:
            # Assume flat structure: {cell_type: [genes]}
            return all_markers

    # Option 2: Bundled JSON in configs/ directory (look relative to script)
    script_dir = Path(__file__).resolve().parent
    bundled_json = script_dir.parent / "configs" / "marker_genes.json"
    if bundled_json.exists():
        print(f"Loading markers from bundled config: {bundled_json}")
        with open(bundled_json) as f:
            all_markers = json.load(f)
        if species in all_markers and tissue_type in all_markers[species]:
            markers = all_markers[species][tissue_type]
            print(f"  Using {species}/{tissue_type}: {len(markers)} cell types")
            return markers
        elif species in all_markers:
            available = list(all_markers[species].keys())
            print(f"  WARNING: tissue_type '{tissue_type}' not found. Available: {available}")
            print(f"  Falling back to '{available[0]}'")
            return all_markers[species][available[0]]

    raise ValueError(
        f"No marker genes available for species={species}, tissue_type={tissue_type}. "
        "Provide a --marker_file JSON or place marker_genes.json in configs/."
    )


def score_cell_type(cluster_expr, markers, genes_available):
    """
    Score a cluster for a cell type based on mean expression of markers.
    Returns normalized score (0-1 range).
    """
    available_markers = [g for g in markers if g in genes_available]

    if len(available_markers) == 0:
        return 0.0

    scores = [cluster_expr[g] for g in available_markers if g in cluster_expr.index]

    if len(scores) == 0:
        return 0.0

    return float(np.mean(scores))


def annotate_clusters(cluster_scores_df, marker_dict):
    """Assign cell types to clusters based on marker expression."""
    annotations = {}
    genes_available = set(cluster_scores_df.columns)

    for cluster_id in cluster_scores_df.index:
        cluster_expr = cluster_scores_df.loc[cluster_id]

        scores = {}
        for cell_type, markers in marker_dict.items():
            scores[cell_type] = score_cell_type(cluster_expr, markers, genes_available)

        if max(scores.values()) > 0:
            best_type = max(scores, key=scores.get)
            confidence = scores[best_type]

            annotations[str(cluster_id)] = {
                "cell_type": best_type,
                "confidence": confidence,
                "scores": {k: round(v, 4) for k, v in scores.items()},
            }
        else:
            annotations[str(cluster_id)] = {
                "cell_type": "Unknown",
                "confidence": 0.0,
                "scores": scores,
            }

    return annotations


def main():
    args = parse_args()

    # Load cluster scores
    cluster_scores = pd.read_csv(args.cluster_scores, index_col=0)
    print(f"Loaded scores for {len(cluster_scores)} clusters x {len(cluster_scores.columns)} genes")

    # Get marker genes (tissue-aware)
    markers = get_marker_genes(args.species, args.tissue_type, args.marker_file)
    print(f"Using {len(markers)} cell type definitions for {args.species}/{args.tissue_type}")
    for ct, genes in markers.items():
        available = [g for g in genes if g in cluster_scores.columns]
        print(f"  {ct}: {len(available)}/{len(genes)} markers found in data")

    # Annotate
    annotations = annotate_clusters(cluster_scores, markers)

    # Write output
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w") as f:
        json.dump(annotations, f, indent=2)

    print(f"\nCell type annotations saved to: {out_path}")
    print("\nSummary:")
    for cluster_id, info in annotations.items():
        print(f"  Cluster {cluster_id}: {info['cell_type']} (confidence: {info['confidence']:.3f})")


if __name__ == "__main__":
    main()
