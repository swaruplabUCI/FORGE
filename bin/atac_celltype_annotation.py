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


# ---- Inline fallback markers (used when no JSON config is available) ----

FALLBACK_MARKERS = {
    "human": {
        "pbmc": {
            "CD4_T_cells": ["CD3D", "CD3E", "CD4", "IL7R", "LEF1"],
            "CD8_T_cells": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMK"],
            "NK_cells": ["NKG7", "GNLY", "KLRD1", "NCAM1", "KLRF1"],
            "B_cells": ["CD79A", "MS4A1", "CD19", "PAX5", "BANK1"],
            "Plasma_cells": ["JCHAIN", "MZB1", "IGHG1", "SDC1", "XBP1"],
            "CD14_Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
            "CD16_Monocytes": ["FCGR3A", "MS4A7", "ITGAL", "LST1"],
            "Dendritic_cells": ["FCER1A", "CST3", "CLEC10A", "CD1C"],
            "Plasmacytoid_DC": ["IL3RA", "CLEC4C", "IRF7", "TCF4"],
            "Platelets": ["PPBP", "PF4", "GP9", "TUBB1"],
        },
        "brain": {
            "Oligodendrocytes": ["PLP1", "MBP", "MOG", "MOBP"],
            "OPC": ["PDGFRA", "CSPG4", "SOX10", "OLIG1"],
            "Astrocytes": ["GFAP", "AQP4", "SLC1A2", "ALDH1L1"],
            "Microglia": ["CX3CR1", "TMEM119", "P2RY12", "CSF1R"],
            "Excitatory_neurons": ["SLC17A7", "CAMK2A", "NEUROD6"],
            "Inhibitory_neurons": ["GAD1", "GAD2", "SLC32A1"],
            "Endothelial": ["CLDN5", "FLT1", "PECAM1"],
            "Pericytes": ["PDGFRB", "RGS5", "ABCC9"],
        },
    },
    "mouse": {
        "pbmc": {
            "CD4_T_cells": ["Cd3d", "Cd3e", "Cd4", "Il7r", "Lef1"],
            "CD8_T_cells": ["Cd3d", "Cd3e", "Cd8a", "Cd8b1", "Gzmk"],
            "NK_cells": ["Nkg7", "Gzma", "Klrb1c", "Ncam1", "Klrd1"],
            "B_cells": ["Cd79a", "Ms4a1", "Cd19", "Pax5", "Bank1"],
            "Plasma_cells": ["Jchain", "Mzb1", "Ighg1", "Sdc1", "Xbp1"],
            "CD14_Monocytes": ["Cd14", "Lyz2", "S100a8", "S100a9", "Vcan"],
            "CD16_Monocytes": ["Fcgr4", "Ms4a7", "Itgal", "Lst1"],
            "Dendritic_cells": ["Fcer1a", "Cst3", "Cd209a", "Itgax"],
            "Plasmacytoid_DC": ["Siglech", "Bst2", "Irf7", "Tcf4"],
            "Platelets": ["Ppbp", "Pf4", "Gp9", "Tubb1"],
        },
        "brain": {
            "Oligodendrocytes": ["Plp1", "Mbp", "Mog", "Mobp"],
            "OPC": ["Pdgfra", "Cspg4", "Sox10", "Olig1"],
            "Astrocytes": ["Gfap", "Aqp4", "Slc1a2", "Aldh1l1"],
            "Microglia": ["Cx3cr1", "Tmem119", "P2ry12", "Csf1r"],
            "Excitatory_neurons": ["Slc17a7", "Camk2a", "Neurod6"],
            "Inhibitory_neurons": ["Gad1", "Gad2", "Slc32a1"],
            "Endothelial": ["Cldn5", "Flt1", "Pecam1"],
            "Pericytes": ["Pdgfrb", "Rgs5", "Abcc9"],
        },
    },
}


def get_marker_genes(species, tissue_type, marker_file=None):
    """
    Load marker genes for the given species and tissue type.
    Priority: 1) user-provided marker_file, 2) bundled JSON, 3) inline fallback.
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

    # Option 3: Inline fallback
    print(f"Using inline fallback markers for {species}/{tissue_type}")
    if species in FALLBACK_MARKERS and tissue_type in FALLBACK_MARKERS[species]:
        return FALLBACK_MARKERS[species][tissue_type]
    elif species in FALLBACK_MARKERS:
        first_key = list(FALLBACK_MARKERS[species].keys())[0]
        print(f"  WARNING: tissue_type '{tissue_type}' not in fallback, using '{first_key}'")
        return FALLBACK_MARKERS[species][first_key]
    else:
        raise ValueError(f"No markers available for species={species}, tissue_type={tissue_type}")


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
