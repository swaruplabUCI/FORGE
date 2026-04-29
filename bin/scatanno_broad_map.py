#!/usr/bin/env python3
"""
scatanno_broad_map.py

Builds SCATANNO_BROAD_MAP : dict[subclass_label -> class_label] from the Allen
Brain Atlas (mouse) subclass→class taxonomy CSV.

Used by merge_annotations.py to populate 'cell_type_broad' on peak_matrix
after scATAnno reference-based annotation. Mirrors the celltypist_broad_map.py
pattern used on the RNA side.

The CSV covers 1005 Allen subclasses collapsing to 35 classes (e.g.
  'Oligo NN' -> 'OPC-Oligo'
  'Astro-TE NN' -> 'Astro-Epen'
  'Microglia NN' -> 'Immune'
  'L2/3 IT CTX Glut' -> 'IT-ET Glut'
).

Empty class_label (non-brain / missing) rows are skipped; downstream code
should fill missing predictions with 'Unknown'.
"""

import csv
from pathlib import Path

ALLEN_SUBCLASS_CSV = Path(
    "/dfs7/swaruplab/lesolano/ref/ALLEN_mouse_10xv3/subclassanno.3.csv"
)


def load_broad_map(csv_path: Path = ALLEN_SUBCLASS_CSV) -> dict:
    """Return {subclass_label: class_label} for all rows with a non-empty class_label."""
    mapping: dict = {}
    if not Path(csv_path).exists():
        raise FileNotFoundError(
            f"Allen subclass→class CSV not found at {csv_path}. "
            "Set SCATANNO_BROAD_MAP_CSV env var or pass csv_path explicitly."
        )
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            sub = row.get("subclass_label", "").strip()
            cls = row.get("class_label", "").strip()
            if sub and cls:
                mapping[sub] = cls
    return mapping


# Module-level singleton — import this from merge_annotations.py
SCATANNO_BROAD_MAP = load_broad_map()


if __name__ == "__main__":
    print(f"Loaded {len(SCATANNO_BROAD_MAP)} subclass -> class mappings")
    classes = sorted(set(SCATANNO_BROAD_MAP.values()))
    print(f"Unique broad classes: {len(classes)}")
    for c in classes:
        n = sum(1 for v in SCATANNO_BROAD_MAP.values() if v == c)
        print(f"  [{n:4d}] {c}")
