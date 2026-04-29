#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import anndata as ad
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peak_matrix", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--condition_key", default="condition")
    ap.add_argument("--cell_type_key", default="cell_type")
    ap.add_argument("--sample_key", default="sample")
    ap.add_argument("--output_dir", default=".")
    args = ap.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    pm = ad.read_h5ad(args.peak_matrix)
    obs = pm.obs.copy()

    # Write cell types list (if available)
    cell_types = []
    if args.cell_type_key in obs.columns:
        cell_types = sorted([str(x) for x in obs[args.cell_type_key].dropna().unique()])
    (out / "cell_types.txt").write_text("\n".join(cell_types) + ("\n" if cell_types else ""))

    # Basic counts
    def safe_counts(col):
        if col not in obs.columns:
            return pd.DataFrame()
        return obs[col].astype(str).value_counts(dropna=False).rename_axis(col).reset_index(name="n_cells")

    safe_counts(args.sample_key).to_csv(out / "cells_by_sample.csv", index=False)
    safe_counts(args.condition_key).to_csv(out / "cells_by_condition.csv", index=False)
    safe_counts(args.cell_type_key).to_csv(out / "cells_by_cell_type.csv", index=False)

    # cell_type x condition
    if args.cell_type_key in obs.columns and args.condition_key in obs.columns:
        ct_cond = (
            obs[[args.cell_type_key, args.condition_key]]
            .astype(str)
            .value_counts()
            .reset_index(name="n_cells")
        )
        ct_cond.to_csv(out / "cells_by_celltype_and_condition.csv", index=False)

    summary = {
        "n_cells": int(pm.n_obs),
        "n_peaks": int(pm.n_vars),
        "has_cell_type": args.cell_type_key in obs.columns,
        "has_condition": args.condition_key in obs.columns,
        "n_cell_types": len(cell_types),
        "cell_types": cell_types[:50],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

if __name__ == "__main__":
    main()
