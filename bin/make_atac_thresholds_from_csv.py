#!/usr/bin/env python3
import argparse
import json
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", required=True,
                        help="sample_statistics.csv from ATAC_INITIAL_QC")
    parser.add_argument("--min_cells", type=int, default=100,
                        help="Minimum cells required for strict thresholds")
    parser.add_argument("--output", required=True,
                        help="Output JSON file (e.g. sampthreshconfig.json)")
    args = parser.parse_args()

    df = pd.read_csv(args.stats)

    # FIX-F5: Validate expected columns before processing
    required_cols = ['sample', 'fragments_q05', 'fragments_q95', 'tsse_q05',
                     'fragments_min', 'fragments_max', 'tsse_min',
                     'n_cells_after_frag_5_95_and_tsse_q05']
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"sample_statistics.csv missing required columns: {missing_cols}\n"
            f"Expected: {required_cols}\n"
            f"Found:    {list(df.columns)}"
        )

    thresholds = {}
    for row in df.itertuples():
        sample = row.sample

        # Proposed thresholds
        min_counts = row.fragments_q05
        max_counts = row.fragments_q95
        min_tsse   = row.tsse_q05

        if row.n_cells_after_frag_5_95_and_tsse_q05 < args.min_cells:
            # Too few cells would survive: relax to keep all cells.
            min_counts = row.fragments_min
            max_counts = row.fragments_max
            min_tsse   = row.tsse_min

        thresholds[sample] = {
            "min_counts": float(min_counts),
            "max_counts": float(max_counts),
            "min_tsse": float(min_tsse),
        }

    with open(args.output, "w") as f:
        json.dump(thresholds, f, indent=2)

    print(f"Wrote thresholds for {len(thresholds)} samples to {args.output}")

if __name__ == "__main__":
    main()
