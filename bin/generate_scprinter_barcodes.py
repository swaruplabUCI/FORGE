#!/usr/bin/env python3
"""
generate_scprinter_barcodes.py

Extract barcode whitelists from SnapATAC QC AnnData files.

For each (sample_name, qc_anndata.h5ad):

  - Reads obs_names as the set of QC-passed barcodes
  - Optionally strips leading 'sampleX_' prefixes when the remainder looks
    like a 10x-style barcode (ACGTN...-1) — using your existing logic
  - Writes <sample_name>_barcodes.txt with one barcode per line.

This script is multi-sample: you provide parallel lists of --qc-anndatas
and --sample-names.
"""

import argparse
import re
from pathlib import Path

import anndata as ad
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate scPrinter barcode whitelists from QC AnnData files"
    )
    p.add_argument(
        "--qc-anndatas",
        nargs="+",
        required=True,
        help="List of QC AnnData .h5ad files (one per sample)",
    )
    p.add_argument(
        "--sample-names",
        nargs="+",
        required=True,
        help="List of sample names matching the QC AnnData files",
    )
    p.add_argument(
        "--barcode-dir",
        required=True,
        help="Output directory to write <sample>_barcodes.txt files",
    )
    p.add_argument(
        "--strip-sample-prefix",
        action="store_true",
        default=True,
        help="Strip leading 'sampleX_' prefix when remainder looks like a 10x barcode "
             "(ACGTN...-1). Enabled by default.",
    )
    return p.parse_args()


# Your existing pattern and logic
barcode_pattern = re.compile(r"^[ACGTN]+-\d+$")


def strip_sample_prefix_if_any(bc: str) -> str:
    """
    If barcode looks like 'sampleX_ACGT...-1' and tail is a 10x-style barcode,
    return 'ACGT...-1'; otherwise return bc unchanged.
    """
    if "_" in bc:
        head, tail = bc.split("_", 1)
        if barcode_pattern.match(tail):
            return tail
    return bc


def main():
    args = parse_args()

    qc_paths = [Path(p) for p in args.qc_anndatas]
    sample_names = list(args.sample_names)
    barcode_dir = Path(args.barcode_dir)
    barcode_dir.mkdir(parents=True, exist_ok=True)

    if len(qc_paths) != len(sample_names):
        raise ValueError(
            f"--qc-anndatas ({len(qc_paths)}) and --sample-names ({len(sample_names)}) "
            "must have the same length."
        )

    print("QC AnnData files:")
    for p, s in zip(qc_paths, sample_names):
        print(f"  {s}: {p}")

    print(f"\nBarcode output directory: {barcode_dir}")
    print(f"Strip sample prefix: {args.strip_sample_prefix}")

    for qc_path, sample in zip(qc_paths, sample_names):
        if not qc_path.exists():
            raise FileNotFoundError(f"QC AnnData not found: {qc_path}")

        print(f"\nProcessing sample {sample} from {qc_path}")
        # Open read-only (backed mode) to avoid loading data into memory
        adata = ad.read_h5ad(qc_path, backed="r")

        # obs_names are the QC-passed barcodes
        barcodes = pd.Index(adata.obs_names).astype(str)

        if args.strip_sample_prefix:
            barcodes_clean = barcodes.to_series().map(strip_sample_prefix_if_any)
        else:
            barcodes_clean = pd.Series(barcodes)

        out_path = barcode_dir / f"{sample}_barcodes.txt"
        barcodes_clean.to_csv(out_path, index=False, header=False)
        print(f"  Wrote {len(barcodes_clean)} barcodes to {out_path}")

        # Close backed AnnData
        adata.file.close()

    print("\nAll barcode whitelist files generated.")


if __name__ == "__main__":
    main()
