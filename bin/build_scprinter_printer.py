#!/usr/bin/env python3
"""
build_scprinter_printer.py

Create a scPrinter "printer" object from fragment files and barcode whitelists
using scprinter.pp.import_fragments, and save it as a .h5ad file.

All paths (fragments, barcodes, sample_names) are provided via CLI arguments.
"""

import os
import argparse
from pathlib import Path

import scprinter as scp


def init_scprinter(cache_dir: str):
    """
    Initialize scprinter the same way you do in your notebooks:
      - set datasets.cache path
    """
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp


def parse_args():
    p = argparse.ArgumentParser(
        description="Build scPrinter printer object from fragments + barcode whitelists"
    )

    p.add_argument(
        "--cache-dir",
        required=True,
        help="scPrinter cache directory (SCPRINTER_CACHE_DIR)",
    )

    p.add_argument(
        "--output",
        required=True,
        help="Output .h5ad path for the printer",
    )

    p.add_argument(
        "--genome",
        required=True,
        help="Genome key for scPrinter (e.g. 'mm10')",
    )

    p.add_argument(
        "--fragment-files",
        nargs="+",
        required=True,
        help="List of fragment.tsv.gz files (one per sample)",
    )

    p.add_argument(
        "--barcode-lists",
        nargs="+",
        required=True,
        help="List of barcode whitelist text files (one per sample)",
    )

    p.add_argument(
        "--sample-names",
        nargs="+",
        required=True,
        help="Sample names corresponding to fragment & barcode lists",
    )

    p.add_argument(
        "--min-num-fragments",
        type=int,
        default=0,
        help="Minimum number of fragments per cell [default: %(default)s]",
    )

    p.add_argument(
        "--min-tsse",
        type=float,
        default=0.0,
        help="Minimum TSSE threshold per cell [default: %(default)s]",
    )

    p.add_argument(
        "--sorted-by-barcode",
        action="store_true",
        help="Set if fragments are already sorted by barcode",
    )

    p.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of parallel jobs for import_fragments [default: %(default)s]",
    )

    return p.parse_args()


def main():
    args = parse_args()

    fragment_files = [Path(f) for f in args.fragment_files]
    barcode_lists = [Path(b) for b in args.barcode_lists]
    sample_names = list(args.sample_names)

    if not (len(fragment_files) == len(barcode_lists) == len(sample_names)):
        raise ValueError(
            f"fragment_files ({len(fragment_files)}), barcode_lists ({len(barcode_lists)}), "
            f"and sample_names ({len(sample_names)}) must have the same length."
        )

    print("Using fragment files:")
    for f in fragment_files:
        print("  ", f)

    print("\nUsing barcode lists:")
    for b in barcode_lists:
        print("  ", b)

    print("\nSample names:")
    print("  ", sample_names)

    # Check that all input files exist
    for f in fragment_files:
        if not f.exists():
            raise FileNotFoundError(f"Fragment file not found: {f}")
    for b in barcode_lists:
        if not b.exists():
            raise FileNotFoundError(f"Barcode list not found: {b}")

    # Configure environment for scPrinter
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir
    print(f"\nSCPRINTER_CACHE_DIR set to: {args.cache_dir}")

    # Initialize scprinter in the standard way
    global scp
    scp = init_scprinter(args.cache_dir)

    # Resolve genome object
    try:
        genome_obj = getattr(scp.genome, args.genome)
    except AttributeError:
        raise ValueError(f"Genome '{args.genome}' not found in scprinter.genome")

    # Prepare output path
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nBuilding printer with scprinter.pp.import_fragments")
    print(f"Output printer will be saved to: {out_path}")

    printer = scp.pp.import_fragments(
        path_to_frags=[str(f) for f in fragment_files],
        barcodes=[str(b) for b in barcode_lists],
        savename=str(out_path),
        genome=genome_obj,
        sample_names=sample_names,
        min_num_fragments=args.min_num_fragments,
        min_tsse=args.min_tsse,
        sorted_by_barcode=args.sorted_by_barcode,
        n_jobs=args.n_jobs,
    )

    print("\nPrinter created successfully.")
    print("Printer summary:")
    print(printer)
    print("\nDone.")


if __name__ == "__main__":
    main()
