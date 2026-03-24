#!/usr/bin/env python3
"""
run_scprinter_motif_scan.py

- Loads DA peaks for a given cell type.
- Filters significant peaks by FDR.
- Normalizes barcodes and builds control/treatment groupings from the peak matrix.
- Loads an existing scPrinter printer object from --printer-path
  (created via scprinter.pp.import_fragments).
- Computes TF binding scores for all significant peaks using scp.tl.get_binding_score.
- Writes:
    - motif_enrichment_<cell_type>.csv
    - tfbs_<cell_type>.h5ad
    - motif_plots/tf_enrichment.png
"""

import os
import builtins

# ---------------------------------------------------------------------------
# Work around Python 3 / Cython 'long' issue in sorted_nearest
# ---------------------------------------------------------------------------
if not hasattr(builtins, "long"):
    builtins.long = int

# Disable Numba JIT before importing scprinter (matches shared scripts)
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import argparse
from pathlib import Path

import scprinter as scp
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def init_scprinter(cache_dir: str):
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp


def parse_args():
    p = argparse.ArgumentParser(description="Run scPrinter motif scan")
    p.add_argument(
        "--peak-matrix",
        required=True,
        help="Path to filtered peak matrix .h5ad (e.g. filtered_peak_matrix.h5ad)",
    )
    p.add_argument(
        "--da-peaks",
        required=True,
        help="Path to DA peaks CSV for this cell type (e.g. DA_peaks_<cell_type>.csv)",
    )
    p.add_argument(
        "--cell-type",
        required=True,
        help="Cell type name to subset from adata.obs['cell_type'] (for logging)",
    )
    p.add_argument(
        "--cache-dir",
        required=True,
        help="Path to scPrinter cache directory (SCPRINTER_CACHE_DIR)",
    )
    p.add_argument(
        "--genome",
        required=True,
        help="Genome key for scPrinter (e.g. 'mm10')",
    )
    p.add_argument(
        "--fdr",
        type=float,
        required=True,
        help="FDR threshold to define significant peaks",
    )
    p.add_argument(
        "--printer-path",
        required=True,
        help="Path to existing scPrinter printer .h5ad (built by import_fragments)",
    )
    # ADD THESE TWO ARGUMENTS:
    p.add_argument(
        "--control-condition",
        required=True,
        help="Control condition label (e.g. '5xFAD')",
    )
    p.add_argument(
        "--treatment-condition",
        required=True,
        help="Treatment condition label (e.g. 'SREBF1_OE')",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=1,
        help="Number of jobs/threads for TFBS computation",
    )
    return p.parse_args()

def normalize_peak_barcode(bc: str) -> str:
    """Normalize peak-matrix barcode to match printer.obs_names.

    FIX-37: Peak matrix uses colon delimiter (e.g. 'L1_Donor_0_july:1000602')
    while the printer uses underscore (e.g. 'L1_Donor_0_july_1000602').
    Simply replace the first colon with underscore.

    The old logic (`if "_" not in head`) was always False because sample names
    like 'L1_Donor_0_july' contain underscores.
    """
    if ":" in bc:
        return bc.replace(":", "_", 1)
    return bc


def tail_only(bc: str) -> str:
    """Strip sample prefix, keeping only the cell barcode portion."""
    if "_" in bc:
        return bc.split("_", 1)[-1]
    if ":" in bc:
        return bc.split(":", 1)[-1]
    return bc


def choose_normalization(peak_barcodes, printer_barcodes: set) -> str:
    """FIX-42: Auto-detect barcode normalization mode.

    Same logic as run_scprinter_footprinting.py — probes head_tail vs tail_only
    on a sample of barcodes and returns the strategy with more matches.
    """
    sample = peak_barcodes[:min(500, len(peak_barcodes))]
    n_head_tail = sum(1 for bc in sample if normalize_peak_barcode(bc) in printer_barcodes)
    n_tail_only = sum(1 for bc in sample if tail_only(bc) in printer_barcodes)
    chosen = "tail_only" if n_tail_only > n_head_tail else "head_tail"
    print(f"  FIX-42 choose_normalization: head_tail={n_head_tail}, tail_only={n_tail_only} → '{chosen}'")
    return chosen


def apply_normalization(bc: str, mode: str) -> str:
    """Apply the chosen barcode normalization strategy."""
    if mode == "tail_only":
        return tail_only(bc)
    return normalize_peak_barcode(bc)


def main():
    args = parse_args()

    # Configure scPrinter cache
    os.environ["SCPRINTER_CACHE_DIR"] = args.cache_dir

    # Initialize scprinter using your standard pattern
    global scp
    scp = init_scprinter(args.cache_dir)

    # Resolve genome object from string
    try:
        genome_obj = getattr(scp.genome, args.genome)
    except AttributeError:
        raise ValueError(f"Genome '{args.genome}' not found in scprinter.genome")

    # -------------------------------------------------------------------------
    # Load DA peaks and filter by FDR
    # -------------------------------------------------------------------------
    print(f"Loading DA peaks from {args.da_peaks}")
    da_df = pd.read_csv(args.da_peaks)

    if "padj" not in da_df.columns:
        raise ValueError("DA peaks file does not contain 'padj' column")

    if "peak" not in da_df.columns:
        raise ValueError("DA peaks file does not contain 'peak' column")

    sig_peaks = da_df[da_df["padj"] < args.fdr].copy()
    print(f"Total peaks: {len(da_df)}, significant (padj < {args.fdr}): {len(sig_peaks)}")

    if sig_peaks.empty:
        print("No significant peaks at given FDR; writing empty enrichment file.")
        enrichment_df = pd.DataFrame(
            columns=["tf", "cell_type", "n_up_peaks", "n_down_peaks"]
        )
        out_csv = f"motif_enrichment_{args.cell_type}.csv"
        enrichment_df.to_csv(out_csv, index=False)
        print(f"Saved empty enrichment file to {out_csv}")
        return

    # -------------------------------------------------------------------------
    # Convert peak names to regions (Chromosome, Start, End)
    # -------------------------------------------------------------------------
    print("Converting peak names to genomic regions")
    regions = []
    for peak in sig_peaks["peak"]:
        # Expecting something like "chrX:start-end"
        parts = peak.replace(":", "-").split("-")
        if len(parts) != 3:
            print(f"  Warning: could not parse peak '{peak}' -> skipping")
            continue
        try:
            regions.append(
                {
                    "Chromosome": parts[0],
                    "Start": int(parts[1]),
                    "End": int(parts[2]),
                }
            )
        except ValueError:
            print(f"  Warning: non-integer coordinates for peak '{peak}' -> skipping")

    if not regions:
        raise ValueError("No valid regions could be parsed from peaks")

    regions_df = pd.DataFrame(regions)
    print(f"Parsed {len(regions_df)} regions from significant peaks")

    # -------------------------------------------------------------------------
    # Load peak matrix and build cell-type-specific groups
    # -------------------------------------------------------------------------
    print(f"Loading peak matrix from {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)

    if "cell_type" not in adata.obs.columns:
        raise ValueError("adata.obs does not contain a 'cell_type' column")

    adata_ct = adata[adata.obs["cell_type"] == args.cell_type].copy()
    n_cells = adata_ct.n_obs
    print(f"Found {n_cells} cells for cell_type = {args.cell_type}")
    if n_cells == 0:
        raise ValueError(f"No cells found for cell_type '{args.cell_type}'")

    if "condition" not in adata_ct.obs.columns:
        raise ValueError("adata_ct.obs does not contain a 'condition' column")

    conditions = adata_ct.obs["condition"].astype(str).unique()
    print(f"Conditions found in this cell type: {conditions}")

    control_label = args.control_condition
    treatment_label = args.treatment_condition

    group_labels = [lab for lab in [control_label, treatment_label] if lab in conditions]
    if len(group_labels) != 2:
        raise ValueError(
            f"Expected both {control_label} and {treatment_label} in adata_ct.obs['condition'], "
            f"but found: {list(conditions)}"
        )

    # -------------------------------------------------------------------------
    # Load existing printer
    # -------------------------------------------------------------------------
    printer_path = Path(args.printer_path)
    if not printer_path.exists():
        raise FileNotFoundError(f"Printer file not found: {printer_path}")

    print(f"Loading existing printer from {printer_path}")
    printer = scp.load_printer(str(printer_path), genome_obj)

    # Load TFBS model once
    print("Loading TF binding score model into printer")
    printer.load_bindingscore_model("TF", scp.datasets.pretrained_TFBS_model)

    # Barcodes present in the printer
    printer_barcodes = np.asarray(printer.obs_names).astype(str)
    printer_bcs_set = set(printer_barcodes)

    # FIX-42: Auto-detect barcode normalization mode
    all_bcs_raw = adata_ct.obs_names.astype(str)
    bc_norm_mode = choose_normalization(np.asarray(all_bcs_raw), printer_bcs_set)

    # -------------------------------------------------------------------------
    # Build grouping (control vs treatment) with normalized barcodes
    # -------------------------------------------------------------------------
    grouping = []
    for lab in group_labels:
        mask = adata_ct.obs["condition"].astype(str) == lab
        bcs_ct_raw = adata_ct.obs_names[mask].astype(str)
        bcs_ct = np.array([apply_normalization(b, bc_norm_mode) for b in bcs_ct_raw], dtype=str)

        print(f"\n  [{lab}] example peak-matrix barcodes (raw):        {list(bcs_ct_raw[:3])}")
        print(f"  [{lab}] example peak-matrix barcodes (normalized): {list(bcs_ct[:3])}")
        print(f"  [{lab}] example printer barcodes:                  {list(printer_barcodes[:3])}")

        bcs_common = np.intersect1d(bcs_ct, printer_barcodes)
        n_common = len(bcs_common)
        print(f"  Condition {lab}: {n_common} cells in printer (out of {int(mask.sum())} in peak matrix)")
        if n_common == 0:
            raise ValueError(
                f"No cells for condition '{lab}' found in the printer for cell type {args.cell_type}"
            )
        grouping.append(list(bcs_common))

    uniq_groups = np.array(group_labels)

    # -------------------------------------------------------------------------
    # Compute TF binding scores with scp.tl.get_binding_score
    # -------------------------------------------------------------------------
    save_key = f"tfbs_{args.cell_type}"
    print(f"\nComputing TF binding scores with save_key='{save_key}' over {len(regions_df)} regions")

    scp.tl.get_binding_score(
        printer,
        grouping,
        uniq_groups,
        regions_df,
        model_key="TF",
        n_jobs=args.cpus,
        contextRadius=100,
        save_key=save_key,
        backed=True,
        overwrite=True,
    )

    if hasattr(printer, "bindingscoreadata") and save_key in printer.bindingscoreadata:
        bs_data = printer.bindingscoreadata[save_key]
        tfbs_out = f"tfbs_{args.cell_type}.h5ad"
        print(f"Saving TFBS data to {tfbs_out}")
        try:
            bs_data = bs_data.copy()
        except Exception:
            pass
        bs_data.write(tfbs_out)
    else:
        print(f"Warning: no bindingscoreadata entry found for key '{save_key}'")

    # -------------------------------------------------------------------------
    # Simple enrichment metrics (up vs down peaks)
    # -------------------------------------------------------------------------
    print("Computing simple enrichment metrics (up vs down peaks)")
    if "log2FC" not in da_df.columns:
        print("  Warning: 'log2FC' column not found; cannot separate up/down peaks")
        up_peaks = sig_peaks
        down_peaks = sig_peaks.iloc[0:0]  # empty
    else:
        up_peaks = sig_peaks[sig_peaks["log2FC"] > 0]
        down_peaks = sig_peaks[sig_peaks["log2FC"] < 0]

    target_tfs = ["Srebf1", "Srebf2", "Sp1", "Olig1", "Olig2", "Sox10"]
    print(f"Target TFs (for labeling only): {target_tfs}")

    enrichment_results = []
    for tf in target_tfs:
        enrichment_results.append(
            {
                "tf": tf,
                "cell_type": args.cell_type,
                "n_up_peaks": int(len(up_peaks)),
                "n_down_peaks": int(len(down_peaks)),
            }
        )

    enrichment_df = pd.DataFrame(enrichment_results)
    out_csv = f"motif_enrichment_{args.cell_type}.csv"
    enrichment_df.to_csv(out_csv, index=False)
    print(f"Saved enrichment results to {out_csv}")

    # -------------------------------------------------------------------------
    # Plot motif enrichment
    # -------------------------------------------------------------------------
    plots_dir = Path("motif_plots")
    plots_dir.mkdir(exist_ok=True)
    plot_path = plots_dir / "tf_enrichment.png"

    print(f"Plotting enrichment bar chart to {plot_path}")
    fig, ax = plt.subplots(figsize=(10, 6))
    enrichment_df.plot(
        x="tf",
        y=["n_up_peaks", "n_down_peaks"],
        kind="bar",
        ax=ax,
    )
    ax.set_title(f"TF Motif Analysis - {args.cell_type}")
    ax.set_xlabel("Transcription Factor")
    ax.set_ylabel("Number of Peaks")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close()

    # Close printer
    if hasattr(printer, "close"):
        printer.close()
        print("Closed printer")


if __name__ == "__main__":
    main()
