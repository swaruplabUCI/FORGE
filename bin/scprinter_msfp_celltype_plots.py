#!/usr/bin/env python3
import os
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import scprinter as scp


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--printer", required=True)
    p.add_argument("--peak-matrix", required=True)
    p.add_argument("--coords", required=True, help="gene_coordinates.json from RESOLVE_GENE_COORDINATES")
    p.add_argument("--genes", default="Srebf1,Hspa8,Hmgb1")
    p.add_argument("--cell-type-key", default="cell_type")
    p.add_argument("--outdir", required=True)
    p.add_argument("--cpus", type=int, default=8)
    p.add_argument("--scale", type=int, default=50)
    return p.parse_args()

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def tail_only(barcode: str) -> str:
    barcode = str(barcode)
    return barcode.split(":", 1)[1] if ":" in barcode else barcode


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    printer = scp.load_printer(args.printer, scp.genome.mm10)
    printer.load_disp_model()

    pm = ad.read_h5ad(args.peak_matrix)

    ct_key = args.cell_type_key
    if ct_key not in pm.obs.columns:
        raise ValueError(f"Missing {ct_key} in peak_matrix.obs. Found: {pm.obs.columns.tolist()}")

    # Build barcode->cell_type table, normalize to match printer barcodes
    df = pd.DataFrame({
        "barcode": [tail_only(b) for b in pm.obs_names.astype(str)],
        "cell_type": pm.obs[ct_key].astype(str).values
    })

    # Keep only barcodes present in printer
    printer_bcs = set(map(str, printer.obs_names))
    df = df[df["barcode"].isin(printer_bcs)].copy()

    # scprinter expects first col barcodes, second col group labels
    barcodeGroups = df[["barcode", "cell_type"]].copy()

    # Convert to grouping
    grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcodeGroups)

    # Sort groups by name for stable plots
    order = np.argsort(uniq_groups)
    grouping = [grouping[i] for i in order]
    uniq_groups = uniq_groups[order]

    # Load coordinates
    with open(args.coords, "r") as f:
        coords = json.load(f)

    genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    for gene in genes:
        if gene not in coords or not coords[gene].get("found"):
            print(f"[WARN] Gene not found in coords: {gene}; skipping")
            continue
        promoter = coords[gene]["promoter"]
        region_str = f"{promoter['chr']}:{promoter['start']}-{promoter['end']}"

        regions_df = pd.DataFrame([{
            "Chromosome": promoter["chr"],
            "Start": int(promoter["start"]),
            "End": int(promoter["end"]),
        }])

        save_key = f"celltype_footprints_{gene}"

        # Compute footprints once per gene across all cell types
        scp.tl.get_footprint_score(
            printer,
            grouping,
            uniq_groups,
            regions_df,
            modes=np.arange(2, 101),
            n_jobs=args.cpus,
            save_key=save_key,
            backed=False,
            overwrite=True,
        )

        # Plot stacked footprints (one scale) across all cell types
        fig = plt.figure(figsize=(10, max(6, 0.25 * len(uniq_groups))))
        scp.pl.plot_footprints(
            printer,
            save_key=save_key,
            group_names=uniq_groups,
            row_label=uniq_groups,
            region=region_str,
            stack=True,
            scales=[args.scale],
            cmap="Blues",
            vmin=0.5,
            vmax=2.0,
        )
        plt.tight_layout()
        figpath = outdir / f"{gene}_footprints_all_celltypes_scale{args.scale}.png"
        plt.savefig(figpath, dpi=200)
        plt.close()
        print(f"[OK] wrote {figpath}")

        # Low-hanging fruit: insertion profile for first ~8 cell types
        n_show = min(8, len(uniq_groups))
        fig, axs = plt.subplots(n_show, 1, figsize=(10, 1.5 * n_show), sharex=True)
        if n_show == 1:
            axs = [axs]
        scp.pl.plot_group_atac(
            printer,
            grouping[:n_show],
            np.arange(n_show),
            region_str,
            ax=list(axs),
            smooth=5,
            color="red",
        )
        for i in range(n_show):
            axs[i].set_title(str(uniq_groups[i]))
        plt.tight_layout()
        figpath2 = outdir / f"{gene}_insertions_first{n_show}_celltypes.png"
        plt.savefig(figpath2, dpi=200)
        plt.close()
        print(f"[OK] wrote {figpath2}")

    if hasattr(printer, "close"):
        printer.close()


if __name__ == "__main__":
    main()
