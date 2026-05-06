#!/usr/bin/env python3
"""Per-(cell_type, gene) promoter MSFP compute, stratified by condition.

One task per cell type; loads printer + peak matrix once and iterates every
gene in --genes. Writes one standalone h5ad per gene shaped for
render_promoter_msfp_overlay.py:

    obs.index           = condition labels
    obs["name"]         = condition labels
    uns["scales"]       = np.arange(2, 101)
    obsm[<chr:start-end>] = (n_conditions, n_modes, n_positions)

Cells are subset to (peak_matrix.obs[cell_type_col] == cell_type) AND
(peak_matrix.obs[condition_col] in {control, treatment}); other conditions
are dropped. Empty subsets exit cleanly with a warning (no h5ad emitted).
"""

import os
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scprinter as scp

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--printer", required=True)
    p.add_argument("--peak-matrix", required=True)
    p.add_argument("--coords", required=True, help="gene_coordinates.json")
    p.add_argument("--genes", required=True, help="comma-separated gene list")
    p.add_argument("--cell-type", required=True)
    p.add_argument("--cell-type-col", required=True)
    p.add_argument("--condition-col", required=True)
    p.add_argument("--control-condition", required=True)
    p.add_argument("--treatment-condition", required=True)
    p.add_argument("--genome", default="mm10")
    p.add_argument("--cpus", type=int, default=8)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def tail_only(barcode):
    barcode = str(barcode)
    return barcode.split(":", 1)[1] if ":" in barcode else barcode


def safe_label(s):
    return str(s).replace("/", "_").replace(" ", "_").replace("(", "_").replace(")", "_")


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome_obj = getattr(scp.genome, args.genome)
    printer = scp.load_printer(args.printer, genome_obj)
    printer.load_disp_model()

    pm = ad.read_h5ad(args.peak_matrix)
    for col in (args.cell_type_col, args.condition_col):
        if col not in pm.obs.columns:
            raise ValueError(
                f"Missing '{col}' in peak_matrix.obs. Found: {pm.obs.columns.tolist()}"
            )

    keep_conds = {args.control_condition, args.treatment_condition}
    mask = (pm.obs[args.cell_type_col].astype(str) == args.cell_type) & \
           (pm.obs[args.condition_col].astype(str).isin(keep_conds))
    sub = pm[mask].copy()
    if sub.n_obs == 0:
        print(f"[WARN] No cells match cell_type={args.cell_type!r} "
              f"AND {args.condition_col} in {keep_conds}; emitting nothing.")
        return

    df = pd.DataFrame({
        "barcode": [tail_only(b) for b in sub.obs_names.astype(str)],
        "group":   sub.obs[args.condition_col].astype(str).values,
    })
    printer_bcs = set(map(str, printer.obs_names))
    df = df[df["barcode"].isin(printer_bcs)].copy()
    if df.empty:
        print(f"[WARN] After printer-barcode intersect, 0 cells remain for "
              f"ct={args.cell_type}; emitting nothing.")
        return

    barcodeGroups = df[["barcode", "group"]].copy()
    grouping, uniq_groups = scp.utils.df2cell_grouping(printer, barcodeGroups)
    order = np.argsort(uniq_groups)
    grouping = [grouping[i] for i in order]
    uniq_groups = uniq_groups[order]
    print(f"[INFO] ct={args.cell_type} groups={list(uniq_groups)} n_cells={len(df)}",
          flush=True)

    if len(uniq_groups) < 2:
        print(f"[WARN] Only one condition present for ct={args.cell_type}: "
              f"{list(uniq_groups)}. Emitting single-row panels (no diff).",
              flush=True)

    with open(args.coords, "r") as f:
        coords = json.load(f)

    safe_ct = safe_label(args.cell_type)
    modes = np.arange(2, 101)
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]

    for gene in genes:
        if gene not in coords or not coords[gene].get("found"):
            print(f"[WARN] Gene not found in coords: {gene}; skipping", flush=True)
            continue
        promoter = coords[gene]["promoter"]
        chrom, start, end = promoter["chr"], int(promoter["start"]), int(promoter["end"])
        regions_df = pd.DataFrame([{"Chromosome": chrom, "Start": start, "End": end}])

        save_key = f"promoter_msfp_{gene}_{safe_ct}"
        try:
            scp.tl.get_footprint_score(
                printer, grouping, uniq_groups, regions_df,
                modes=modes, n_jobs=args.cpus,
                save_key=save_key, backed=False, overwrite=True,
            )
        except Exception as e:
            print(f"[WARN] {gene}: get_footprint_score failed: {e}; skipping",
                  flush=True)
            continue

        try:
            fp_data = printer.footprintsadata[save_key]
        except Exception as e:
            print(f"[WARN] {gene}: could not read save_key={save_key}: {e}; skipping",
                  flush=True)
            continue

        obsm_keys = list(fp_data.obsm_keys())
        if not obsm_keys:
            print(f"[WARN] {gene}: empty obsm in fp_data; skipping", flush=True)
            continue

        arr = np.array(fp_data.obsm[obsm_keys[0]])
        if arr.ndim == 2:
            n_modes = len(modes)
            n_pos = arr.shape[1] // n_modes
            if n_pos == 0:
                print(f"[WARN] {gene}: degenerate obsm shape {arr.shape}; skipping",
                      flush=True)
                continue
            arr = arr.reshape(arr.shape[0], n_modes, n_pos)

        cond_labels = list(fp_data.obs_names.astype(str))
        locus_str = f"{chrom}:{start}-{end}"
        out_h5 = ad.AnnData(
            X=np.zeros((arr.shape[0], 1), dtype=np.float32),
            obs=pd.DataFrame({"name": cond_labels}, index=cond_labels),
        )
        out_h5.obsm[locus_str] = arr
        out_h5.uns["scales"] = modes
        out_h5.uns["gene"] = gene
        out_h5.uns["cell_type"] = args.cell_type

        out_path = outdir / f"promoter_fp_{gene}__{safe_ct}.h5ad"
        out_h5.write_h5ad(out_path)
        print(f"[OK] wrote {out_path} arr={arr.shape}", flush=True)

    if hasattr(printer, "close"):
        printer.close()


if __name__ == "__main__":
    main()
