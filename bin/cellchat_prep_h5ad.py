#!/usr/bin/env python3
import argparse
import anndata as ad
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene_matrix", required=True, help="SnapATAC2 gene activity h5ad")
    ap.add_argument("--meta_h5ad", required=True, help="Peak matrix (annotated) h5ad containing obs[cell_type,condition]")
    ap.add_argument("--cell_type_key", default="cell_type")
    ap.add_argument("--condition_key", default="condition")
    ap.add_argument("--control", required=True)
    ap.add_argument("--treatment", required=True)
    ap.add_argument("--out_global", required=True)
    ap.add_argument("--out_control", required=True)
    ap.add_argument("--out_treatment", required=True)
    args = ap.parse_args()

    g = ad.read_h5ad(args.gene_matrix)
    m = ad.read_h5ad(args.meta_h5ad, backed="r")

    # Align on shared barcodes
    common = pd.Index(g.obs_names).intersection(pd.Index(m.obs_names))
    if len(common) == 0:
        raise RuntimeError("No shared obs_names between gene_matrix and meta_h5ad.")

    g = g[common, :].copy()
    meta = m.obs.loc[common].copy()

    # Validate keys
    for k in [args.cell_type_key, args.condition_key]:
        if k not in meta.columns:
            raise KeyError(f"Missing '{k}' in meta_h5ad.obs. Found columns: {list(meta.columns)}")

    # Transfer minimal metadata needed by CellChat (plus sample if present)
    keep_cols = [args.cell_type_key, args.condition_key]
    if "sample" in meta.columns:
        keep_cols.append("sample")
    if "batch" in meta.columns:
        keep_cols.append("batch")

    for c in keep_cols:
        g.obs[c] = meta[c].values

    # Write global
    g.write_h5ad(args.out_global, compression="gzip")

    # Write per-condition
    cond = g.obs[args.condition_key].astype(str)
    g_ctrl = g[cond == str(args.control), :].copy()
    g_trt  = g[cond == str(args.treatment), :].copy()

    if g_ctrl.n_obs == 0 or g_trt.n_obs == 0:
        raise RuntimeError(
            f"Condition split produced empty set(s): "
            f"control={g_ctrl.n_obs}, treatment={g_trt.n_obs}"
        )

    g_ctrl.write_h5ad(args.out_control, compression="gzip")
    g_trt.write_h5ad(args.out_treatment, compression="gzip")

    try:
        m.file.close()
    except Exception:
        pass

if __name__ == "__main__":
    main()
