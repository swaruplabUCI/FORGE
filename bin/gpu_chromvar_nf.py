#!/usr/bin/env python

import argparse
import time
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import os

import cupy as cp
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

# GPU allocator setup
rmm.reinitialize(managed_memory=True, pool_allocator=True)
cp.cuda.set_allocator(rmm_cupy_allocator)

def init_scprinter(cache_dir):
    import scprinter as scp
    ds = scp.datasets.datasets()
    ds.path = cache_dir
    return scp

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peak-matrix", required=True,
                   help="Cell x peak matrix h5ad (staged into work dir by Nextflow)")
    p.add_argument('--metadata', required=False, default=None, 
                   help='Path to metadata CSV (optional if metadata in peak_matrix.obs)')
    p.add_argument("--metadata-id-col", default="cell_id",
                   help="Column in metadata matching AnnData obs_names [default: cell_id]")
    p.add_argument("--cache-dir", required=True,
                   help="scPrinter cache dir (e.g. /dfs7/.../ref/scprinter)")
    p.add_argument("--pfms", required=True,
                   help="Path to JASPAR PFMs (e.g. JASPAR2022_core_nonredundant.jaspar)")
    p.add_argument("--genome", required=True,
                   help="Genome key for scprinter (e.g. mm10)")
    p.add_argument("--out-prefix", required=True,
                   help="Output prefix (will write <prefix>.h5ad and <prefix>_deviations.h5ad)")

    # Optional hooks
    p.add_argument("--chunk-size", type=int, default=10000,
                   help="Cells per chunk for GPU deviation computation [default: 10000]")
    p.add_argument("--da-peaks", default="", help="Optional SnapATAC DA peaks TSV/CSV")
    p.add_argument("--cicero-conns", default="", help="Optional Cicero connections tsv.gz")
    p.add_argument("--ccan-assignments", default="", help="Optional CCAN assignments tsv.gz")

    return p.parse_args()

def main():
    args = parse_args()
    scp = init_scprinter(args.cache_dir)

    # Load peak matrix
    print(f"Loading peak matrix from {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)
    print(f"  adata: n_cells={adata.n_obs}, n_peaks={adata.n_vars}")
    print(f"  adata.obs columns: {list(adata.obs.columns)}")

    if not sp.issparse(adata.X):
        from scipy.sparse import csr_matrix
        adata.X = csr_matrix(adata.X)
    adata.X = adata.X.astype(np.float32)

    # Optionally merge external metadata
    if args.metadata and os.path.exists(args.metadata):
        print(f"Loading external metadata from {args.metadata}")
        meta = pd.read_csv(args.metadata)
        
        if args.metadata_id_col in meta.columns:
            # Cell-level metadata: merge by cell ID
            print(f"  Merging by cell-level ID column: {args.metadata_id_col}")
            meta = meta.set_index(args.metadata_id_col)
            common_ids = adata.obs_names.intersection(meta.index)
            print(f"  Overlap: {len(common_ids)} cells")
            
            if len(common_ids) > 0:
                adata = adata[adata.obs_names.isin(common_ids)].copy()
                meta_aligned = meta.reindex(adata.obs_names)
                for col in meta_aligned.columns:
                    if col not in adata.obs.columns:
                        adata.obs[col] = meta_aligned[col].values
                print(f"  Added metadata columns: {list(meta_aligned.columns)}")
        
        elif 'sample_id' in meta.columns:
            # Sample-level metadata: merge by sample
            print(f"  Merging sample-level metadata (no '{args.metadata_id_col}' column)")
            
            # Extract sample from cell barcode (format: sample1:BARCODE-1)
            if 'sample' not in adata.obs.columns:
                adata.obs['sample'] = adata.obs_names.str.split(':').str[0]
            
            meta_dict = meta.set_index('sample_id').to_dict('index')
            for col in meta.columns:
                if col != 'sample_id' and col not in adata.obs.columns:
                    adata.obs[col] = adata.obs['sample'].map(
                        {k: v[col] for k, v in meta_dict.items()}
                    )
            print(f"  Merged sample-level metadata for {adata.obs['sample'].nunique()} samples")
        else:
            print(f"  WARNING: Could not merge metadata (no '{args.metadata_id_col}' or 'sample_id' column)")
    else:
        print("  No external metadata provided; using peak matrix .obs")
        if args.metadata:
            print(f"  WARNING: Metadata file '{args.metadata}' does not exist")

    # Verify required columns for chromVAR grouping
    print("\nVerifying metadata columns:")
    if 'cell_type' in adata.obs.columns:
        print(f"  ✓ 'cell_type' found: {adata.obs['cell_type'].nunique()} types")
    else:
        print("  ⚠ 'cell_type' not found in .obs; chromVAR grouping may be limited")
    
    if 'sample' in adata.obs.columns:
        print(f"  ✓ 'sample' found: {adata.obs['sample'].nunique()} samples")

    # ChromVAR background
    genome_obj = getattr(scp.genome, args.genome)
    fa = genome_obj.fetch_fa()

    print("\nSampling chromVAR background peaks...")
    t0 = time.time()
    scp.chromvar.sample_bg_peaks(
        adata,
        genome=genome_obj,
        method="chromvar",
        niterations=250
    )
    print(f"Bias+BG took: {time.time() - t0:.2f} sec")

    # Motif matching
    print(f"\nLoading PFMs from {args.pfms}")
    t0 = time.time()
    motifs = scp.motifs.Motifs(
        args.pfms,
        fa,
        bg=list(adata.uns["bg_freq"]),
        n_jobs=32,
        pvalue=5e-5,
        mode="motifmatchr",
    )
    motifs.prep_scanner(pvalue=5e-5)
    motifs.chromvar_scan(adata)
    print(f"Motif scan took: {time.time() - t0:.2f} sec")

    motif_out = args.out_prefix + ".h5ad"
    adata.write_h5ad(motif_out, compression="gzip")
    print(f"Wrote motif-scanned matrix: {motif_out}")

    # GPU chromVAR deviations
    print("\nComputing GPU chromVAR deviations...")
    t0 = time.time()
    print(f"  chunk_size={args.chunk_size}, n_cells={adata.n_obs}, estimated chunks={-(-adata.n_obs // args.chunk_size)}")
    dev = scp.chromvar.compute_deviations(
        adata,
        chunk_size=args.chunk_size
    )
    print(f"GPU deviations took: {time.time() - t0:.2f} sec")

    dev_out = args.out_prefix + "_deviations.h5ad"
    if isinstance(dev, ad.AnnData):
        dev.write_h5ad(dev_out, compression="gzip")
        print(f"Wrote chromVAR deviations: {dev_out}")
    else:
        print(f"WARNING: compute_deviations returned {type(dev)}, not AnnData")

    print("\nDone.")

if __name__ == "__main__":
    main()
