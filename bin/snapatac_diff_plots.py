#!/usr/bin/env python3
import snapatac2 as snap
import anndata as ad
import numpy as np
import argparse
import sys
import pandas as pd
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--peak_matrix', required=True, help="The h5ad file from merge_annotations.py")
    p.add_argument('--cell_type', required=True, help="e.g., 'Excitatory_neurons'")
    p.add_argument('--control', dest='control_condition', required=True)
    p.add_argument('--treatment', dest='treatment_condition', required=True)
    p.add_argument('--condition_key', default='condition')
    p.add_argument('--cell_type_key', default='cell_type')
    p.add_argument('--genome', default='mm10', choices=['mm10', 'hg38'])
    p.add_argument('--fdr', type=float, default=0.05)
    p.add_argument('--output_prefix', required=True)
    args = p.parse_args()

    # Create output directory
    plot_dir = Path(f"plots_{args.output_prefix}")
    plot_dir.mkdir(exist_ok=True, parents=True)

    print(f"Loading annotated matrix: {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)

    # Subset to the specific cell type
    adata_ct = adata[adata.obs[args.cell_type_key] == args.cell_type].copy()
    if adata_ct.n_obs == 0:
        print(f"ERROR: Cell type '{args.cell_type}' not found.")
        sys.exit(1)

    # 1. PRE-TEST PLOTTING: UMAP (From tutorial section 1)
    # This confirms your treatment/control groups are distributed correctly
    if 'X_umap' in adata_ct.obsm:
        print("Generating UMAP plot...")
        snap.pl.umap(adata_ct, color=args.condition_key, 
                    out_file=str(plot_dir / "umap_by_condition.html"), show=False)

    # 2. PRE-TEST PLOTTING: Marker Regions (The 'quick-and-dirty' z-score method)
    # This identifies regions specifically enriched in treatment vs control
    print("Identifying marker regions for heatmap...")
    marker_peaks = snap.tl.marker_regions(adata_ct, groupby=args.condition_key)
    snap.pl.regions(adata_ct, groupby=args.condition_key, peaks=marker_peaks, 
                   out_file=str(plot_dir / "marker_heatmap.html"), show=False)

    # 3. REGRESSION-BASED TEST (The section you wanted to include/reach)
    print(f"Running differential test: {args.treatment_condition} vs {args.control_condition}")
    ctrl_mask = adata_ct.obs[args.condition_key] == args.control_condition
    treat_mask = adata_ct.obs[args.condition_key] == args.treatment_condition
    
    result = snap.tl.diff_test(
        data=adata_ct,
        cell_group1=treat_mask,
        cell_group2=ctrl_mask
    )
    
    df = result.to_pandas()
    df.columns = ['peak', 'log2FC', 'pval', 'padj']
    df.to_csv(f"DA_peaks_{args.output_prefix}.csv", index=False)

    # 4. POST-TEST PLOTTING: Motif Enrichment
    # This visualizes TFs enriched in the significantly open regions
    sig_peaks = df[df['padj'] < args.fdr]['peak'].unique().tolist()
    
    if len(sig_peaks) > 50: # Enrichment requires a reasonable number of peaks
        print(f"Running motif enrichment on {len(sig_peaks)} significant peaks...")
        # Select genome fasta based on argument
        genome_fasta = snap.genome.mm10 if args.genome == 'mm10' else snap.genome.hg38
        
        motifs = snap.tl.motif_enrichment(
            motifs=snap.datasets.cis_bp(unique=True),
            regions=sig_peaks,
            genome_fasta=genome_fasta
        )
        snap.pl.motif_enrichment(motifs, out_file=str(plot_dir / "motif_enrichment.html"), show=False)
    else:
        print("Insufficient significant peaks for motif enrichment plotting.")

    print(f"Analysis complete. Results in DA_peaks_{args.output_prefix}.csv and {plot_dir}/")

if __name__ == "__main__":
    main()
