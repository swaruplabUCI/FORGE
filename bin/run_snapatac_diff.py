#!/usr/bin/env python3
import snapatac2 as snap
import pandas as pd
import anndata as ad
import numpy as np
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--peak_matrix', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--cell_type', required=True)
    parser.add_argument('--control_condition', required=True)  
    parser.add_argument('--treatment_condition', required=True)  
    parser.add_argument('--min_cells', type=int, default=50)
    parser.add_argument('--fdr', type=float, default=0.05)
    parser.add_argument('--output_prefix', required=True)
    parser.add_argument('--cell_type_key', default='celltypist_prediction')
    parser.add_argument('--condition_key', default='condition_group')
    args = parser.parse_args()

    # Sanitize output_prefix: cell type names may contain / or spaces
    args.output_prefix = args.output_prefix.replace('/', '_').replace(' ', '_')

    control_condition = args.control_condition
    treatment_condition = args.treatment_condition
    cond_key = args.condition_key
    
    print(f"=" * 60)
    print(f"Running differential accessibility analysis")
    print(f"Cell type: {args.cell_type}")
    print(f"Comparison: {treatment_condition} vs {control_condition}")
    print(f"=" * 60)
    
    # Load data
    print(f"\nLoading peak matrix: {args.peak_matrix}")
    adata = ad.read_h5ad(args.peak_matrix)
    print(f"Loaded {adata.n_obs} cells × {adata.n_vars} peaks")
    
    #  Validate required columns
    ct_key = args.cell_type_key
    if ct_key not in adata.obs.columns:
        print(f"ERROR: '{ct_key}' column not found")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    if cond_key not in adata.obs.columns:
        print(f"ERROR: '{cond_key}' column not found")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)
    
    # Dynamic cell type filtering 
    cell_mask = adata.obs[ct_key] == args.cell_type
    adata_ct = adata[cell_mask].copy()
    
    print(f"\nCells in {args.cell_type}: {adata_ct.n_obs}")
    
    if adata_ct.n_obs == 0:
        print(f"ERROR: No cells found for cell type '{args.cell_type}'")
        print(f"Available cell types: {adata.obs[ct_key].unique()}")
        sys.exit(1)
    
    # Get cell indices for each condition
    control_mask = adata_ct.obs[cond_key] == control_condition
    treatment_mask = adata_ct.obs[cond_key] == treatment_condition
    
    ctrl_idx = np.where(control_mask)[0]
    treat_idx = np.where(treatment_mask)[0]
    
    print(f"\nControl cells ({control_condition}): {len(ctrl_idx)}")
    print(f"Treatment cells ({treatment_condition}): {len(treat_idx)}")
    
    # Check minimum cells
    if len(ctrl_idx) < args.min_cells or len(treat_idx) < args.min_cells:
        print(f"\nWARNING: Insufficient cells for reliable analysis")
        print(f"Minimum required: {args.min_cells} per group")
        # Still proceed but flag it
    
    # Create output directory
    plot_dir = f"plots_{args.output_prefix}"
    Path(plot_dir).mkdir(exist_ok=True)
    
    # Run differential test
    print(f"\nRunning SnapATAC2 differential test...")
    result = snap.tl.diff_test(
        data=adata_ct,
        cell_group1=treat_idx,  # Treatment group
        cell_group2=ctrl_idx    # Control group (reference)
    )
    
    # Convert to pandas and save
    result_df = result.to_pandas()
    result_df.columns = ['peak', 'log2FC', 'pval', 'padj']
    
    output_file = f"DA_peaks_{args.output_prefix}.csv"
    result_df.to_csv(output_file, index=False)
    print(f"\nSaved results to: {output_file}")

    # Generate visualizations
    print(f"\nGenerating visualizations...")
    try:
        # Create volcano plot using SnapATAC2's plotting function
        fig = snap.pl.diff_test(adata_ct, result, show=False)
        
        # Save as HTML (interactive)
        html_path = f"{plot_dir}/volcano_{args.output_prefix}.html"
        fig.write_html(html_path)
        print(f"  - HTML saved: {html_path}")
        
        # Save as PNG (static)
        png_path = f"{plot_dir}/volcano_{args.output_prefix}.png"
        try:
            fig.write_image(png_path, scale=2)
            print(f"  - PNG saved: {png_path}")
        except Exception as e:
            print(f"  ! PNG generation failed (kaleido might be missing): {e}")

        # Save as PDF (vector)
        pdf_path = f"{plot_dir}/volcano_{args.output_prefix}.pdf"
        try:
            fig.write_image(pdf_path)
            print(f"  - PDF saved: {pdf_path}")
        except Exception as e:
            print(f"  ! PDF generation failed: {e}")
            
    except Exception as e:
        print(f"WARNING: Visualization generation failed: {e}")
    
    # Summary statistics
    n_sig = (result_df['padj'] < args.fdr).sum()
    print(f"\n{'=' * 60}")
    print(f"RESULTS SUMMARY")
    print(f"{'=' * 60}")
    print(f"Total peaks tested: {len(result_df)}")
    print(f"Significant peaks (FDR < {args.fdr}): {n_sig}")
    
    if n_sig > 0:
        n_up = ((result_df['padj'] < args.fdr) & (result_df['log2FC'] > 0)).sum()
        n_down = ((result_df['padj'] < args.fdr) & (result_df['log2FC'] < 0)).sum()
        print(f"  ↑ Up in {treatment_condition}: {n_up}")
        print(f"  ↓ Down in {treatment_condition}: {n_down}")
    
    print(f"{'=' * 60}\n")

if __name__ == "__main__":
    main()
