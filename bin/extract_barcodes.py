#!/usr/bin/env python3
import argparse
from pathlib import Path
import anndata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5ad_files', nargs='+', required=True)
    parser.add_argument('--output_dir', default='.')
    args = parser.parse_args()
    
    for h5ad_path in args.h5ad_files:
        # For ATAC data from snapATAC2
        import anndata
        adata = anndata.read_h5ad(h5ad_path)
        sample_name = Path(h5ad_path).stem.replace('_filtered', '')
        
        # Extract barcodes
        barcodes = adata.obs_names.tolist()
        
        # Write barcode file
        output_file = f"{args.output_dir}/{sample_name}_barcodes.txt"
        with open(output_file, 'w') as f:
            for bc in barcodes:
                f.write(f"{bc}\n")
        
        print(f"Wrote {len(barcodes)} QC-filtered barcodes to {output_file}")
        
        # DIAGNOSTIC: Print first few barcodes to verify format
        print(f"   Sample barcodes (first 3): {barcodes[:3]}")

if __name__ == "__main__":
    main()
