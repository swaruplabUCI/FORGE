#!/usr/bin/env python3
"""
Concatenate multiple h5mu/cellismo files
"""

import mudata as md
from pathlib import Path
import argparse
import sys
import json

def parse_args():
    parser = argparse.ArgumentParser(
        description='Concatenate h5mu or cellismo files'
    )
    
    parser.add_argument('--input_files', nargs='+', required=True,
                       help='List of h5mu/cellismo files to concatenate')
    parser.add_argument('--output_file', type=str, required=True,
                       help='Output filename')
    parser.add_argument('--output_dir', type=str, default='.',
                       help='Output directory')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("H5MU CONCATENATION")
    print("="*80)
    
    # Load files
    print(f"\nLoading {len(args.input_files)} files...")
    mdatas = {}
    
    for idx, filepath in enumerate(args.input_files, 1):
        filepath = Path(filepath)
        
        if not filepath.exists():
            print(f"⚠ WARNING: {filepath} not found, skipping")
            continue
        
        print(f"  [{idx}] Loading {filepath.name}...")
        mdata = md.read_h5mu(filepath)
        
        # Add metadata
        mdata.obs['batch'] = filepath.stem
        mdata.obs['file_source'] = filepath.name
        
        # Make cell barcodes unique
        new_obs_names = [f"{filepath.stem}_{bc}" for bc in mdata.obs_names]
        mdata.obs_names = new_obs_names
        
        for mod_name in mdata.mod.keys():
            mdata.mod[mod_name].obs_names = new_obs_names
        
        mdatas[filepath.stem] = mdata
        print(f"      Shape: {mdata.shape}")
        print(f"      Modalities: {list(mdata.mod.keys())}")
    
    print(f"\n✓ Loaded {len(mdatas)} files")
    
    # Concatenate
    print("\nConcatenating...")
    mdata_combined = md.concat(
        mdatas,
        join='outer',
        merge='same',
        uns_merge='same',
        label='file_source',
        index_unique=None
    )
    
    print(f"\n✓ Combined MuData:")
    print(f"  Cells: {mdata_combined.n_obs:,}")
    print(f"  Features: {mdata_combined.n_var:,}")
    print(f"  Modalities: {list(mdata_combined.mod.keys())}")
    
    # Save
    output_file = OUTPUT_DIR / args.output_file
    mdata_combined.write_h5mu(output_file)
    print(f"\n✓ Saved: {output_file}")
    
    # Save report
    report = {
        'n_input_files': len(mdatas),
        'total_cells': int(mdata_combined.n_obs),
        'total_features': int(mdata_combined.n_var),
        'modalities': list(mdata_combined.mod.keys()),
        'files_processed': list(mdatas.keys())
    }
    
    report_file = OUTPUT_DIR / "integration_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✓ Saved report: {report_file}")
    
    print("\n" + "="*80)
    print("✅ CONCATENATION COMPLETE!")
    print("="*80)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
