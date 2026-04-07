#!/usr/bin/env python

import argparse
import muon as mu
from h5ad_compat import sanitize_adata

def main():
    p = argparse.ArgumentParser(description="Export RNA modality from MuData for DORC")
    p.add_argument("--mudata", required=True, help="Input MuData file (.h5mu)")
    p.add_argument("--out", required=True, help="Output RNA h5ad for DORC")
    args = p.parse_args()

    mdata = mu.read_h5mu(args.mudata)
    rna = mdata.mod["rna"]
    # Propagate cell type annotations from global obs to RNA modality
    for col in ['cell_type', 'scanvi_prediction', 'celltypist_prediction',
                'cell_type_prediction']:
        if col in mdata.obs.columns and col not in rna.obs.columns:
            rna.obs[col] = mdata.obs.loc[rna.obs.index, col].values
    # FIX-77: Always overwrite sample_id from global obs (RNA modality may have
    # library-level names from concatenation, while global obs has donor-level IDs)
    if 'sample_id' in mdata.obs.columns:
        rna.obs['sample_id'] = mdata.obs.loc[rna.obs.index, 'sample_id'].values
    sanitize_adata(rna, args.out)
    rna.write_h5ad(args.out)

if __name__ == "__main__":
    main()
