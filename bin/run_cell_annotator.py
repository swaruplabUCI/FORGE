#!/usr/bin/env python3
import argparse
import scanpy as sc
import pandas as pd
import logging
from cell_annotator import CellAnnotator  # Assuming this is the main class or function
from h5ad_compat import sanitize_adata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--species', default='human')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info(f"Loading data from {args.input}")
    adata = sc.read_h5ad(args.input)

    logger.info("Running CellAnnotator...")
    # NOTE: Adjust this call based on the actual API of the cell-annotator package
    # Example assumption:
    annotator = CellAnnotator(species=args.species)
    adata = annotator.annotate(adata) 

    logger.info(f"Saving annotated data to {args.output}")
    sanitize_adata(adata, args.output)
    adata.write_h5ad(args.output)
    
    with open("cell_type_annotation_report.txt", "w") as f:
        f.write("Cell Annotator Run Complete\\n")
        f.write(f"Species: {args.species}\\n")
        if 'cell_type' in adata.obs:
             f.write(str(adata.obs['cell_type'].value_counts()))

if __name__ == "__main__":
    main()
