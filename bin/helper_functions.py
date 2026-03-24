#!/usr/bin/env python
"""Helper functions for RNA-seq processing pipeline"""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import tables
from scipy.stats import median_abs_deviation
from typing import Dict

def is_outlier(adata, metric: str, nmads: int):
    """Check for outliers using MAD method"""
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""
    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d

def anndata_from_h5(file: str, analyzed_barcodes_only: bool = True) -> ad.AnnData:
    """Load an output h5 file into an AnnData object."""
    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            print('Assuming we are loading a "filtered" file that contains only cells.')
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]

    adata = ad.AnnData(X=X,
                      obs={'barcode': d.pop('barcodes').astype(str)},
                      var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                         else d.pop('name')).astype(str)},
                      dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    if 'genes' in d.keys():
        d['id'] = d.pop('genes')
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    _fill_adata_slots_automatically(adata, d)
    
    return adata
