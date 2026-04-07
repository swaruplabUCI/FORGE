"""
h5ad_compat.py — Cross-container anndata compatibility sanitization.

Newer anndata (>=0.10, in scgpu_extended) writes certain uns keys with
encoding types (e.g. 'null 0.1.0') that older anndata (in snapatac_extended,
scenicplus) cannot read. This module strips those keys in-place before
writing so that all downstream containers can read the output.

Stripped keys are saved to a JSON sidecar file (<output>.uns_stripped.json)
so they can be re-attached later if needed.

Usage:
    from h5ad_compat import sanitize_adata
    sanitize_adata(adata, "output.h5ad")
    adata.write_h5ad("output.h5ad")

For MuData:
    from h5ad_compat import sanitize_mudata
    sanitize_mudata(mdata, "output.h5mu")
    mdata.write_h5mu("output.h5mu")

To restore stripped keys:
    from h5ad_compat import restore_uns
    restore_uns(adata, "output.h5ad")
"""

import json
from pathlib import Path


def _extract_from_uns(uns, label=""):
    """Extract and remove problematic keys from an uns dict.
    Returns a dict of what was stripped (empty if nothing)."""
    stripped = {}

    # uns/log1p/base — written as null encoding by newer scanpy/anndata
    if "log1p" in uns and isinstance(uns["log1p"], dict):
        if "base" in uns["log1p"]:
            val = uns["log1p"].pop("base")
            stripped["log1p.base"] = _serialize(val)
            if not uns["log1p"]:
                del uns["log1p"]

    # uns/rank_genes_groups/params/layer — null encoding
    if "rank_genes_groups" in uns:
        rgg = uns["rank_genes_groups"]
        if isinstance(rgg, dict) and "params" in rgg and isinstance(rgg["params"], dict):
            if "layer" in rgg["params"]:
                val = rgg["params"].pop("layer")
                stripped["rank_genes_groups.params.layer"] = _serialize(val)

    return stripped


def sanitize_adata(adata, output_path=None):
    """Remove uns keys incompatible with older anndata, saving them to a sidecar.

    Args:
        adata: AnnData object to sanitize in-place.
        output_path: If provided, stripped keys are saved to
                     <output_path>.uns_stripped.json. If None, keys are
                     still stripped but no sidecar is written.
    """
    stripped = _extract_from_uns(adata.uns)
    if stripped and output_path is not None:
        _write_sidecar(stripped, output_path)


def sanitize_mudata(mdata, output_path=None):
    """Sanitize a MuData object and all its modalities."""
    all_stripped = {}

    # Top-level uns
    top = _extract_from_uns(mdata.uns)
    if top:
        all_stripped["global"] = top

    # Each modality
    if hasattr(mdata, "mod"):
        for mod_name in mdata.mod:
            mod_stripped = _extract_from_uns(mdata.mod[mod_name].uns)
            if mod_stripped:
                all_stripped[mod_name] = mod_stripped

    if all_stripped and output_path is not None:
        _write_sidecar(all_stripped, output_path)


def restore_uns(adata, output_path):
    """Re-attach stripped keys from a sidecar file into adata.uns.

    Args:
        adata: AnnData object to restore into.
        output_path: Path to the h5ad file whose sidecar to read.
    """
    sidecar = Path(str(output_path) + ".uns_stripped.json")
    if not sidecar.exists():
        return

    with open(sidecar) as f:
        stripped = json.load(f)

    for dotkey, val in stripped.items():
        _set_nested(adata.uns, dotkey, val)


def _write_sidecar(stripped, output_path):
    """Write stripped keys to a JSON sidecar next to the output file."""
    sidecar = Path(str(output_path) + ".uns_stripped.json")
    with open(sidecar, "w") as f:
        json.dump(stripped, f, indent=2)


def _serialize(val):
    """Make a value JSON-serializable."""
    if val is None:
        return None
    try:
        json.dumps(val)
        return val
    except (TypeError, ValueError):
        return str(val)


def _set_nested(d, dotkey, val):
    """Set a nested dict value from a dot-separated key path.
    e.g. _set_nested(uns, 'log1p.base', None) → uns['log1p']['base'] = None
    """
    parts = dotkey.split(".")
    for part in parts[:-1]:
        if part not in d:
            d[part] = {}
        d = d[part]
    d[parts[-1]] = val
