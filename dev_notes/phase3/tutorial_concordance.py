#!/usr/bin/env python3
"""
tutorial_concordance.py — RNA vs ATAC cell-type concordance for the T2 tutorial run.

Follows the production method used for the four published datasets
(oneOff/20260602/{pbmc,brain,kidney}_concordance.py): read both modalities' labels
off the MultiVI h5mu, map each tool's vocabulary onto a shared broad vocabulary,
then report L1 raw agreement plus a confusion matrix.

Self-contained on purpose — it embeds the CellTypist->broad map rather than
importing the oneOff module, so it stays runnable from a clone of this repo.

WHAT THIS MEASURES, AND WHAT IT DOES NOT
----------------------------------------
The two arms are annotated independently: RNA by CellTypist on the gene matrix,
ATAC by the marker path on the peak/gene-activity matrix. Neither sees the other's
labels. So this is a genuine join across modalities.

On the tutorial subset it is a WIRING check, not a biological result. ATAC is
restricted to chr21+chr22 (~2.6% of the genome), so the ATAC marker panel has
almost nothing to score against and its labels collapse. The number this prints
is therefore expected to be near the chance floor. That is the documented,
intended behaviour of the subset — see the CAVEAT block in
configs/datasets/tutorial_pbmc.config. Do not read it as pipeline failure, and do
not quote it as a quality metric.

Run inside snapatac_extended.sif:

    singularity exec --bind /dfs7,/tmp \
        singularity_cache/snapatac_extended.sif \
        python3 dev_notes/phase3/tutorial_concordance.py \
            --h5mu results_tutorial_remeasure/multiome/multivi/multivi_integrated.h5mu \
            --outdir dev_notes/phase3/concordance
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Optional

import h5py
import numpy as np
import pandas as pd

# ── shared broad vocabulary ───────────────────────────────────────────────────
# Same 10 classes as the published-dataset analysis, so the tutorial number is
# comparable to the ones in oneOff/20260602/pbmc_concordance_combined.csv.
BROAD_CLASSES = (
    'CD4+ T', 'CD8+ T', 'Regulatory T', 'Unconventional T',
    'Memory B', 'Naive B', 'NK', 'Monocytes',
    'Dendritic Cells', 'Progenitors/Other',
)

_SKIP = {'Unknown', 'unknown', 'Unclassified', 'nan', 'None', ''}

_CELLTYPIST_MAP: dict[str, Optional[str]] = {
    # CD4+ T
    'Tcm/Naive helper T cells': 'CD4+ T',
    'Tem/Effector helper T cells': 'CD4+ T',
    'Tem/Effector helper T cells PD1+': 'CD4+ T',
    'Follicular helper T cells': 'CD4+ T',
    'Type 1 helper T cells': 'CD4+ T',
    'Type 17 helper T cells': 'CD4+ T',
    'Memory CD4+ cytotoxic T cells': 'CD4+ T',
    'Cycling T cells': 'CD4+ T',
    'T(agonist)': 'CD4+ T',
    # CD8+ T
    'Tcm/Naive cytotoxic T cells': 'CD8+ T',
    'Tem/Trm cytotoxic T cells': 'CD8+ T',
    'Tem/Temra cytotoxic T cells': 'CD8+ T',
    'Trm cytotoxic T cells': 'CD8+ T',
    'CD8a/a': 'CD8+ T',
    # Regulatory T
    'Regulatory T cells': 'Regulatory T',
    'Treg(diff)': 'Regulatory T',
    # Unconventional T
    'MAIT cells': 'Unconventional T',
    'NKT cells': 'Unconventional T',
    'gamma-delta T cells': 'Unconventional T',
    'CRTAM+ gamma-delta T cells': 'Unconventional T',
    'Double-negative thymocytes': 'Unconventional T',
    # B lineage
    'Memory B cells': 'Memory B',
    'Age-associated B cells': 'Memory B',
    'Germinal center B cells': 'Memory B',
    'Proliferative germinal center B cells': 'Memory B',
    'Plasmablasts': 'Memory B',
    'Plasma cells': 'Memory B',
    'Follicular B cells': 'Memory B',
    'Naive B cells': 'Naive B',
    'B cells': None,          # subtype unresolved -> excluded from denominator
    # NK
    'NK cells': 'NK',
    'CD16+ NK cells': 'NK',
    'CD16- NK cells': 'NK',
    'Cycling NK cells': 'NK',
    # Monocytes / macrophages
    'Classical monocytes': 'Monocytes',
    'Non-classical monocytes': 'Monocytes',
    'Monocytes': 'Monocytes',
    'Mono-mac': 'Monocytes',
    'Macrophages': 'Monocytes',
    'Alveolar macrophages': 'Monocytes',
    # DC
    'DC1': 'Dendritic Cells',
    'DC2': 'Dendritic Cells',
    'DC3': 'Dendritic Cells',
    'pDC': 'Dendritic Cells',
    # Progenitors / other
    'HSC/MPP': 'Progenitors/Other',
    'ILC3': 'Progenitors/Other',
    'MEMP': 'Progenitors/Other',
    'ELP': 'Progenitors/Other',
    'ETP': 'Progenitors/Other',
    'Megakaryocyte precursor': 'Progenitors/Other',
    'Megakaryocytes/platelets': 'Progenitors/Other',
    'Mast cells': 'Progenitors/Other',
    'Fibroblasts': 'Progenitors/Other',
}

# ATAC marker-path / scATAnno-reference vocabulary -> broad.
_ATAC_MAP: dict[str, Optional[str]] = {
    'Monocyte': 'Monocytes', 'Monocytes': 'Monocytes',
    'Memory CD4 T': 'CD4+ T', 'Naive CD4 T': 'CD4+ T',
    'Naive CD8 T': 'CD8+ T', 'Effector memory CD8 T': 'CD8+ T',
    'Central memory CD8 T': 'CD8+ T',
    'Treg': 'Regulatory T', 'MAIT': 'Unconventional T',
    'Memory B': 'Memory B', 'Memory_B': 'Memory B',
    'Naive B': 'Naive B', 'Naive_B': 'Naive B',
    'Plasma cell': 'Memory B', 'Plasma_cells': 'Memory B', 'Plasma cells': 'Memory B',
    'NK': 'NK', 'cDC': 'Dendritic Cells', 'pDC': 'Dendritic Cells',
    'Dendritic cell': 'Dendritic Cells',
}


def map_rna_broad(label: str) -> Optional[str]:
    if not label or label in _SKIP:
        return None
    return _CELLTYPIST_MAP.get(label)


def map_atac_broad(label: str) -> Optional[str]:
    if not label or label in _SKIP:
        return None
    return _ATAC_MAP.get(label)


# ── h5 readers (obs only — never materialise X) ───────────────────────────────
def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def _read_categorical(node):
    if isinstance(node, h5py.Group) and 'categories' in node and 'codes' in node:
        cats = _decode(node['categories'][:])
        codes = node['codes'][:]
        return pd.Categorical.from_codes(codes, cats).astype(object)
    return _decode(node[:])


def load_labels(h5mu_path: str) -> pd.DataFrame:
    """Both modalities already share the MultiVI obs index — no barcode join needed."""
    with h5py.File(h5mu_path, 'r') as f:
        idx = _decode(f['obs/_index'][:])
        cols = {
            'rna_cell_type': _read_categorical(f['mod/rna/obs/cell_type']),
            'atac_cell_type': _read_categorical(f['mod/atac/obs/cell_type']),
        }
        for mod in ('rna', 'atac'):
            key = f'mod/{mod}/obs/cell_type_source'
            if key in f:
                cols[f'{mod}_source'] = _read_categorical(f[key])
    return pd.DataFrame(cols, index=idx)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--h5mu', required=True, help='multivi_integrated.h5mu')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    obs = load_labels(args.h5mu)

    obs['rna_broad'] = [map_rna_broad(x) for x in obs['rna_cell_type']]
    obs['atac_broad'] = [map_atac_broad(x) for x in obs['atac_cell_type']]

    # Unmapped labels are excluded from the denominator, not counted as disagreement.
    unmapped_rna = sorted({r for r, b in zip(obs['rna_cell_type'], obs['rna_broad'])
                           if b is None})
    unmapped_atac = sorted({a for a, b in zip(obs['atac_cell_type'], obs['atac_broad'])
                            if b is None})

    both = obs['rna_broad'].notna() & obs['atac_broad'].notna()
    n_both = int(both.sum())
    n_agree = int((obs['rna_broad'] == obs['atac_broad'])[both].sum())
    l1 = n_agree / n_both if n_both else float('nan')

    # Degeneracy diagnostic — the reason this number is a wiring check, not biology.
    n_rna_lab = obs['rna_cell_type'].nunique()
    n_atac_lab = obs['atac_cell_type'].nunique()
    n_rna_broad = obs.loc[both, 'rna_broad'].nunique()
    n_atac_broad = obs.loc[both, 'atac_broad'].nunique()

    # Chance floor: agreement expected if ATAC labels were drawn from the RNA
    # marginal. With a constant ATAC label this reduces to that class's RNA share.
    p_rna = obs.loc[both, 'rna_broad'].value_counts(normalize=True)
    if n_atac_broad == 1:
        only = obs.loc[both, 'atac_broad'].iloc[0]
        chance = float(p_rna.get(only, 0.0))
    else:
        p_atac = obs.loc[both, 'atac_broad'].value_counts(normalize=True)
        chance = float(sum(p_rna.get(c, 0.0) * p_atac.get(c, 0.0) for c in BROAD_CLASSES))

    conf = pd.crosstab(obs.loc[both, 'rna_broad'], obs.loc[both, 'atac_broad'],
                       rownames=['RNA_broad'], colnames=['ATAC_broad'])
    conf.to_csv(os.path.join(args.outdir, 'L1_confusion_matrix.csv'))
    obs.to_csv(os.path.join(args.outdir, 'per_cell_results.csv'))

    summary = {
        'n_cells_multivi': int(len(obs)),
        'n_cells_scored': n_both,
        'n_agree': n_agree,
        'L1_concordance': round(l1, 4),
        'chance_floor': round(chance, 4),
        'n_rna_labels_raw': int(n_rna_lab),
        'n_atac_labels_raw': int(n_atac_lab),
        'n_rna_broad_classes': int(n_rna_broad),
        'n_atac_broad_classes': int(n_atac_broad),
        'atac_degenerate': bool(n_atac_label_degenerate := n_atac_lab == 1),
        'atac_only_label': (obs['atac_cell_type'].iloc[0] if n_atac_label_degenerate
                            else None),
        'rna_source': (obs['rna_source'].iloc[0] if 'rna_source' in obs else None),
        'atac_source': (obs['atac_source'].iloc[0] if 'atac_source' in obs else None),
        'unmapped_rna_labels': unmapped_rna,
        'unmapped_atac_labels': unmapped_atac,
    }
    with open(os.path.join(args.outdir, 'summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary, indent=2))
    print('\n=== L1 confusion (RNA broad x ATAC broad) ===')
    print(conf.to_string())
    if n_atac_label_degenerate:
        print(f"\nNOTE: ATAC produced a single label ({summary['atac_only_label']}) for "
              f"all {len(obs)} cells. L1 here is a JOIN/WIRING check, not a biological\n"
              f"      agreement rate. Expected on a chr21+chr22 subset — the ATAC marker\n"
              f"      panel is genome-wide and has almost nothing to score against.")


if __name__ == '__main__':
    main()
