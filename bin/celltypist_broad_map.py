#!/usr/bin/env python3
"""
Shared CELLTYPIST_BROAD_MAP for hierarchical condensation of
CellTypist Immune_All_Low predictions to ~10 broad cell types.

Used by both RNA (plot_postscanvi.py) and ATAC (merge_annotations.py)
annotation pipelines.
"""

# Broad cell type mapping for CellTypist Immune_All_Low predictions.
# Used adaptively: only when <50% of fine-grained types pass min_cells threshold.
CELLTYPIST_BROAD_MAP = {
    # CD4+ T cells
    "Tcm/Naive helper T cells": "CD4+ T cells",
    "Tem/Effector helper T cells": "CD4+ T cells",
    "Tem/Effector helper T cells PD1+": "CD4+ T cells",
    "Follicular helper T cells": "CD4+ T cells",
    "Type 17 helper T cells": "CD4+ T cells",
    "Type 1 helper T cells": "CD4+ T cells",
    "Memory CD4+ cytotoxic T cells": "CD4+ T cells",
    "Cycling T cells": "CD4+ T cells",
    # Regulatory T cells
    "Regulatory T cells": "Regulatory T cells",
    "Treg(diff)": "Regulatory T cells",
    # CD8+ T cells
    "Tcm/Naive cytotoxic T cells": "CD8+ T cells",
    "Tem/Trm cytotoxic T cells": "CD8+ T cells",
    "Tem/Temra cytotoxic T cells": "CD8+ T cells",
    "Trm cytotoxic T cells": "CD8+ T cells",
    "CD8a/a": "CD8+ T cells",
    # Unconventional T cells
    "MAIT cells": "Unconventional T cells",
    "gamma-delta T cells": "Unconventional T cells",
    "CRTAM+ gamma-delta T cells": "Unconventional T cells",
    "NKT cells": "Unconventional T cells",
    "Double-negative thymocytes": "Unconventional T cells",
    "T(agonist)": "Unconventional T cells",
    # B cells
    "Naive B cells": "B cells",
    "Memory B cells": "B cells",
    "Age-associated B cells": "B cells",
    "Plasma cells": "B cells",
    "Plasmablasts": "B cells",
    "Germinal center B cells": "B cells",
    "Proliferative germinal center B cells": "B cells",
    "B cells": "B cells",
    # Monocytes & Macrophages
    "Classical monocytes": "Monocytes",
    "Non-classical monocytes": "Monocytes",
    "Macrophages": "Monocytes",
    "Alveolar macrophages": "Monocytes",
    "Intermediate macrophages": "Monocytes",
    "Mono-mac": "Monocytes",
    "Monocyte precursor": "Monocytes",
    "Monocytes": "Monocytes",
    # Dendritic Cells
    "DC1": "Dendritic Cells",
    "DC2": "Dendritic Cells",
    "DC3": "Dendritic Cells",
    "pDC": "Dendritic Cells",
    "Migratory DCs": "Dendritic Cells",
    # NK cells
    "CD16+ NK cells": "NK cells",
    "CD16- NK cells": "NK cells",
    "NK cells": "NK cells",
    "Cycling NK cells": "NK cells",
    # Progenitors/Other
    "ILC3": "Progenitors/Other",
    "ILC precursor": "Progenitors/Other",
    "HSC/MPP": "Progenitors/Other",
    "Megakaryocyte precursor": "Progenitors/Other",
    "Megakaryocytes/platelets": "Progenitors/Other",
    "Late erythroid": "Progenitors/Other",
    "Endothelial cells": "Progenitors/Other",
    "Epithelial cells": "Progenitors/Other",
    "Fibroblasts": "Progenitors/Other",
    "ELP": "Progenitors/Other",
}
