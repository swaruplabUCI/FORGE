// ============================================================================
// EXTRACT_CHROMVAR_MOTIFS — FIX-23b
// Per-cell-type TF motif extraction from GPU ChromVAR deviation matrix.
//
// For each cell type in the deviation matrix, computes mean |z-score| per
// TF motif, selects top N, and outputs a structured JSON with per-cell-type
// TF lists.  This drives the DISCOVERY mode of scPRINTER footprinting:
// each cell type gets its own set of enriched TFs as footprinting targets.
//
// ChromVAR var_names are tab-separated: "MA0001.2\tSPI1".
// The TF gene symbol is the LAST field after splitting on "\t".
// ============================================================================

process EXTRACT_CHROMVAR_MOTIFS {

    tag "extract_chromvar_motifs_per_celltype"
    label 'process_low'

    publishDir "${params.outdir}/chromvar/motif_extraction", mode: 'copy'

    input:
    path chromvar_dev   // chromvar_matrix_deviations.h5ad from GPU_CHROMVAR
    val top_n           // top N TFs per cell type
    val min_zscore      // minimum mean |z-score| threshold
    val cell_type_col   // obs column for cell type (marker: 'cell_type', celltypist: 'celltypist_prediction')

    output:
    path "per_celltype_motifs.json", emit: motif_list
    path "motif_report.txt",         emit: report

    script:
    """
    #!/usr/bin/env python3
    \"\"\"
    FIX-23b: Per-cell-type TF motif extraction from ChromVAR deviations.

    For each cell type, identifies the top N enriched TF motifs by
    mean |z-score|.  Outputs JSON with per-cell-type TF lists and a
    union of all unique TFs for coordinate resolution.
    \"\"\"
    import anndata as ad
    import numpy as np
    import pandas as pd
    import json
    import sys

    print("FIX-23b: Per-cell-type ChromVAR motif extraction", flush=True)

    # ---- Load ChromVAR deviations ----
    adata = ad.read_h5ad('${chromvar_dev}')
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} motifs", flush=True)

    # ---- Verify cell type column ----
    ct_col = '${cell_type_col}'
    if ct_col not in adata.obs.columns:
        raise ValueError(
            f"ChromVAR deviation matrix missing '{ct_col}' in .obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    cell_types = sorted(adata.obs[ct_col].unique().tolist())
    print(f"  Cell types ({len(cell_types)}): {cell_types}", flush=True)

    # ---- Parse TF names from tab-separated var_names ----
    raw_names = adata.var_names.tolist()
    tf_names = []
    motif_ids = []
    for name in raw_names:
        parts = str(name).split('\\t')
        if len(parts) >= 2:
            motif_ids.append(parts[0])
            tf_names.append(parts[-1])
        else:
            motif_ids.append(name)
            tf_names.append(name)

    print(f"  Example var_names: {raw_names[:3]}", flush=True)
    print(f"  Parsed TF names:  {tf_names[:3]}", flush=True)

    # ---- Dense data matrix ----
    if hasattr(adata.X, 'toarray'):
        data = adata.X.toarray()
    else:
        data = np.array(adata.X)

    top_n = int(${top_n})
    min_zscore = float(${min_zscore})

    # ---- Per-cell-type extraction ----
    result = {"cell_types": {}, "all_unique_tfs": [], "params": {}}
    all_tfs_union = set()
    report_lines = []

    report_lines.append("ChromVAR Per-Cell-Type Motif Extraction (FIX-23b)")
    report_lines.append("=" * 65)
    report_lines.append(f"Total motifs: {len(raw_names)}")
    report_lines.append(f"Cell types:   {len(cell_types)}")
    report_lines.append(f"Top N per CT: {top_n}")
    report_lines.append(f"Min z-score:  {min_zscore}")
    report_lines.append("")

    for ct in cell_types:
        ct_mask = (adata.obs[ct_col] == ct).values
        n_cells = int(ct_mask.sum())

        if n_cells == 0:
            print(f"  WARNING: No cells for cell_type '{ct}', skipping", flush=True)
            continue

        # Subset data for this cell type
        ct_data = data[ct_mask, :]

        # Compute per-motif mean |z-score| for this cell type
        mean_abs_z = np.nanmean(np.abs(ct_data), axis=0)
        variances = np.nanvar(ct_data, axis=0)

        # Build motif table
        motif_table = pd.DataFrame({
            'raw_name': raw_names,
            'motif_id': motif_ids,
            'tf_name': tf_names,
            'mean_abs_zscore': mean_abs_z,
            'variance': variances,
        })

        # Deduplicate: keep highest-variance motif per TF
        motif_dedup = (
            motif_table
            .sort_values('variance', ascending=False)
            .drop_duplicates(subset='tf_name', keep='first')
            .sort_values('mean_abs_zscore', ascending=False)
        )

        # Filter by z-score threshold and take top N
        filtered = motif_dedup[motif_dedup['mean_abs_zscore'] >= min_zscore].head(top_n)

        ct_tfs = filtered['tf_name'].tolist()
        all_tfs_union.update(ct_tfs)

        result["cell_types"][ct] = {
            "n_cells": n_cells,
            "top_tfs": ct_tfs,
            "motif_details": [
                {
                    "tf": row['tf_name'],
                    "motif_id": row['motif_id'],
                    "mean_abs_zscore": round(float(row['mean_abs_zscore']), 4),
                    "variance": round(float(row['variance']), 4),
                }
                for _, row in filtered.iterrows()
            ]
        }

        # Report
        report_lines.append(f"--- {ct} ({n_cells} cells) ---")
        for i, (_, row) in enumerate(filtered.iterrows(), 1):
            report_lines.append(
                f"  {i:3d}. {row['tf_name']:20s}  ({row['motif_id']:12s})  "
                f"z={row['mean_abs_zscore']:.4f}  var={row['variance']:.4f}"
            )
        report_lines.append("")

        print(f"  {ct}: {n_cells} cells, top {len(ct_tfs)} TFs: {ct_tfs[:3]}...", flush=True)

    # ---- Union of all TFs ----
    result["all_unique_tfs"] = sorted(all_tfs_union)
    result["params"] = {
        "top_n_per_celltype": top_n,
        "min_zscore": min_zscore,
        "n_cell_types": len(cell_types),
        "total_unique_tfs": len(all_tfs_union),
    }

    # ---- Write outputs ----
    with open('per_celltype_motifs.json', 'w') as f:
        json.dump(result, f, indent=2)

    report_lines.append(f"Total unique TFs across all cell types: {len(all_tfs_union)}")
    report_lines.append(f"Union TF list: {sorted(all_tfs_union)}")

    with open('motif_report.txt', 'w') as f:
        f.write('\\n'.join(report_lines) + '\\n')

    print(f"FIX-23b: Extraction complete — {len(all_tfs_union)} unique TFs across {len(cell_types)} cell types", flush=True)
    """
}
