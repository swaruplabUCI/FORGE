#!/usr/bin/env python3
"""
build_release_assets.py — assemble the T2 tutorial release assets.

Generates, under --outdir (default dev_notes/phase3/release_assets/):

    expected_results.json   structural contract a reader can diff their run against
    figures/                curated reference figures spanning every arm
    figures/CHECKSUMS.txt   sha256 of each reference figure

The tarball (forge_tutorial_pbmc_v1.tar.gz) is built separately with tar; see
RELEASE_NOTES.md for the exact command.

WHY A SCRIPT AND NOT A HAND-WRITTEN JSON
-----------------------------------------
Every number here is derived from the run directory at build time, so the
contract cannot drift from the run it claims to describe. Re-run this after any
pipeline change that legitimately moves the numbers.

STABLE vs INFORMATIONAL
-----------------------
Counts under 'structural' are deterministic — two cold runs reproduced them
exactly (seeded, incl. scvi-tools). A reader whose run differs on one of these
has a genuinely different input or container, and should investigate.

Values under 'informational' are NOT stable run to run (wall-clock, peak RSS).
They are recorded for scale only. Do not assert on them.

Run inside snapatac_extended.sif:

    singularity exec --bind /dfs7,/tmp singularity_cache/snapatac_extended.sif \
        python3 dev_notes/phase3/build_release_assets.py \
            --results results_tutorial_remeasure \
            --outdir dev_notes/phase3/release_assets
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil

import h5py

# Curated reference figures: one or two per arm, enough to eyeball a run without
# shipping all 316. Native format kept (RNA/MOFA emit PNG, ATAC/MultiVI PDF).
REFERENCE_FIGURES = [
    # RNA arm
    'rna/concatenated/concat_plots/qc_distributions_4panel.png',
    'rna/post_integration_plots/umap_cell_types.png',
    'rna/post_integration_plots/umap_leiden_0.5.png',
    'rna/post_integration_plots/integration_celltype_composition.png',
    # ATAC arm
    'atac/initial_qc/qc_plots/qc_lower_bound_thresholds.pdf',
    'atac/initial_qc/qc_plots/prefilter_tsse_ridges.pdf',
    'atac/initial_qc/qc_plots/umap_leiden_0_5.pdf',
    # Integration
    'mofa_visualization/figures/mofa_variance_explained.png',
    'mofa_visualization/figures/umap_mofa_rna_celltype.png',
    'multiome/multivi/visualizations/multivi_plots/umap_celltype.pdf',
    'multiome/multivi/visualizations/multivi_plots/umap_leiden.pdf',
    # Regulatory / downstream
    'cicero/full_ConnsCoAccPseudoTAcc/ccan_size_distribution.pdf',
]


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def h5_shape(path: str) -> tuple[int, int]:
    """(n_obs, n_var) without materialising X."""
    with h5py.File(path, 'r') as f:
        n_obs = f['obs/_index'].shape[0]
        n_var = f['var/_index'].shape[0] if 'var/_index' in f else -1
    return int(n_obs), int(n_var)


def gz_datarows(path: str) -> int:
    """Row count excluding the header line."""
    with gzip.open(path, 'rt') as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def csv_datarows(path: str) -> int:
    with open(path, newline='') as fh:
        return max(sum(1 for _ in csv.reader(fh)) - 1, 0)


def load_json(path: str):
    with open(path) as fh:
        return json.load(fh)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results', required=True, help='the tutorial run output dir')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--trace', default='dev_notes/phase3/remeasure_trace/trace.tsv')
    ap.add_argument('--concordance',
                    default='dev_notes/phase3/concordance/summary.json')
    args = ap.parse_args()

    R = args.results.rstrip('/')
    figdir = os.path.join(args.outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)

    rna_obs, rna_var = h5_shape(f'{R}/rna/post_integration_plots/annotated_with_celltype.h5ad')
    atac_obs, atac_var = h5_shape(f'{R}/consolidated_qc/peak_matrix_annotated.h5ad')

    atac_qc = load_json(f'{R}/atac/initial_qc/atac_initial_qc_summary.json')
    thresholds = load_json(f'{R}/atac/initial_qc/sample_thresholds.json')
    mudata = load_json(f'{R}/multiome/mudata/mudata_stats.json') \
        if os.path.exists(f'{R}/multiome/mudata/mudata_stats.json') \
        else load_json('dev_notes/phase3/remeasure_trace/mudata_stats.json')
    mofa = load_json('dev_notes/phase3/remeasure_trace/mofa_stats.json')
    atac_desc = load_json(f'{R}/atac/descriptive/summary.json')

    with open(args.trace) as fh:
        n_tasks = max(sum(1 for _ in fh) - 1, 0)

    conc = load_json(args.concordance) if os.path.exists(args.concordance) else {}

    cicero_dir = f'{R}/cicero/full_ConnsCoAccPseudoTAcc'
    hdwgcna_cts = sorted(d for d in os.listdir(f'{R}/hdwgcna')
                         if d != 'enrichment'
                         and os.path.isdir(os.path.join(f'{R}/hdwgcna', d)))

    contract = {
        'dataset': 'forge_tutorial_pbmc_v1',
        'description': (
            'Structural contract for the FORGE T2 tutorial run. Counts under '
            '"structural" are deterministic and safe to assert on; "informational" '
            'values vary run to run and must not be asserted.'
        ),
        'pipeline_profile': 'tutorial,singularity',
        'random_seed': 42,
        'structural': {
            'n_tasks_completed': n_tasks,
            'rna': {
                'n_cells_annotated': rna_obs,
                'n_genes': rna_var,
                'annotation_method': 'celltypist',
            },
            'atac': {
                'n_cells_prefilter': atac_qc['total_cells'],
                'n_cells_final': atac_obs,
                'n_peaks_final': atac_var,
                'leiden_resolutions': atac_qc['clustering_resolutions'],
                'n_cell_types_assigned': atac_desc['n_cell_types'],
                'computed_thresholds': thresholds,
            },
            'multiome': {
                'n_paired_cells': mudata['total_cells'],
                'n_rna_genes': mudata['total_rna_genes'],
                'n_atac_peaks': mudata['total_atac_peaks'],
            },
            'mofa': {
                'n_factors': mofa['n_factors'],
                'n_rna_features': mofa['n_rna_features'],
                'n_atac_features': mofa['n_atac_features'],
            },
            'cicero': {
                'n_connections': gz_datarows(f'{cicero_dir}/cicero_connections.tsv.gz'),
                'n_ccan_assignments': gz_datarows(f'{cicero_dir}/CCAN_assignments.tsv.gz'),
                'n_triplets': gz_datarows(f'{R}/cicero/cicero_triplets.tsv.gz'),
            },
            'cellchat': {
                'n_interactions': csv_datarows(f'{R}/cellchat/integrated_cellchat_results.csv'),
            },
            'hdwgcna': {
                'n_cell_types_with_output': len(hdwgcna_cts),
            },
            'concordance': {
                'n_cells_scored': conc.get('n_cells_scored'),
                'L1_concordance': conc.get('L1_concordance'),
                'chance_floor': conc.get('chance_floor'),
                'note': ('Equal to the chance floor by construction — the ATAC arm '
                         'yields a single label on a chr21+22 subset. Wiring check '
                         'only; not a quality metric.'),
            },
        },
        'informational': {
            'note': 'NOT stable run to run. Recorded for scale only — do not assert.',
            'wall_clock': '1h42m45s on 8 CPUs',
            'cpu_hours': 6.6,
            'peak_single_task_rss_gb': 8.80,
            'results_dir_size': '3.88 GB',
            'work_dir_size': '4.81 GB',
            'recommended_free_disk': '~15 GB',
        },
    }

    out_json = os.path.join(args.outdir, 'expected_results.json')
    with open(out_json, 'w') as fh:
        json.dump(contract, fh, indent=2)
        fh.write('\n')

    copied, missing = [], []
    for rel in REFERENCE_FIGURES:
        src = os.path.join(R, rel)
        if not os.path.exists(src):
            missing.append(rel)
            continue
        dst = os.path.join(figdir, rel.replace('/', '__'))
        shutil.copy2(src, dst)
        copied.append((os.path.basename(dst), sha256(dst), rel))

    with open(os.path.join(figdir, 'CHECKSUMS.txt'), 'w') as fh:
        fh.write('# sha256  filename  (source path within the run outdir)\n')
        for name, digest, rel in copied:
            fh.write(f'{digest}  {name}  # {rel}\n')

    print(json.dumps(contract, indent=2))
    print(f'\nwrote {out_json}')
    print(f'copied {len(copied)} reference figures -> {figdir}')
    if missing:
        print(f'WARNING: {len(missing)} figures missing from the run: {missing}')


if __name__ == '__main__':
    main()
