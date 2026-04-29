#!/usr/bin/env python3
"""
Post-hoc QC report — re-reads peak_matrix_annotated.h5ad and renders QC
visualizations not produced inline by the main pipeline.

V1 plots:
  1. tsse_vs_nfrag_marginal.pdf      — TSS vs log10(nFragment) joint scatter
                                        with KDE marginals (FORGE-style).
  2. tsse_vs_nfrag_marginal_per_sample.pdf — same, faceted by sample if present.
  3. tsse_vs_nfrag_hexbin.pdf        — hexbin variant for very-large-N viewing.

Outputs a `post_qc_manifest.json` describing the produced plots, the obs
columns observed, and the thresholds drawn — viz-only mode reads this to
know what's available without re-deriving from the h5ad.
"""

from __future__ import annotations
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def _read_obs_only(path: str) -> pd.DataFrame:
    """Materialize ONLY the obs DataFrame (with `sample` stamped) from one of:
      - a single .h5ad file
      - a directory containing per-sample *.h5ad files
      - a glob pattern matching per-sample *.h5ad files
    The `.h5ads` AnnDataSet (snapatac2 bundle) is intentionally NOT supported
    as a direct input — its full materialization via to_adata() peaks at >100 GB
    on multisample runs. Feed the per-sample h5ad glob instead (each sample's
    h5ad already carries tsse/n_fragment from ATAC_INITIAL_QC)."""
    import glob as _glob
    p = str(path)

    if p.endswith('.h5ads'):
        sys.exit(
            f"[post_qc_report] {p} is an AnnDataSet bundle. Pass the per-sample "
            f"*.h5ad files instead (or their parent directory)."
        )

    if os.path.isdir(p):
        h5ad_files = sorted(_glob.glob(os.path.join(p, '*.h5ad')))
    elif any(c in p for c in '*?['):
        h5ad_files = sorted(_glob.glob(p))
    else:
        h5ad_files = [p]

    if not h5ad_files:
        sys.exit(f"[post_qc_report] No .h5ad files found at {p}")

    parts = []
    for f in h5ad_files:
        # backed='r' keeps .X out of memory; we only consume .obs
        ad = sc.read_h5ad(f, backed='r')
        obs = ad.obs.copy()
        if 'sample' not in obs.columns:
            obs['sample'] = os.path.splitext(os.path.basename(f))[0]
        parts.append(obs)
        ad.file.close()
    return pd.concat(parts, axis=0, ignore_index=False)


def _resolve_obs_cols(obs: pd.DataFrame) -> tuple[str, str, str | None]:
    tsse_col = next((c for c in ('tsse', 'tss_enrichment', 'TSSEnrichment') if c in obs.columns), None)
    nfrag_col = next((c for c in ('n_fragment', 'nFrags', 'n_fragments') if c in obs.columns), None)
    if not tsse_col or not nfrag_col:
        sys.exit(f"[post_qc_report] Required obs columns not found. "
                 f"Have: {list(obs.columns)}")
    sample_col = next((c for c in ('sample', 'Sample', 'batch') if c in obs.columns), None)
    return tsse_col, nfrag_col, sample_col


def _plot_marginal(df: pd.DataFrame, x_col: str, y_col: str, tsse_min: float,
                   nfrag_min: float, nfrag_max: float, out_pdf: str,
                   hue: str | None = None) -> None:
    g = sns.JointGrid(data=df, x=x_col, y=y_col, height=8, ratio=4)
    if hue and df[hue].nunique() <= 30:
        g.plot_joint(sns.scatterplot, s=4, alpha=0.35, hue=df[hue], legend='brief', edgecolor='none')
    else:
        g.plot_joint(sns.scatterplot, s=3, alpha=0.3, color='steelblue', edgecolor='none')
    g.plot_marginals(sns.kdeplot, fill=True, color='steelblue', alpha=0.55)
    g.refline(x=np.log10(max(nfrag_min, 1)), y=tsse_min, color='red', linestyle='--', linewidth=1)
    g.refline(x=np.log10(max(nfrag_max, 1)), color='red', linestyle='--', linewidth=1)
    g.set_axis_labels('log10(n_fragments)', 'TSS enrichment')
    g.figure.suptitle('Per-cell QC: TSS enrichment vs library size', y=1.02)
    g.figure.savefig(out_pdf, bbox_inches='tight')
    plt.close(g.figure)


def _plot_hexbin(df: pd.DataFrame, x_col: str, y_col: str, tsse_min: float,
                 nfrag_min: float, nfrag_max: float, out_pdf: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    hb = ax.hexbin(df[x_col], df[y_col], gridsize=80, mincnt=1, cmap='viridis')
    fig.colorbar(hb, ax=ax, label='cells / hex')
    ax.axhline(tsse_min, color='red', linestyle='--', linewidth=1)
    ax.axvline(np.log10(max(nfrag_min, 1)), color='red', linestyle='--', linewidth=1)
    ax.axvline(np.log10(max(nfrag_max, 1)), color='red', linestyle='--', linewidth=1)
    ax.set_xlabel('log10(n_fragments)')
    ax.set_ylabel('TSS enrichment')
    ax.set_title('Per-cell QC density (hexbin)')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--h5ad', required=True)
    p.add_argument('--tsse-min', type=float, default=6.0)
    p.add_argument('--nfrag-min', type=float, default=5000.0)
    p.add_argument('--nfrag-max', type=float, default=100000.0)
    p.add_argument('--outdir', default='.')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    obs = _read_obs_only(args.h5ad)
    tsse_col, nfrag_col, sample_col = _resolve_obs_cols(obs)
    n_cells = len(obs)

    df = obs[[tsse_col, nfrag_col]].copy()
    df = df.rename(columns={tsse_col: 'tsse', nfrag_col: 'n_fragment'})
    df['log10_nfrag'] = np.log10(df['n_fragment'].clip(lower=1).astype(float))
    if sample_col:
        df['_sample'] = obs[sample_col].astype(str).values

    out_paths = {}

    marg_pdf = os.path.join(args.outdir, 'tsse_vs_nfrag_marginal.pdf')
    _plot_marginal(df, 'log10_nfrag', 'tsse', args.tsse_min, args.nfrag_min,
                   args.nfrag_max, marg_pdf,
                   hue='_sample' if sample_col else None)
    out_paths['marginal'] = os.path.basename(marg_pdf)

    if sample_col:
        per_sample_pdf = os.path.join(args.outdir, 'tsse_vs_nfrag_marginal_per_sample.pdf')
        samples = sorted(df['_sample'].unique())
        ncol = min(4, len(samples))
        nrow = int(np.ceil(len(samples) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow), squeeze=False)
        for ax, s in zip(axes.flat, samples):
            sub = df[df['_sample'] == s]
            ax.scatter(sub['log10_nfrag'], sub['tsse'], s=2, alpha=0.3, color='steelblue', edgecolor='none')
            ax.axhline(args.tsse_min, color='red', linestyle='--', linewidth=0.8)
            ax.axvline(np.log10(max(args.nfrag_min, 1)), color='red', linestyle='--', linewidth=0.8)
            ax.axvline(np.log10(max(args.nfrag_max, 1)), color='red', linestyle='--', linewidth=0.8)
            ax.set_title(f'{s}  (n={len(sub):,})', fontsize=10)
            ax.set_xlabel('log10(n_fragments)'); ax.set_ylabel('TSS enrichment')
        for ax in axes.flat[len(samples):]:
            ax.set_visible(False)
        fig.tight_layout()
        fig.savefig(per_sample_pdf, bbox_inches='tight')
        plt.close(fig)
        out_paths['per_sample'] = os.path.basename(per_sample_pdf)

    hex_pdf = os.path.join(args.outdir, 'tsse_vs_nfrag_hexbin.pdf')
    _plot_hexbin(df, 'log10_nfrag', 'tsse', args.tsse_min, args.nfrag_min,
                 args.nfrag_max, hex_pdf)
    out_paths['hexbin'] = os.path.basename(hex_pdf)

    manifest = {
        'h5ad_input': os.path.basename(args.h5ad),
        'n_cells': int(n_cells),
        'obs_columns_seen': list(obs.columns),
        'tsse_column': tsse_col,
        'nfrag_column': nfrag_col,
        'sample_column': sample_col,
        'thresholds': {
            'tsse_min': args.tsse_min,
            'nfrag_min': args.nfrag_min,
            'nfrag_max': args.nfrag_max,
        },
        'plots': out_paths,
    }
    with open(os.path.join(args.outdir, 'post_qc_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"[post_qc_report] Wrote {len(out_paths)} plots + manifest to {args.outdir}")


if __name__ == '__main__':
    main()
