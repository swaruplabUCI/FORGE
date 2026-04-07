#!/usr/bin/env python3
"""
Initial ATAC-seq QC with UNIFORM filtering on list of AnnData objects
Corrected version based on successful test_simple.py workflow
"""

import snapatac2 as snap
import numpy as np
import pandas as pd
import json
import argparse
import os
import logging
import contextlib
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import re

# Configure logging
logging.getLogger().setLevel(logging.WARNING)
for name in ("kaleido", "choreographer"):
    logging.getLogger(name).setLevel(logging.ERROR)

@contextlib.contextmanager
def quiet_plotly_export():
    """Silence plotly export noise"""
    devnull = open(os.devnull, "w")
    try:
        with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
            yield
    finally:
        devnull.close()

def ridge_plot(df, value_col, group_col="Sample", bw_adjust=0.9, height_per_group=1.0):
    """Create ArchR-style ridge plots"""
    groups = list(pd.unique(df[group_col]))
    n = len(groups)
    fig, axes = plt.subplots(n, 1, figsize=(10, height_per_group*n), sharex=True)
    if n == 1:
        axes = [axes]
    
    x_min = float(df[value_col].min())
    x_max = float(df[value_col].max())
    
    for ax, g in zip(axes, groups):
        sub = df[df[group_col] == g]
        sns.kdeplot(data=sub, x=value_col, fill=True, bw_adjust=bw_adjust, 
                   clip=(x_min, x_max), ax=ax, color="#4C72B0")
        ax.set_ylabel(g)
        ax.set_yticks([])
        ax.grid(False)
    
    axes[-1].set_xlabel(value_col)
    plt.tight_layout()
    return fig

def build_qc_df(ad_list, sample_names, max_cells=20000, random_state=0):
    """Build long-form DataFrame for visualization"""
    rows = []
    rng = np.random.default_rng(random_state)
    
    for name, ad in zip(sample_names, ad_list):
        ts = ad.obs["tsse"].to_numpy()
        nf = ad.obs["n_fragment"].to_numpy()
        
        mask = np.isfinite(ts) & np.isfinite(nf) & (nf > 0)
        ts = ts[mask]
        nf = nf[mask]
        
        if ts.size > max_cells:
            idx = rng.choice(ts.size, size=max_cells, replace=False)
            ts = ts[idx]
            nf = nf[idx]
        
        df_sample = pd.DataFrame({
            "Sample": name,
            "TSSEnrichment": ts,
            "log10_nFrags": np.log10(nf),
        })
        rows.append(df_sample)
    
    return pd.concat(rows, axis=0, ignore_index=True)

def main():
    parser = argparse.ArgumentParser(description='Initial ATAC-seq QC pipeline')
    parser.add_argument('--fragment_files', nargs='+', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--species', default='human', choices=['human', 'mouse'])
    parser.add_argument('--genome_build', default=None, help='Genome build (e.g. hg38, mm10, mm39). Auto-detected from species if not set.')
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--min_fragments', type=int, default=1000)
    parser.add_argument('--min_counts', type=int, default=5000)
    parser.add_argument('--min_tsse', type=float, default=6)
    parser.add_argument('--max_counts', type=int, default=100000)
    parser.add_argument('--n_features', type=int, default=50000)
    parser.add_argument('--batch_correction', default='mnn', choices=['mnn', 'harmony', 'none'])
    parser.add_argument('--batch_key', default='sample')
    parser.add_argument('--clustering_resolutions', nargs='+', type=float, 
                       default=[0.5, 5.0])
    parser.add_argument('--peak_fdr', type=float, default=0.01)
    parser.add_argument('--tempdir', default='/tmp')
    parser.add_argument('--n_jobs', type=int, default=1)  # Changed default to 1
    args = parser.parse_args()
    
    # CRITICAL FIX: Disable HDF5 file locking for DFS/NFS
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    
    # Setup
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True)
    plot_dir = out_dir / "qc_plots"
    plot_dir.mkdir(exist_ok=True)
    
    # Load metadata
    metadata = pd.read_csv(args.metadata)

    # Create a unique QC sample ID per (sample_id, batch)
    # e.g. L1_Donor_0_july, L1_Donor_0_november
    metadata['qc_sample'] = (
        metadata['sample_id'].astype(str) + '_' + metadata['batch'].astype(str)
    )

    sample_names = metadata['qc_sample'].tolist()
    
    # Define genome — use explicit build if provided, else default per species
    genome_map = {'hg38': snap.genome.hg38, 'hg19': snap.genome.hg19,
                  'mm10': snap.genome.mm10, 'mm39': snap.genome.mm39,
                  'GRCm39': snap.genome.mm39, 'GRCh38': snap.genome.hg38}
    if args.genome_build and args.genome_build in genome_map:
        genome = genome_map[args.genome_build]
    else:
        genome = snap.genome.hg38 if args.species == 'human' else snap.genome.mm10
    
    # Match fragment files to metadata order
    frag_bases = metadata['fragment_file'].tolist()
    file_map = {}
    for f in args.fragment_files:
        base = Path(f).name
        # Key 1: strip .bed.gz or .tsv.gz  →  keeps _barcode_sorted
        # FIX-37b: Support both BD (.bed.gz) and 10x (.tsv.gz) fragment formats
        key_full = re.sub(r'\.(bed|tsv)\.gz$', '', base)
        file_map[key_full] = f
        # Key 2: also strip _(barcode|coord)_sorted suffix  →  legacy compat
        key_short = re.sub(r'_(barcode|coord)_sorted$', '', key_full)
        if key_short != key_full:
            file_map[key_short] = f

    ordered_fragment_files = []
    for b in frag_bases:
        # Strip .bed.gz or .tsv.gz from manifest value too, if present
        # FIX-37b: Support both BD (.bed.gz) and 10x (.tsv.gz) fragment formats
        b_key = re.sub(r'\.(bed|tsv)\.gz$', '', b)
        if b_key not in file_map:
            available = sorted(file_map.keys())[:10]
            raise FileNotFoundError(
                f"Fragment file for '{b}' not found. "
                f"Available keys (first 10): {available}"
            )
        ordered_fragment_files.append(file_map[b_key])
    
    out_files = [out_dir / f"{sample}.h5ad" for sample in sample_names]
    dataset_path = out_dir / "atac_initial_qc.h5ads"
    
    # ========================================================================
    # STEP 1: IMPORT FRAGMENTS (returns list of AnnData objects)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 1: Importing fragments")
    print(f"{'='*70}")
    
    GMWM_all = snap.pp.import_fragments(
        ordered_fragment_files,
        file=[str(p) for p in out_files],  # Backed mode
        chrom_sizes=genome,
        min_num_fragments=args.min_fragments,
        tempdir=args.tempdir,
        n_jobs=args.n_jobs,
        sorted_by_barcode=False
    )
    
    print(f"Imported {len(GMWM_all)} samples")
    
    # ========================================================================
    # STEP 2: CALCULATE QC METRICS
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 2: Calculating QC metrics")
    print(f"{'='*70}")
    
    print("Calculating fragment size distribution...")
    snap.metrics.frag_size_distr(GMWM_all, add_key='preQCfrag_size_distr', n_jobs=args.n_jobs)
    
    print("Calculating TSS enrichment...")
    snap.metrics.tsse(GMWM_all, genome, n_jobs=args.n_jobs)
    
    # ========================================================================
    # STEP 3: GENERATE PRE-FILTER VISUALIZATIONS
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 3: Generating pre-filter visualizations")
    print(f"{'='*70}")
    
    # snapATAC2 native plots
    for i, ad in enumerate(GMWM_all, start=1):
        fig = snap.pl.tsse(ad, show=False, out_file=None, interactive=False)
        with quiet_plotly_export():
            fig.write_image(str(plot_dir / f"prefilter_tsse_s{i}.pdf"))
        
        fig = snap.pl.frag_size_distr(ad, use_rep='preQCfrag_size_distr', 
                                      show=False, out_file=None, interactive=False)
        with quiet_plotly_export():
            fig.write_image(str(plot_dir / f"prefilter_fsd_s{i}.pdf"))
        
        fig.update_yaxes(type="log")
        with quiet_plotly_export():
            fig.write_image(str(plot_dir / f"prefilter_fsd_log_s{i}.pdf"))
    
    # ArchR-style plots
    qc_df = build_qc_df(GMWM_all, sample_names)
    
    fig = ridge_plot(qc_df, "TSSEnrichment")
    fig.savefig(plot_dir / "prefilter_tsse_ridges.pdf", dpi=300)
    plt.close()
    
    fig = ridge_plot(qc_df, "log10_nFrags")
    fig.savefig(plot_dir / "prefilter_fragments_ridges.pdf", dpi=300)
    plt.close()

    # ---- Upper bound QC threshold plots ----
    print("Generating upper/lower bound QC threshold plots...")

    # Compute nucleosome signal per sample (mono-nucleosome / nucleosome-free ratio)
    nuc_signal_rows = []
    for name, ad in zip(sample_names, GMWM_all):
        if 'preQCfrag_size_distr' in ad.uns:
            fsd = ad.uns['preQCfrag_size_distr']
            # fsd is typically a 2D array: rows=cells, cols=fragment sizes (1bp bins)
            if hasattr(fsd, 'shape') and len(fsd.shape) == 2:
                nf_band = fsd[:, :147].sum(axis=1)      # nucleosome-free (<147bp)
                mono_band = fsd[:, 147:294].sum(axis=1)  # mono-nucleosome (147-294bp)
                nf_band = np.where(nf_band == 0, 1, nf_band)
                nuc_sig = mono_band / nf_band
            else:
                nuc_sig = np.full(ad.n_obs, np.nan)
        else:
            nuc_sig = np.full(ad.n_obs, np.nan)
        nuc_signal_rows.append(pd.DataFrame({
            'Sample': name,
            'nucleosome_signal': nuc_sig,
            'tsse': ad.obs['tsse'].to_numpy(),
            'n_fragment': ad.obs['n_fragment'].to_numpy(),
        }))
    nuc_df = pd.concat(nuc_signal_rows, ignore_index=True)
    nuc_df = nuc_df[np.isfinite(nuc_df['tsse']) & np.isfinite(nuc_df['n_fragment']) & (nuc_df['n_fragment'] > 0)]
    nuc_df['log10_nFrags'] = np.log10(nuc_df['n_fragment'])

    # Upper bound: Scatter (fragments vs TSS enrichment, colored by TSS)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sc_ax = axes[0]
    scatter = sc_ax.scatter(nuc_df['log10_nFrags'], nuc_df['tsse'],
                            c=nuc_df['tsse'], cmap='RdYlGn', s=1, alpha=0.3, rasterized=True)
    sc_ax.axhline(y=args.min_tsse, color='red', linestyle='--', label=f'min TSS={args.min_tsse}')
    sc_ax.axvline(x=np.log10(args.min_counts), color='blue', linestyle='--', label=f'min frags={args.min_counts}')
    sc_ax.axvline(x=np.log10(args.max_counts), color='blue', linestyle=':', label=f'max frags={args.max_counts}')
    sc_ax.set_xlabel('log10(Fragment Count)')
    sc_ax.set_ylabel('TSS Enrichment')
    sc_ax.set_title('Upper Bound QC: Fragments vs TSS')
    sc_ax.legend(fontsize=7)
    plt.colorbar(scatter, ax=sc_ax, label='TSS Enrichment')

    # Upper bound: Violin of TSS enrichment
    sns.violinplot(data=nuc_df, x='Sample', y='tsse', ax=axes[1], inner='quartile', cut=0)
    axes[1].axhline(y=args.min_tsse, color='red', linestyle='--')
    axes[1].set_title('TSS Enrichment Distribution')
    axes[1].tick_params(axis='x', rotation=45)

    # Upper bound: Violin of nucleosome signal
    if nuc_df['nucleosome_signal'].notna().sum() > 0:
        sns.violinplot(data=nuc_df, x='Sample', y='nucleosome_signal', ax=axes[2], inner='quartile', cut=0)
        axes[2].axhline(y=2.0, color='red', linestyle='--', label='threshold=2')
        axes[2].set_title('Nucleosome Signal Distribution')
        axes[2].tick_params(axis='x', rotation=45)
        axes[2].legend(fontsize=7)
    else:
        axes[2].text(0.5, 0.5, 'Nucleosome signal\nnot available', ha='center', va='center', transform=axes[2].transAxes)
        axes[2].set_title('Nucleosome Signal (N/A)')

    plt.tight_layout()
    fig.savefig(plot_dir / "qc_upper_bound_thresholds.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: qc_upper_bound_thresholds.pdf")

    # ---- Lower bound QC threshold plots ----
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Joint density: TSS vs log-fragments
    axes[0].scatter(nuc_df['log10_nFrags'], nuc_df['tsse'], s=1, alpha=0.1, c='steelblue', rasterized=True)
    try:
        sns.kdeplot(data=nuc_df, x='log10_nFrags', y='tsse', ax=axes[0],
                    levels=5, color='darkblue', linewidths=0.8)
    except Exception:
        pass
    axes[0].axhline(y=args.min_tsse, color='red', linestyle='--', linewidth=0.8)
    axes[0].axvline(x=np.log10(args.min_counts), color='red', linestyle='--', linewidth=0.8)
    axes[0].set_xlabel('log10(Fragment Count)')
    axes[0].set_ylabel('TSS Enrichment')
    axes[0].set_title('Lower Bound QC: Density')

    # Histogram of low-fragment tail
    low_frags = nuc_df[nuc_df['n_fragment'] < np.percentile(nuc_df['n_fragment'], 25)]
    axes[1].hist(low_frags['n_fragment'], bins=50, color='steelblue', alpha=0.7, edgecolor='white')
    axes[1].axvline(x=args.min_counts, color='red', linestyle='--', label=f'min_counts={args.min_counts}')
    axes[1].set_xlabel('Fragment Count')
    axes[1].set_ylabel('Number of Cells')
    axes[1].set_title('Low-Fragment Tail Distribution')
    axes[1].legend(fontsize=7)

    # Scatter of cells below threshold, colored by TSS
    below = nuc_df[nuc_df['n_fragment'] < args.min_counts * 2]
    if len(below) > 0:
        scatter = axes[2].scatter(below['n_fragment'], below['tsse'],
                                  c=below['tsse'], cmap='RdYlGn', s=3, alpha=0.5, rasterized=True)
        axes[2].axvline(x=args.min_counts, color='red', linestyle='--')
        axes[2].axhline(y=args.min_tsse, color='red', linestyle='--')
        axes[2].set_xlabel('Fragment Count')
        axes[2].set_ylabel('TSS Enrichment')
        axes[2].set_title('Cells Near Lower Threshold')
        plt.colorbar(scatter, ax=axes[2], label='TSS Enrichment')
    else:
        axes[2].text(0.5, 0.5, 'No cells near\nlower threshold', ha='center', va='center', transform=axes[2].transAxes)

    plt.tight_layout()
    fig.savefig(plot_dir / "qc_lower_bound_thresholds.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: qc_lower_bound_thresholds.pdf")

    # ========================================================================
    # STEP 4: APPLY UNIFORM FILTERING TO LIST
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 4: Applying uniform filtering to list of AnnData objects")
    print(f"{'='*70}")
    
    print(f"Applying uniform thresholds:")
    print(f"  min_counts: {args.min_counts}")
    print(f"  min_tsse: {args.min_tsse}")
    print(f"  max_counts: {args.max_counts}")
    
    # Apply filtering directly to the list (operates on all samples)
    snap.pp.filter_cells(GMWM_all, min_counts=args.min_counts, min_tsse=args.min_tsse, max_counts=args.max_counts)
    
    total_cells = sum(ad.n_obs for ad in GMWM_all)
    print(f"\nTotal cells after filtering: {total_cells}")
    
    # ========================================================================
    # STEP 5: ADD TILE MATRIX (operates on list directly)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 5: Adding tile matrix")
    print(f"{'='*70}")

    snap.pp.add_tile_matrix(GMWM_all, bin_size=5000)

    # ========================================================================
    # STEP 6: FEATURE SELECTION (operates on list directly)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 6: Selecting features")
    print(f"{'='*70}")

    snap.pp.select_features(GMWM_all, n_features=args.n_features)

    # ========================================================================
    # STEP 7: SKIP DOUBLET DETECTION (INITIAL QC ONLY)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 7: Skipping doublet detection for initial QC")
    print(f"{'='*70}")
    print("Note: Doublet detection will be performed after sample-specific filtering")

    # ========================================================================
    # STEP 8: FILES ALREADY SAVED IN BACKED MODE
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 8: Individual sample files already saved in backed mode")
    print(f"{'='*70}")
    
    for i, (adata, out_file) in enumerate(zip(GMWM_all, out_files)):
        print(f"  {sample_names[i]}: {adata.n_obs} cells -> {out_file}")
        # Close file handles
        if hasattr(adata, 'file') and adata.file is not None:
            adata.file.close()

    # ========================================================================
    # STEP 9: CREATE ANNDATASET FROM SAVED FILES
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 9: Creating AnnDataSet")
    print(f"{'='*70}")
    
    data = snap.AnnDataSet(
        adatas=[(name, str(path)) for name, path in zip(sample_names, out_files)],
        filename=str(dataset_path),
        add_key='sample'
    )
    
    # CRITICAL: Create unique cell IDs and verify (from lab code Page 1)
    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids
    
    # CRITICAL: Assert uniqueness (from lab code Page 1)
    assert data.n_obs == np.unique(data.obs_names).size, "Cell IDs are not unique!"
    
    print(f"AnnDataSet created: {data.n_obs} cells")
    print(f"Verified unique cell IDs: {np.unique(data.obs_names).size}")

    # --------------------------------------------------------------------
    # Attach 'batch' (month) from metadata so it lines up with AnnDataSet
    # --------------------------------------------------------------------
    # metadata['qc_sample'] was constructed as sample_id + '_' + batch and
    # sample_names = metadata['qc_sample'], which are used as 'sample' in data.obs.
    meta_map = metadata.set_index('qc_sample')

    if 'batch' not in meta_map.columns:
        raise KeyError("'batch' column not found in metadata; cannot annotate AnnDataSet with batch.")

    # Map the month-level batch label (e.g. june/july/nov) onto every cell via its sample ID.
    # `data.obs['sample']` is not a pandas Series, so we cannot use `.map()` directly.
    batch_map = meta_map['batch'].to_dict()
    data.obs['batch'] = [batch_map.get(s, None) for s in data.obs['sample']]

    # ========================================================================
    # STEP 10: FEATURE SELECTION ON ANNDATASET (CRITICAL!)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 10: Selecting features on AnnDataSet")
    print(f"{'='*70}")
    
    snap.pp.select_features(data, n_features=args.n_features)

    # ========================================================================
    # STEP 11: SPECTRAL EMBEDDING (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 11: Computing spectral embedding")
    print(f"{'='*70}")
    
    snap.tl.spectral(data)  # Uses AnnDataSet
    
    # ========================================================================
    # STEP 12: BATCH CORRECTION (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print(f"STEP 12: Batch correction ({args.batch_correction})")
    print(f"{'='*70}")
    
    n_samples = len(sample_names)
    if n_samples <= 1:
        print(f"  Single sample detected ({n_samples}) — automatically skipping batch correction")
        use_rep = "X_spectral"
    elif args.batch_correction == 'none':
        print("  Skipping batch correction (batch_correction='none')")
        use_rep = "X_spectral"
    elif args.batch_correction == 'harmony':
        snap.pp.harmony(data, batch=args.batch_key, max_iter_harmony=20)
        use_rep = "X_spectral_harmony"
    else:  # mnn
        snap.pp.mnc_correct(data, batch=args.batch_key)
        use_rep = "X_spectral_mnn"
    
    # ========================================================================
    # STEP 13: UMAP (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 13: Computing UMAP")
    print(f"{'='*70}")
    
    snap.tl.umap(data, use_rep=use_rep)
    
    # ========================================================================
    # STEP 14: KNN AND CLUSTERING (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 14: Building KNN graph and clustering")
    print(f"{'='*70}")
    
    snap.pp.knn(data, use_rep=use_rep)
    
    for res in args.clustering_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        snap.tl.leiden(data, resolution=res, key_added=key)
        n_clusters = len(set(data.obs[key]))
        snap.pl.umap(data, color=key, show=False, 
                     out_file=str(plot_dir / f"umap_{key}.pdf"), 
                     height=500, interactive=False)
        print(f"  Resolution {res}: {n_clusters} clusters")
    
    # ========================================================================
    # STEP 17: GENERATE POST-PROCESSING VISUALIZATIONS
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 17: Generating final visualizations")
    print(f"{'='*70}")
    
    # UMAP plots
    for res in args.clustering_resolutions:
        key = f"leiden_{str(res).replace('.', '_')}"
        fig = snap.pl.umap(data, color=key, show=False, out_file=None, 
                          interactive=False, height=500)
        with quiet_plotly_export():
            fig.write_image(str(plot_dir / f"umap_{key}.pdf"))
    
    # UMAP by sample
    fig = snap.pl.umap(data, color='sample', show=False, out_file=None, 
                      interactive=False, height=500)
    with quiet_plotly_export():
        fig.write_image(str(plot_dir / "umap_by_sample.pdf"))

    # UMAP by batch (from metadata 'batch' column)
    fig = snap.pl.umap(
        data,
        color="batch",
        show=False,
        out_file=None,
        interactive=False,
        height=500,
    )
    with quiet_plotly_export():
        fig.write_image(str(plot_dir / "umap_by_batch.pdf"))
    
    # ========================================================================
    # STEP 18: SAVE SUMMARY AND CLOSE
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 18: Saving summary stats and closing")
    print(f"{'='*70}")
    
    sample_stats = []
    for name, ad in zip(sample_names, GMWM_all):
        n_frag = ad.obs['n_fragment'].to_numpy()
        tsse   = ad.obs['tsse'].to_numpy()

        frag_q05 = float(np.quantile(n_frag, 0.05))
        frag_q95 = float(np.quantile(n_frag, 0.95))
        tsse_q05 = float(np.quantile(tsse, 0.05))
        tsse_q95 = float(np.quantile(tsse, 0.95))

        # Simulate your intended per-sample filter:
        #   n_fragment in [frag_q05, frag_q95] AND tsse >= tsse_q05
        mask = (
            np.isfinite(n_frag) & np.isfinite(tsse) &
            (n_frag >= frag_q05) & (n_frag <= frag_q95) &
            (tsse >= tsse_q05)
        )
        n_cells_after_filter = int(mask.sum())

        stats = {
            'sample': name,
            'n_cells': int(ad.n_obs),
            'median_fragments': float(np.median(n_frag)),
            'median_tsse': float(np.median(tsse)),
            'fragments_q05': frag_q05,
            'fragments_q95': frag_q95,
            'tsse_q05': tsse_q05,
            'tsse_q95': tsse_q95,
            # NEW: how many cells would survive the proposed thresholds
            'n_cells_after_frag_5_95_and_tsse_q05': n_cells_after_filter,
            # NEW: min/max in case we need “keep all cells” thresholds
            'fragments_min': float(n_frag.min()),
            'fragments_max': float(n_frag.max()),
            'tsse_min': float(tsse.min()),
        }
        sample_stats.append(stats)

    stats_df = pd.DataFrame(sample_stats)
    stats_df.to_csv(out_dir / 'sample_statistics.csv', index=False)

    summary = {
        'n_samples': len(sample_names),
        'total_cells': data.n_obs,
        'samples': sample_names,
        'species': args.species,
        'batch_correction': args.batch_correction,
        'n_features': args.n_features,
        'filtering_thresholds': {
            'min_counts': args.min_counts,
            'min_tsse': args.min_tsse,
            'max_counts': args.max_counts
        },
        'clustering_resolutions': args.clustering_resolutions
    }
    
    with open(out_dir / 'atac_initial_qc_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Close the AnnDataSet
    data.close()
    
    print(f"Results saved to: {out_dir}")
    print(f"AnnDataSet: {dataset_path}")
    print(f"QC plots: {plot_dir}")
    print(f"Sample statistics: {out_dir / 'sample_statistics.csv'}")
    print(f"\n  NEXT STEPS:")
    print(f"  1. Review QC plots to determine appropriate thresholds")
    print(f"  2. Annotate cell types using clustering results")
    print(f"  3. Run comprehensive pipeline (atac_comprehensive.py) with:")
    print(f"     - Cell type annotations for peak calling")
    print(f"     - Sample-specific thresholds (optional)")

if __name__ == "__main__":
    main()
