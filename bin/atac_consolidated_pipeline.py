#!/usr/bin/env python3
"""
Initial ATAC-seq QC with UNIFORM filtering on list of AnnData objects
Corrected version based on successful test_simple.py workflow
"""
import os
import sys

# Fix LD_LIBRARY_PATH for snapatac2
# os.environ['LD_LIBRARY_PATH'] = '/dfs7/swaruplab/lesolano/tools/pysnATAC_nf/lib:' + os.environ.get('LD_LIBRARY_PATH', '')
import snapatac2 as snap
import numpy as np
import pandas as pd
import json
import argparse
import logging
import contextlib
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import re
import gc
import shutil
import scanpy as sc

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
    parser.add_argument('--thresholds_file', type=str, default=None,
                   help='JSON file with per-sample thresholds')
    # FIX-29b: Local GTF annotation to avoid internet download on HPC compute nodes
    parser.add_argument('--gtf', type=str, default=None,
                    help='Local GTF/GFF3 annotation file (avoids internet download)')
    parser.add_argument('--cisbp_meme', type=str, default=None,
                    help='Local CIS-BP MEME motif file for motif enrichment (avoids internet download)')
    # FIX-44: Configurable genome build (was hardcoded mm10 for mouse)
    parser.add_argument('--genome_build', default=None,
                    help='Genome build (e.g. hg38, mm10, mm39). Auto-detected from species if not set.')
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

    # FIX-29b: Handle FIX-28 manifest where sample_id already includes _{batch} suffix.
    # Avoid double-suffix (e.g. L1_Donor_0_july_july) by checking first.
    def make_qc_sample(row):
        sid = str(row['sample_id'])
        batch = str(row['batch'])
        if sid.endswith('_' + batch):
            return sid
        return sid + '_' + batch
    metadata['qc_sample'] = metadata.apply(make_qc_sample, axis=1)

    sample_names = metadata['qc_sample'].tolist()

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

    # Use the unique QC sample names for output file names
    out_files = [out_dir / f"{sample}.h5ad" for sample in sample_names]
    # Final combined AnnDataSet filename
    dataset_path = out_dir / "atac_complete.h5ads"

    # FIX-44: Define genome — use explicit build if provided, else default per species
    genome_map = {'hg38': snap.genome.hg38, 'hg19': snap.genome.hg19,
                  'mm10': snap.genome.mm10, 'mm39': snap.genome.mm39,
                  'GRCm39': snap.genome.mm39, 'GRCh38': snap.genome.hg38}
    if args.genome_build and args.genome_build in genome_map:
        genome = genome_map[args.genome_build]
    else:
        genome = snap.genome.hg38 if args.species == 'human' else snap.genome.mm10

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
    # FIX-29b: Use local GTF if provided (HPC compute nodes have no internet)
    tsse_anno = args.gtf if args.gtf else genome
    if args.gtf:
        print(f"  Using local annotation: {args.gtf}")
    snap.metrics.tsse(GMWM_all, tsse_anno, n_jobs=args.n_jobs)

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

    # ========================================================================
    # STEP 4: APPLY SAMPLE-SPECIFIC FILTERING
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 4: Applying sample-specific filtering")
    print(f"{'='*70}")

    # Load thresholds from JSON config
    if args.thresholds_file and os.path.exists(args.thresholds_file):
        with open(args.thresholds_file, 'r') as f:
            thr_map = json.load(f)
        print(f" Loaded sample-specific thresholds from {args.thresholds_file}")
        print(f"   Samples in config: {list(thr_map.keys())}")
    else:
        print(" No thresholds file provided - using uniform thresholds")
        thr_map = {name: {
            "min_counts": args.min_counts,
            "min_tsse": args.min_tsse,
            "max_counts": args.max_counts
        } for name in sample_names}

    # Close any open backed file handles from import
    try:
        for obj in GMWM_all:
            if hasattr(obj, "flush"):
                obj.flush()
            if hasattr(obj, "file") and obj.file is not None:
                obj.file.close()
    except:
        pass

    import gc
    gc.collect()

    # Sequential per-sample processing (load -> filter -> save)
    filtered_files = []
    for sname, out_path in zip(sample_names, out_files):
        print(f"\n{'='*50}")
        print(f"Processing {sname}...")
        print(f"{'='*50}")

        # Load into memory (NOT backed mode)
        ad_mem = snap.read(str(out_path), backed=None)

        # Get thresholds for this sample
        thr = thr_map.get(sname, {
            "min_counts": args.min_counts,
            "min_tsse": args.min_tsse,
            "max_counts": args.max_counts
        })

        n_before = ad_mem.n_obs
        print(f"  Cells before filtering: {n_before}")
        print(f"  Thresholds:")
        print(f"    min_counts: {thr.get('min_counts', args.min_counts)}")
        print(f"    min_tsse: {thr.get('min_tsse', args.min_tsse)}")
        print(f"    max_counts: {thr.get('max_counts', args.max_counts)}")

        # Filter cells (returns boolean index)
        keep_idx = snap.pp.filter_cells(
            ad_mem,
            min_counts=thr.get("min_counts", args.min_counts),
            min_tsse=thr.get("min_tsse", args.min_tsse),
            max_counts=thr.get("max_counts", args.max_counts),
            inplace=False
        )

        # Subset and copy
        ad_mem = ad_mem[keep_idx, :].copy()
        n_after = ad_mem.n_obs
        pct_retained = (n_after / n_before * 100) if n_before > 0 else 0

        print(f"  Cells after filtering: {n_after} ({pct_retained:.1f}% retained)")

        # Skip downstream steps if no cells remain
        if n_after == 0:
            print(f"  WARNING: No cells retained for {sname}! Skipping tile matrix/scrublet.")
            # Write empty file to maintain consistency
            ad_mem.write_h5ad(str(out_path), compression='gzip')
            filtered_files.append(out_path)
            del ad_mem
            gc.collect()
            continue

        # Add tile matrix
        snap.pp.add_tile_matrix(ad_mem, bin_size=5000)
        print(f"  Tile matrix added")

        # Feature selection for scrublet
        snap.pp.select_features(ad_mem, n_features=args.n_features)
        print(f"  Feature selection complete")

        # Scrublet doublet detection
        feat_mask = None
        if "selected" in ad_mem.var.columns:
            feat_mask = np.asarray(ad_mem.var["selected"], dtype=bool)

        # Guard against very small libraries that break Scrublet
        if ad_mem.n_obs < 50:
            print(f"  WARNING: {sname} has only {ad_mem.n_obs} cells; skipping Scrublet/doublet filtering.")
        else:
            try:
                snap.pp.scrublet(ad_mem, features=feat_mask)
                print(f"  Scrublet scoring complete")

                # Keep a backup before doublet filtering in case it removes everything
                ad_mem_backup = ad_mem.copy()
                n_before_dbl = ad_mem.n_obs
                snap.pp.filter_doublets(ad_mem, probability_threshold=0.5, inplace=True)
                # Guard: if ALL cells were removed, scrublet scores are unreliable --
                # skip doublet filtering entirely for this sample.
                if ad_mem.n_obs == 0:
                    print(f"  WARNING: Doublet filter removed ALL {n_before_dbl} cells -- "
                          f"scrublet scores unreliable; skipping doublet removal for {sname}")
                    ad_mem = ad_mem_backup
                del ad_mem_backup
                print(f"  After doublet removal: {ad_mem.n_obs} cells")
            except Exception as e:
                print(f"  WARNING: Scrublet failed for {sname} ({e}); skipping doublet filtering for this sample.")

        # Overwrite file atomically
        out_tmp = str(Path(out_path).with_suffix(".tmp.h5ad"))
        ad_mem.write_h5ad(out_tmp, compression='gzip')
        os.replace(out_tmp, out_path)
        print(f" Saved to {out_path}")

        filtered_files.append(out_path)

        del ad_mem
        gc.collect()

    print(f"\n{'='*70}")
    print(f"Sample-specific filtering complete")
    print(f"{'='*70}")

    # Validate all h5ad files before proceeding
    print("\nValidating written h5ad files...")
    valid_files = []
    corrupted_files = []

    for f in filtered_files:
        try:
            # Test read with backed mode
            test_adata = snap.read(str(f), backed='r')
            n_cells = test_adata.n_obs
            test_adata.close()  # Close file handle

            print(f"  {f.name}: {n_cells} cells")
            valid_files.append(f)

        except Exception as e:
            print(f"  {f.name}: CORRUPTED - {str(e)}")
            corrupted_files.append(f)

    if corrupted_files:
        print(f"\n  WARNING: {len(corrupted_files)} corrupted files detected:")
        for f in corrupted_files:
            print(f"    - {f.name}")

        print("\n  Proceeding with valid files only. Corrupted files will be skipped.")

    # Update filtered_files to only include valid files
    filtered_files = valid_files

    if not filtered_files:
        raise RuntimeError("No valid h5ad files found after filtering!")

    # FIX-45a: Rename per-sample h5ads from qc_sample → sample_id BEFORE
    # creating AnnDataSet, so that obs_names, sample columns, peak_matrix,
    # and gene_matrix all use sample_id consistently (matching RNA barcodes).
    sid_map = dict(zip(metadata['qc_sample'], metadata['sample_id']))
    renamed_files = []
    for f in filtered_files:
        qc_name = f.stem
        if qc_name in sid_map and qc_name != sid_map[qc_name]:
            dst = f.parent / f"{sid_map[qc_name]}.h5ad"
            if not dst.exists():
                os.rename(str(f), str(dst))
                print(f"  Renamed {f.name} -> {dst.name} (qc_sample → sample_id)")
                renamed_files.append(dst)
            else:
                renamed_files.append(f)
        else:
            renamed_files.append(f)
    filtered_files = renamed_files

    print(f"\n Proceeding with {len(filtered_files)} valid samples")
    # --------------------------------------------------------------------
    # Export per-sample QC h5ad files for scPrinter / downstream use
    # --------------------------------------------------------------------
    scp_qc_dir = out_dir / "scprinter_qc_h5ads"
    scp_qc_dir.mkdir(exist_ok=True)

    scp_qc_files = []
    for f in filtered_files:
        stem = f.stem  # e.g. "sample1_batch1"
        dest = scp_qc_dir / f"{stem}.qc.h5ad"
        if dest.exists():
            dest.unlink()
        # Use copy to be safe on DFS/NFS
        import shutil
        shutil.copy2(f, dest)
        scp_qc_files.append(dest)

    print("\nExported per-sample QC AnnData files for scPrinter:")
    for f in scp_qc_files:
        print(f"  {f}")

    # Safe cell counting with explicit file closing
    total_cells = 0
    for f in filtered_files:
        try:
            adata_temp = snap.read(str(f), backed='r')
            total_cells += adata_temp.n_obs
            if hasattr(adata_temp, 'close'):
                adata_temp.close()
        except Exception as e:
            print(f"  Error reading {f.name}: {e}")
            continue

    print(f"   Total cells retained: {total_cells}")
    print(f"{'='*70}")

    # IMPORTANT: Reload filtered files for downstream steps
    # Tile matrices and features were already added in Step 4
    print(f"\n{'='*70}")
    print("STEP 5: Reloading filtered samples (read-only)")
    print(f"{'='*70}")

    GMWM_all = [snap.read(str(f), backed='r') for f in filtered_files]
    sample_names_final = [Path(f).stem for f in filtered_files]
    out_files_final = filtered_files

    # ========================================================================
    # STEP 9: CREATE ANNDATASET FROM SAVED FILES
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 9: Creating AnnDataSet")
    print(f"{'='*70}")

    data = snap.AnnDataSet(
        adatas=[(name, str(path)) for name, path in zip(sample_names_final, out_files_final)],
        filename=str(dataset_path),
        add_key='sample'
    )

    #CRITICAL: Create unique cell IDs and verify (from lab code Page 1)
    # snapatac2 AnnDataSet.obs is a PyDataFrameElem (polars-backed), not pandas.
    # Check for 'sample' column using try/except to handle both APIs.
    has_sample = False
    try:
        _ = data.obs['sample']
        has_sample = True
    except (KeyError, TypeError):
        pass

    if not has_sample:
        # Build sample labels manually from the constituent AnnData objects
        sample_labels = []
        for name, path in zip(sample_names_final, out_files_final):
            adata_tmp = snap.read(str(path), backed='r')
            sample_labels.extend([name] * adata_tmp.n_obs)
            if hasattr(adata_tmp, 'close'):
                adata_tmp.close()
        data.obs = pd.DataFrame({'sample': sample_labels}, index=data.obs_names)
        print(f"  Built .obs manually ({len(sample_labels)} cells, {len(sample_names_final)} samples)")

    unique_cell_ids = [sa + ':' + bc for sa, bc in zip(data.obs['sample'], data.obs_names)]
    data.obs_names = unique_cell_ids

    # CRITICAL: Assert uniqueness (from lab code Page 1)
    assert data.n_obs == np.unique(data.obs_names).size, "Cell IDs are not unique!"

    print(f"AnnDataSet created: {data.n_obs} cells")
    print(f"Verified unique cell IDs: {np.unique(data.obs_names).size}")

    # FIX-45a: After early rename, sample names in the AnnDataSet are sample_id,
    # so index by sample_id for downstream lookups.
    meta_map = metadata.set_index('sample_id')

    if 'batch' not in meta_map.columns:
        raise KeyError("'batch' column not found in metadata; cannot annotate AnnDataSet with batch.")

    # Map the batch label onto every cell via its sample ID.
    # Instead, convert the metadata Series to a plain dict and look up per sample.
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

    snap.tl.spectral(data)  #  Uses AnnDataSet

    # ========================================================================
    # STEP 12: BATCH CORRECTION (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print(f"STEP 12: Batch correction ({args.batch_correction})")
    print(f"{'='*70}")

    if args.batch_correction == 'none':
        print("  Skipping batch correction (single sample / batch_correction='none')")
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
    # STEP 15: PEAK CALLING (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 15: Calling peaks")
    print(f"{'='*70}")

    # Use last resolution
    clustering_key = f"leiden_{str(args.clustering_resolutions[-1]).replace('.', '_')}"

    snap.tl.macs3(
        data,
        groupby=clustering_key,
        qvalue=args.peak_fdr,
        extsize=150,
        shift=-75
    )

    # Filter to standard chromosomes
    if args.species == 'human':
        standard_chroms = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY'}
    else:
        standard_chroms = {f'chr{i}' for i in range(1, 20)} | {'chrX', 'chrY'}

    filtered_peaks = {}
    for group, peaks_df in data.uns['macs3'].items():
        if isinstance(peaks_df, pd.DataFrame):
            mask = peaks_df['chrom'].isin(standard_chroms)
            filtered = peaks_df[mask].copy()
        else:  # Polars
            import polars as pl
            filtered = peaks_df.filter(pl.col('chrom').is_in(list(standard_chroms)))
        filtered_peaks[group] = filtered
        print(f"  {group}: {len(filtered)} peaks")

    # Merge peaks
    union_peaks = snap.tl.merge_peaks(filtered_peaks, chrom_sizes=genome)

    if hasattr(union_peaks, 'to_pandas'):
        peak_list = union_peaks['Peaks'].to_list()
    elif isinstance(union_peaks, pd.DataFrame):
        peak_list = union_peaks['Peaks'].tolist()
    else:
        peak_list = union_peaks

    print(f"Union peak set: {len(peak_list)} peaks")

    # ========================================================================
    # STEP 16: CREATE PEAK MATRIX (on AnnDataSet)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 16: Creating peak matrix")
    print(f"{'='*70}")

    peak_matrix = snap.pp.make_peak_matrix(
        data,
        use_rep=peak_list,
        counting_strategy='fragment',
        inplace=False
    )

    # Transfer metadata
    try:
        obs_cols = list(data.obs.keys())
    except Exception:
        obs_cols = []

    for col in obs_cols:
        if col not in peak_matrix.obs.columns:
            peak_matrix.obs[col] = data.obs[col].to_list()

    # Transfer UMAP coordinates from AnnDataSet to peak_matrix
    try:
        umap_coords = np.array(data.obsm['X_umap'])
        peak_matrix.obsm['X_umap'] = umap_coords
        print(f"  Transferred UMAP coordinates to peak_matrix: {umap_coords.shape}")
    except (KeyError, Exception) as e:
        print(f"  Warning: Could not transfer UMAP to peak_matrix: {e}")

    peak_matrix.write_h5ad(out_dir / "peak_matrix.h5ad", compression='gzip')
    print(f"Peak matrix saved: {peak_matrix.n_obs} cells x {peak_matrix.n_vars} peaks")

    # ========================================================================
    # STEP 16.5: COMPUTE GENE ACTIVITY SCORES (USE ORIGINAL DATA!)
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 16.5: Computing gene activity scores")
    print(f"{'='*70}")

    # FIX-29b: Use local GTF if provided (avoids internet download on HPC)
    # FIX-44: Use genome_map for correct build (was hardcoded mm10)
    if args.gtf:
        gene_annot = args.gtf
        print(f"  Using local annotation for gene matrix: {args.gtf}")
    else:
        gene_annot = genome  # Uses FIX-44 genome_map result from above

    # Compute gene activity matrix from the ORIGINAL data object (has fragments)
    gene_matrix = snap.pp.make_gene_matrix(
        data,  # AnnDataSet with fragments
        gene_anno=gene_annot,
        upstream=5000,
        downstream=0,
        inplace=False
    )

    print(f"Gene matrix created: {gene_matrix.n_obs} cells x {gene_matrix.n_vars} genes")

    # FIX: Materialize AnnDataSet.obs to Polars DataFrame first
    print("  Materializing AnnDataSet.obs...")
    try:
        obs_materialized = data.obs[:]  # Returns Polars DataFrame
        print(f"  Type: {type(obs_materialized)}")
        print(f"  Columns: {list(obs_materialized.columns)}")

        # Convert Polars to Pandas for compatibility with AnnData
        if hasattr(obs_materialized, 'to_pandas'):
            obs_df = obs_materialized.to_pandas()
        else:
            # Already pandas
            obs_df = obs_materialized

        print(f"  Converted to pandas: {type(obs_df)}")

        # Transfer metadata columns
        for col in obs_df.columns:
            if col not in gene_matrix.obs.columns:
                gene_matrix.obs[col] = obs_df[col].values
                print(f"    Transferred: {col}")

        print(f"  Gene matrix now has {len(gene_matrix.obs.columns)} metadata columns")

    except Exception as e:
        print(f"  Warning: Could not transfer from AnnDataSet: {e}")
        print("  Falling back to peak_matrix metadata...")

    # Transfer UMAP coordinates from AnnDataSet to gene_matrix
    try:
        umap_coords = np.array(data.obsm['X_umap'])
        gene_matrix.obsm['X_umap'] = umap_coords
        print(f"  Transferred UMAP coordinates to gene_matrix: {umap_coords.shape}")
    except (KeyError, Exception) as e:
        print(f"  Warning: Could not transfer UMAP from AnnDataSet: {e}")
        # Fallback: use peak_matrix's UMAP if available (transferred earlier)
        if 'X_umap' in peak_matrix.obsm:
            gene_matrix.obsm['X_umap'] = peak_matrix.obsm['X_umap']
            print(f"  Fallback: transferred UMAP from peak_matrix: {peak_matrix.obsm['X_umap'].shape}")
        else:
            print(f"  Warning: No UMAP coordinates available for gene_matrix")

    # Normalize gene scores (log1p transform)
    gene_matrix.X = gene_matrix.X.astype(np.float32)
    # Save gene matrix BEFORE normalization
    gene_matrix_path = out_dir / "gene_matrix.h5ad"
    #  Use Scanpy for normalization (gene_matrix is an AnnData object)
    print("  Normalizing gene activity scores...")
    # Normalize to 10,000 counts per cell
    sc.pp.normalize_total(gene_matrix, target_sum=10000)

    # Log-transform
    sc.pp.log1p(gene_matrix)

    print(f"  Normalized and log-transformed gene scores")

    # Save gene matrix
    gene_matrix_path = out_dir / "gene_matrix.h5ad"
    gene_matrix.write_h5ad(gene_matrix_path, compression='gzip')
    print(f"Gene matrix saved: {gene_matrix_path}")

    # ========================================================================
    # STEP 16.6: GENERATE MARKER GENE DOTPLOT
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 16.6: Generating marker gene dotplot")
    print(f"{'='*70}")

    # Derive marker genes from data via differential expression
    clustering_key = "leiden_0_5"

    try:
        sc.tl.rank_genes_groups(gene_matrix, clustering_key, method='wilcoxon')

        # Extract top 5 genes per cluster, deduplicate
        top_genes = []
        for group in gene_matrix.obs[clustering_key].unique():
            genes = gene_matrix.uns['rank_genes_groups']['names'][str(group)][:5]
            top_genes.extend(genes)
        top_genes = list(dict.fromkeys(top_genes))  # deduplicate, preserve order

        print(f"Derived {len(top_genes)} DE marker genes across {gene_matrix.obs[clustering_key].nunique()} clusters")

        if top_genes:
            import shutil
            sc.pl.dotplot(
                gene_matrix,
                top_genes,
                groupby=clustering_key,
                standard_scale='var',
                save=f'_{clustering_key}.pdf',
                show=False
            )

            src = Path('./figures') / f'dotplot_{clustering_key}.pdf'
            if src.exists():
                shutil.move(src, plot_dir / f"marker_dotplot_{clustering_key}.pdf")
                print(f"Marker dotplot saved: marker_dotplot_{clustering_key}.pdf")

    except Exception as e:
        print(f"WARNING: Could not generate dotplot: {e}")

    # Export cluster-level aggregated gene scores
    cluster_avg_expr = {}
    for cluster_id in set(gene_matrix.obs[clustering_key]):
        cluster_mask = gene_matrix.obs[clustering_key] == cluster_id
        cluster_cells = gene_matrix[cluster_mask]
        mean_expr = np.array(cluster_cells.X.mean(axis=0)).flatten()
        cluster_avg_expr[cluster_id] = mean_expr

    cluster_expr_df = pd.DataFrame(
        cluster_avg_expr,
        index=gene_matrix.var_names
    ).T

    cluster_expr_path = out_dir / "cluster_avg_gene_scores.csv"
    cluster_expr_df.to_csv(cluster_expr_path)
    print(f"Cluster-averaged gene scores saved: {cluster_expr_path}")

    # ========================================================================
    # STEP 16.7: DATA-DRIVEN CLUSTER ANNOTATION
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 16.7: Data-driven cluster annotation")
    print(f"{'='*70}")

    # Use rank_genes_groups results from STEP 16.6 to label clusters by top DE gene
    clustering_key = "leiden_0_5"
    try:
        top_gene_per_cluster = {}
        for group in gene_matrix.obs[clustering_key].unique():
            top_gene = gene_matrix.uns['rank_genes_groups']['names'][str(group)][0]
            top_gene_per_cluster[str(group)] = top_gene

        # Label clusters as "Cluster_N (TOP_GENE)"
        cluster_labels = {k: f"Cluster_{k} ({v})" for k, v in top_gene_per_cluster.items()}
        peak_matrix.obs['cluster_label'] = peak_matrix.obs[clustering_key].map(cluster_labels)
        data.obs['cluster_label'] = [cluster_labels.get(str(c), f'Cluster_{c}') for c in data.obs[clustering_key]]

        print("\nCluster labels (top DE gene):")
        for k, v in sorted(cluster_labels.items(), key=lambda x: int(x[0])):
            n_cells = (peak_matrix.obs[clustering_key] == k).sum()
            print(f"  {v}: {n_cells} cells")

    except Exception as e:
        print(f"Warning: Could not generate cluster labels: {e}")

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
        color='batch',
        show=False,
        out_file=None,
        interactive=False,
        height=500,
    )
    with quiet_plotly_export():
        fig.write_image(str(plot_dir / "umap_by_batch.pdf"))

    # NOTE: Cell type UMAP and marker regions are generated downstream by
    # MERGE_ANNOTATIONS after CellTypist or marker-based annotation is applied.
    # This script only produces cluster-level and sample/batch UMAPs.

    # ========================================================================
    # STEP 18: SAVE SUMMARY AND CLOSE
    # ========================================================================
    print(f"\n{'='*70}")
    print("STEP 18: Saving summary and closing")
    print(f"{'='*70}")

    summary = {
        'n_samples': len(sample_names_final),
        'total_cells': data.n_obs,
        'total_peaks': len(peak_list),
        'samples': sample_names_final,
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

    with open(out_dir / 'atac_pipeline_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    # Close the AnnDataSet
    data.close()

    print(f"\n{'='*70}")
    print("ATAC CONSOLIDATED PIPELINE COMPLETE")
    print(f"{'='*70}")
    print(f"Results saved to: {out_dir}")
    print(f"AnnDataSet: {dataset_path}")
    print(f"Peak matrix: {out_dir / 'peak_matrix.h5ad'}")
    print(f"QC plots: {plot_dir}")

if __name__ == "__main__":
    main()
