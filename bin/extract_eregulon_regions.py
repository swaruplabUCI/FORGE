#!/usr/bin/env python3
"""
extract_eregulon_regions.py — Recipe Steps B3/C4

Parse SCENIC+ eRegulon TSVs into per-TF enhancer region sets and target gene
lists.  Optionally intersect with DORC significant peaks.

Inputs:
    - eRegulon_direct.tsv (from SCENIC+)
    - region_to_gene_adj.tsv (from SCENIC+)
    - dorc_significant.tsv.gz (optional, from scPrinter DORC)

Outputs:
    - eregulon_region_sets/ (per-TF BED files)
    - eregulon_target_genes.json ({TF: [genes]})
    - eregulon_manifest.json
    - eregulon_dorc_intersection/ (optional)
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Extract eRegulon regions for footprinting")
    p.add_argument("--ereg-direct", required=True, help="eRegulon_direct.tsv from SCENIC+")
    p.add_argument("--r2g", required=True, help="region_to_gene_adj.tsv from SCENIC+")
    p.add_argument("--dorc-sig", default=None, help="DORC significant peaks TSV.gz (optional)")
    p.add_argument("--min-regions", type=int, default=10,
                   help="Minimum regions per TF to write a region set")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def parse_region(region_str):
    """Parse region string to (chr, start, end). Handles 'chr1:100-200' and 'chr1_100_200'."""
    if ':' in region_str and '-' in region_str:
        chrom, rest = region_str.split(':', 1)
        start, end = rest.split('-', 1)
        return chrom, int(start), int(end)
    elif '_' in region_str:
        parts = region_str.split('_')
        if len(parts) >= 3:
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
                return chrom, start, end
            except ValueError:
                pass
    return None, None, None


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    region_sets_dir = outdir / "eregulon_region_sets"
    region_sets_dir.mkdir(parents=True, exist_ok=True)

    # Load eRegulon data
    print(f"Loading eRegulon data from {args.ereg_direct}")
    ereg_df = pd.read_csv(args.ereg_direct, sep='\t')
    print(f"  Loaded {len(ereg_df)} eRegulon entries")

    # SCENIC+ V2 emits eRegulon_direct.tsv with one row per (TF, Region, Gene)
    # triplet — columns: Region, Gene, TF, importance_R2G, rho_R2G, ... .
    # The (TF, Region) pairs are already filtered through SCENIC+'s own
    # quality checks (motif enrichment + region-to-gene importance). We use
    # ereg_df directly when it carries Region+Gene+TF, and only fall back to
    # the gene→region join through r2g_df when ereg_df lacks Region.
    def _resolve(df, predicates, fallback):
        """Return first column matching one of the predicates, else fallback."""
        cols = list(df.columns)
        for pred in predicates:
            for c in cols:
                if pred(c):
                    return c
        return fallback

    tf_col = _resolve(
        ereg_df,
        [lambda c: c.lower() == 'tf',
         lambda c: 'tf' in c.lower() and 'tf2g' not in c.lower()],
        ereg_df.columns[0],
    )
    gene_col = _resolve(
        ereg_df,
        [lambda c: c.lower() == 'gene',
         lambda c: c.lower() == 'target',
         lambda c: 'target' in c.lower() and 'gene' in c.lower(),
         lambda c: 'gene' in c.lower() and 'name' not in c.lower() and 'signature' not in c.lower()],
        ereg_df.columns[-1],
    )
    region_col = _resolve(
        ereg_df,
        [lambda c: c.lower() == 'region',
         lambda c: 'region' in c.lower() and 'name' not in c.lower() and 'signature' not in c.lower()],
        None,
    )

    has_region_in_ereg = region_col is not None and region_col in ereg_df.columns
    print(f"  Detected columns: TF='{tf_col}', gene='{gene_col}', region='{region_col}'")

    # tf_genes is built from ereg_df regardless of source.
    tf_genes = {}
    for _, row in ereg_df.iterrows():
        tf = str(row[tf_col]).strip()
        gene = str(row[gene_col]).strip()
        if tf and gene and tf != 'nan' and gene != 'nan':
            tf_genes.setdefault(tf, set()).add(gene)
    print(f"  Found {len(tf_genes)} TFs with target genes")

    # tf_regions: prefer ereg_df's Region column directly (SCENIC+ V2 path).
    # When unavailable, gracefully fall back to the gene→region join via r2g_df.
    tf_regions = {}
    if has_region_in_ereg:
        for _, row in ereg_df.iterrows():
            tf = str(row[tf_col]).strip()
            region = str(row[region_col]).strip()
            if tf and region and tf != 'nan' and region != 'nan':
                tf_regions.setdefault(tf, set()).add(region)
        print(f"  Built tf→regions directly from ereg_df.{region_col}")
    else:
        print(f"Loading region-to-gene adjacency from {args.r2g}")
        r2g_df = pd.read_csv(args.r2g, sep='\t')
        print(f"  Loaded {len(r2g_df)} region-to-gene entries")

        r2g_region_col = _resolve(
            r2g_df,
            [lambda c: c.lower() == 'region',
             lambda c: 'region' in c.lower()],
            r2g_df.columns[0],
        )
        r2g_gene_col = _resolve(
            r2g_df,
            [lambda c: c.lower() == 'gene',
             lambda c: c.lower() == 'target',
             lambda c: 'target' in c.lower() and 'gene' in c.lower(),
             lambda c: 'gene' in c.lower()],
            None,
        )
        if r2g_gene_col is None or r2g_gene_col == r2g_region_col:
            raise SystemExit(
                f"r2g column auto-detect failed: region={r2g_region_col!r} "
                f"gene={r2g_gene_col!r}, columns={list(r2g_df.columns)}. "
                "Cannot derive tf→regions without an unambiguous gene column."
            )
        print(f"  r2g columns: region='{r2g_region_col}', gene='{r2g_gene_col}'")

        gene_regions = {}
        for _, row in r2g_df.iterrows():
            region = str(row[r2g_region_col]).strip()
            gene = str(row[r2g_gene_col]).strip()
            if region and gene and region != 'nan' and gene != 'nan':
                gene_regions.setdefault(gene, set()).add(region)

        for tf, genes in tf_genes.items():
            regions = set()
            for gene in genes:
                if gene in gene_regions:
                    regions.update(gene_regions[gene])
            if regions:
                tf_regions[tf] = regions

    print(f"  TFs with linked regions: {len(tf_regions)}")

    # Optional DORC intersection
    dorc_peaks = None
    if args.dorc_sig:
        print(f"Loading DORC significant peaks from {args.dorc_sig}")
        try:
            dorc_df = pd.read_csv(args.dorc_sig, sep='\t', compression='gzip')
            dorc_intersection_dir = outdir / "eregulon_dorc_intersection"
            dorc_intersection_dir.mkdir(exist_ok=True)

            # Extract peak coordinates from DORC
            dorc_peaks = set()
            for col in dorc_df.columns:
                if 'peak' in col.lower() or 'region' in col.lower():
                    dorc_peaks.update(dorc_df[col].dropna().astype(str).tolist())
                    break
            if not dorc_peaks:
                # Try first column
                dorc_peaks = set(dorc_df.iloc[:, 0].dropna().astype(str).tolist())

            print(f"  DORC significant peaks: {len(dorc_peaks)}")
        except Exception as e:
            print(f"  WARNING: Failed to load DORC data: {e}")
            dorc_peaks = None

    # Write per-TF region sets
    manifest_entries = []
    target_genes_dict = {}

    for tf in sorted(tf_regions.keys()):
        regions = tf_regions[tf]
        target_genes_dict[tf] = sorted(tf_genes.get(tf, set()))

        # Parse regions to BED format
        bed_records = []
        for region_str in regions:
            chrom, start, end = parse_region(region_str)
            if chrom is not None:
                bed_records.append({
                    'Chromosome': chrom,
                    'Start': start,
                    'End': end,
                })

        if len(bed_records) < args.min_regions:
            continue

        bed_df = pd.DataFrame(bed_records).sort_values(['Chromosome', 'Start'])
        tf_safe = tf.replace(' ', '_').replace('/', '_')
        bed_fname = f"{tf_safe}.bed"
        bed_path = region_sets_dir / bed_fname
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False)

        manifest_entries.append({
            'tf': tf,
            'bed_file': bed_fname,
            'n_regions': len(bed_df),
            'n_target_genes': len(tf_genes.get(tf, set())),
            'source': 'eregulon',
        })

        # DORC intersection if available
        if dorc_peaks is not None:
            dorc_hits = [r for r in regions if r in dorc_peaks]
            if dorc_hits:
                dorc_bed = []
                for r in dorc_hits:
                    c, s, e = parse_region(r)
                    if c:
                        dorc_bed.append({'Chromosome': c, 'Start': s, 'End': e})
                if dorc_bed:
                    dorc_bed_df = pd.DataFrame(dorc_bed).sort_values(['Chromosome', 'Start'])
                    dorc_bed_path = dorc_intersection_dir / f"{tf_safe}_dorc.bed"
                    dorc_bed_df.to_csv(dorc_bed_path, sep='\t', index=False, header=False)

    # Write manifest
    manifest = {
        'region_sets': manifest_entries,
        'total_tfs': len(tf_regions),
        'tfs_with_regions': len(manifest_entries),
        'min_regions': args.min_regions,
    }
    manifest_path = outdir / "eregulon_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    # Write target genes
    tg_path = outdir / "eregulon_target_genes.json"
    with open(tg_path, 'w') as f:
        json.dump(target_genes_dict, f, indent=2)

    # Summary
    summary_lines = [
        "eRegulon Region Extraction Summary",
        "=" * 40,
        f"Total eRegulon entries: {len(ereg_df)}",
        f"TFs found: {len(tf_genes)}",
        f"TFs with linked regions: {len(tf_regions)}",
        f"Region sets written (>= {args.min_regions} regions): {len(manifest_entries)}",
        f"Total target genes: {sum(len(v) for v in target_genes_dict.values())}",
    ]
    if dorc_peaks is not None:
        summary_lines.append(f"DORC peaks intersected: {len(dorc_peaks)}")

    summary_path = outdir / "eregulon_summary.txt"
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
