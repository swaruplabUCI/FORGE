#!/usr/bin/env python3
"""
extract_ccan_enhancers.py — Recipe Step A4

Extract distal enhancer peaks from Cicero CCANs.  Peaks within a CCAN that
do NOT overlap a TSS (+/- promoter window) are labelled as putative enhancers,
linked to genes via the CCAN's promoter anchor(s).

Inputs:
    - Cicero connections TSV (co-accessibility scores)
    - CCAN assignment file (peak → CCAN_id mapping)
    - GTF annotation file (gene → TSS)

Outputs:
    - ccan_enhancer_peaks.bed.gz  (chr, start, end, CCAN_id, linked_genes)
    - ccan_enhancer_gene_links.tsv
    - ccan_enhancer_summary.txt
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Extract CCAN enhancer peaks")
    p.add_argument("--connections", required=True, help="Cicero connections TSV")
    p.add_argument("--ccan", required=True, help="CCAN assignments file")
    p.add_argument("--gtf", required=True, help="GTF file for TSS extraction")
    p.add_argument("--promoter-upstream", type=int, default=2000)
    p.add_argument("--promoter-downstream", type=int, default=500)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def parse_gtf_tss(gtf_path, upstream=2000, downstream=500):
    """Extract TSS positions and promoter windows from GTF."""
    print(f"Parsing TSS from GTF: {gtf_path}")
    tss_records = []

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue

            attrs = {}
            for match in re.finditer(r'(\w+) "([^"]+)"', parts[8]):
                attrs[match.group(1)] = match.group(2)

            gene_name = attrs.get('gene_name', attrs.get('gene_symbol'))
            if not gene_name:
                continue

            chrom = parts[0]
            strand = parts[6]
            start = int(parts[3])
            end = int(parts[4])

            tss = start if strand == '+' else end

            prom_start = max(0, tss - upstream) if strand == '+' else max(0, tss - downstream)
            prom_end = tss + downstream if strand == '+' else tss + upstream

            tss_records.append({
                'chr': chrom,
                'tss': tss,
                'prom_start': prom_start,
                'prom_end': prom_end,
                'gene_name': gene_name,
                'strand': strand,
            })

    tss_df = pd.DataFrame(tss_records)
    print(f"  Extracted {len(tss_df)} gene TSS positions")
    return tss_df


def parse_peak_coords(peak_str):
    """Parse peak string like 'chr1_12345_12500' or 'chr1:12345-12500' to (chr, start, end)."""
    if '_' in peak_str and not peak_str.startswith('chr_'):
        parts = peak_str.split('_')
        if len(parts) >= 3:
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
                return chrom, start, end
            except ValueError:
                pass
    if ':' in peak_str and '-' in peak_str:
        chrom, rest = peak_str.split(':', 1)
        start, end = rest.split('-', 1)
        return chrom, int(start), int(end)
    return None, None, None


def overlaps_promoter(chrom, start, end, tss_df):
    """Check if a peak overlaps any TSS promoter region."""
    mask = (
        (tss_df['chr'] == chrom) &
        (tss_df['prom_start'] < end) &
        (tss_df['prom_end'] > start)
    )
    if mask.any():
        return True, tss_df.loc[mask, 'gene_name'].unique().tolist()
    return False, []


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse TSS
    tss_df = parse_gtf_tss(args.gtf, args.promoter_upstream, args.promoter_downstream)

    # Build a per-chromosome index for faster lookup
    tss_by_chr = {chrom: grp for chrom, grp in tss_df.groupby('chr')}

    # Load CCAN assignments
    print(f"Loading CCAN assignments from {args.ccan}")
    ccan_df = pd.read_csv(args.ccan, sep='\t')

    # Detect column names
    peak_col = None
    ccan_col = None
    for col in ccan_df.columns:
        if 'peak' in col.lower() or 'site' in col.lower():
            peak_col = col
        if 'ccan' in col.lower() or 'module' in col.lower() or 'cluster' in col.lower():
            ccan_col = col

    if peak_col is None or ccan_col is None:
        # Fallback: assume first two columns
        peak_col = ccan_df.columns[0]
        ccan_col = ccan_df.columns[1]

    print(f"  Using columns: peak='{peak_col}', ccan='{ccan_col}'")
    print(f"  Total peaks in CCANs: {len(ccan_df)}")
    print(f"  Unique CCANs: {ccan_df[ccan_col].nunique()}")

    # Classify each peak as promoter or enhancer
    enhancer_records = []
    promoter_records = []
    skipped = 0

    for _, row in ccan_df.iterrows():
        peak_str = str(row[peak_col])
        ccan_id = row[ccan_col]
        chrom, start, end = parse_peak_coords(peak_str)

        if chrom is None:
            skipped += 1
            continue

        chr_tss = tss_by_chr.get(chrom, pd.DataFrame())
        if chr_tss.empty:
            # No TSS on this chromosome — treat as enhancer with no gene link
            enhancer_records.append({
                'chr': chrom, 'start': start, 'end': end,
                'CCAN_id': ccan_id, 'peak': peak_str,
                'is_promoter': False, 'linked_genes': '',
            })
            continue

        is_prom, genes = overlaps_promoter(chrom, start, end, chr_tss)

        record = {
            'chr': chrom, 'start': start, 'end': end,
            'CCAN_id': ccan_id, 'peak': peak_str,
            'is_promoter': is_prom,
            'linked_genes': ','.join(genes) if genes else '',
        }

        if is_prom:
            promoter_records.append(record)
        else:
            enhancer_records.append(record)

    print(f"\n  Promoter peaks: {len(promoter_records)}")
    print(f"  Enhancer peaks: {len(enhancer_records)}")
    print(f"  Skipped (unparseable): {skipped}")

    # Build gene linkage: for each CCAN, find which genes its promoter anchors link to
    promoter_df = pd.DataFrame(promoter_records)
    enhancer_df = pd.DataFrame(enhancer_records)

    if not promoter_df.empty:
        ccan_gene_map = {}
        for _, row in promoter_df.iterrows():
            ccan_id = row['CCAN_id']
            genes = row['linked_genes'].split(',') if row['linked_genes'] else []
            if ccan_id not in ccan_gene_map:
                ccan_gene_map[ccan_id] = set()
            ccan_gene_map[ccan_id].update(genes)

        # Assign linked genes to enhancer peaks via CCAN membership
        if not enhancer_df.empty:
            enhancer_df['linked_genes'] = enhancer_df['CCAN_id'].map(
                lambda x: ','.join(sorted(ccan_gene_map.get(x, set())))
            )

    # Write enhancer BED
    if not enhancer_df.empty:
        bed_df = enhancer_df[['chr', 'start', 'end', 'CCAN_id', 'linked_genes']].copy()
        bed_df = bed_df.sort_values(['chr', 'start']).reset_index(drop=True)

        bed_path = outdir / 'ccan_enhancer_peaks.bed.gz'
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False, compression='gzip')
        print(f"\nWrote {len(bed_df)} enhancer peaks to {bed_path}")
    else:
        # Write empty BED
        bed_path = outdir / 'ccan_enhancer_peaks.bed.gz'
        with gzip.open(bed_path, 'wt') as f:
            pass
        print("WARNING: No enhancer peaks found")

    # Write gene links TSV
    links_records = []
    if not enhancer_df.empty:
        for _, row in enhancer_df.iterrows():
            if row['linked_genes']:
                for gene in row['linked_genes'].split(','):
                    links_records.append({
                        'chr': row['chr'],
                        'start': row['start'],
                        'end': row['end'],
                        'CCAN_id': row['CCAN_id'],
                        'gene': gene.strip(),
                    })

    links_df = pd.DataFrame(links_records) if links_records else pd.DataFrame(
        columns=['chr', 'start', 'end', 'CCAN_id', 'gene']
    )
    links_path = outdir / 'ccan_enhancer_gene_links.tsv'
    links_df.to_csv(links_path, sep='\t', index=False)
    print(f"Wrote {len(links_df)} enhancer-gene links to {links_path}")

    # Summary
    summary_lines = [
        "CCAN Enhancer Extraction Summary",
        "=" * 40,
        f"Total peaks in CCANs: {len(ccan_df)}",
        f"Promoter peaks: {len(promoter_records)}",
        f"Enhancer peaks: {len(enhancer_records)}",
        f"Unique CCANs: {ccan_df[ccan_col].nunique()}",
        f"CCANs with promoter anchors: {len(set(r['CCAN_id'] for r in promoter_records))}",
        f"Enhancer-gene links: {len(links_df)}",
        f"Unique genes linked: {links_df['gene'].nunique() if not links_df.empty else 0}",
    ]
    summary_path = outdir / 'ccan_enhancer_summary.txt'
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
