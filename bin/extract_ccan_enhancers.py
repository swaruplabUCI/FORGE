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
    p.add_argument("--ccan", required=False, default=None, help="CCAN assignments file (required for 'ccan' mode)")
    p.add_argument("--gtf", required=True, help="GTF file for TSS extraction")
    p.add_argument("--promoter-upstream", type=int, default=2000)
    p.add_argument("--promoter-downstream", type=int, default=500)
    p.add_argument("--outdir", default=".")
    p.add_argument("--mode", choices=['ccan', 'pairwise_95'], default='ccan',
                   help="Enhancer filtering mode: 'ccan' (CCAN membership) or "
                        "'pairwise_95' (Shi et al. 95th percentile threshold)")
    p.add_argument("--percentile-threshold", type=float, default=95,
                   help="Percentile threshold for pairwise_95 mode (default: 95)")
    p.add_argument("--condition-label", default=None,
                   help="Optional condition label for output file naming (e.g., '5xFAD')")
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


def run_pairwise_95_mode(args, tss_df, tss_by_chr, outdir):
    """Shi et al. pairwise 95th percentile co-accessibility threshold mode.

    Mirrors PiD_AD_putativeEnhancer.Rmd:
      1. Load pairwise Cicero connections
      2. Classify Peak1/Peak2 as promoter or non-promoter
      3. Keep only promoter-to-non-promoter pairs with coaccess > 0
      4. Apply 95th percentile threshold
      5. Surviving non-promoter peaks = high-confidence enhancers
    """
    percentile = args.percentile_threshold
    cond_label = args.condition_label or ''
    suffix = f"_{cond_label}" if cond_label else ''

    print(f"\n{'='*60}")
    print(f"PAIRWISE_95 MODE — {percentile}th percentile threshold")
    if cond_label:
        print(f"  Condition: {cond_label}")
    print(f"{'='*60}")

    # Load Cicero connections (pairwise)
    print(f"Loading Cicero connections from {args.connections}")
    conns = pd.read_csv(args.connections, sep='\t')
    print(f"  Total connections: {len(conns)}")

    # Detect column names for Peak1, Peak2, coaccess
    peak1_col = peak2_col = coaccess_col = None
    for col in conns.columns:
        cl = col.lower()
        if 'peak1' in cl or col == conns.columns[0]:
            if peak1_col is None:
                peak1_col = col
        if 'peak2' in cl or col == conns.columns[1]:
            if peak2_col is None:
                peak2_col = col
        if 'coaccess' in cl or 'score' in cl or col == conns.columns[2]:
            if coaccess_col is None:
                coaccess_col = col

    # Ensure we have all three
    if peak1_col is None:
        peak1_col = conns.columns[0]
    if peak2_col is None:
        peak2_col = conns.columns[1]
    if coaccess_col is None:
        coaccess_col = conns.columns[2]

    print(f"  Columns: peak1='{peak1_col}', peak2='{peak2_col}', coaccess='{coaccess_col}'")

    # Convert coaccess to numeric, drop NaN
    conns[coaccess_col] = pd.to_numeric(conns[coaccess_col], errors='coerce')
    conns = conns.dropna(subset=[coaccess_col])

    # Step 1: Filter to coaccess > 0
    conns_pos = conns[conns[coaccess_col] > 0].copy()
    print(f"  Connections with coaccess > 0: {len(conns_pos)}")

    # Step 2: Classify Peak1 and Peak2 as promoter or non-promoter
    def classify_peak(peak_str):
        chrom, start, end = parse_peak_coords(str(peak_str))
        if chrom is None:
            return None, None, None, False, []
        chr_tss = tss_by_chr.get(chrom, pd.DataFrame())
        if chr_tss.empty:
            return chrom, start, end, False, []
        is_prom, genes = overlaps_promoter(chrom, start, end, chr_tss)
        return chrom, start, end, is_prom, genes

    print("  Classifying peaks as promoter/non-promoter...")
    p1_class = conns_pos[peak1_col].map(lambda x: classify_peak(x))
    p2_class = conns_pos[peak2_col].map(lambda x: classify_peak(x))

    conns_pos['p1_is_promoter'] = p1_class.map(lambda x: x[3])
    conns_pos['p1_genes'] = p1_class.map(lambda x: x[4])
    conns_pos['p2_is_promoter'] = p2_class.map(lambda x: x[3])
    conns_pos['p2_genes'] = p2_class.map(lambda x: x[4])
    conns_pos['p2_chr'] = p2_class.map(lambda x: x[0])
    conns_pos['p2_start'] = p2_class.map(lambda x: x[1])
    conns_pos['p2_end'] = p2_class.map(lambda x: x[2])
    conns_pos['p1_chr'] = p1_class.map(lambda x: x[0])
    conns_pos['p1_start'] = p1_class.map(lambda x: x[1])
    conns_pos['p1_end'] = p1_class.map(lambda x: x[2])

    # Drop rows where either peak couldn't be parsed
    conns_pos = conns_pos.dropna(subset=['p1_chr', 'p2_chr'])

    # Step 3: Filter to promoter-to-non-promoter pairs (either direction)
    mask_p1_prom_p2_enh = conns_pos['p1_is_promoter'] & ~conns_pos['p2_is_promoter']
    mask_p1_enh_p2_prom = ~conns_pos['p1_is_promoter'] & conns_pos['p2_is_promoter']

    # Normalize: always put promoter as "anchor" and enhancer as "target"
    prom_enh_pairs = []

    # Case 1: Peak1=promoter, Peak2=enhancer
    subset1 = conns_pos[mask_p1_prom_p2_enh]
    for _, row in subset1.iterrows():
        prom_enh_pairs.append({
            'enh_chr': row['p2_chr'], 'enh_start': int(row['p2_start']),
            'enh_end': int(row['p2_end']), 'enh_peak': str(row[peak2_col]),
            'prom_genes': row['p1_genes'],
            'coaccess': row[coaccess_col],
        })

    # Case 2: Peak1=enhancer, Peak2=promoter
    subset2 = conns_pos[mask_p1_enh_p2_prom]
    for _, row in subset2.iterrows():
        prom_enh_pairs.append({
            'enh_chr': row['p1_chr'], 'enh_start': int(row['p1_start']),
            'enh_end': int(row['p1_end']), 'enh_peak': str(row[peak1_col]),
            'prom_genes': row['p2_genes'],
            'coaccess': row[coaccess_col],
        })

    pairs_df = pd.DataFrame(prom_enh_pairs)
    print(f"  Promoter-to-enhancer pairs: {len(pairs_df)}")

    if pairs_df.empty:
        print("WARNING: No promoter-to-enhancer pairs found")
        _write_empty_outputs(outdir, suffix)
        return

    # Step 4: Apply percentile threshold
    threshold = np.percentile(pairs_df['coaccess'].values, percentile)
    print(f"  {percentile}th percentile co-accessibility threshold: {threshold:.4f}")

    pairs_thresh = pairs_df[pairs_df['coaccess'] >= threshold].copy()
    print(f"  Pairs above threshold: {len(pairs_thresh)}")

    # Step 5: Deduplicate enhancers — keep the best coaccess per enhancer peak
    # Each enhancer is linked to its promoter anchor's gene(s)
    enhancer_records = []
    seen_peaks = {}

    for _, row in pairs_thresh.iterrows():
        peak_key = f"{row['enh_chr']}:{row['enh_start']}-{row['enh_end']}"
        genes = row['prom_genes'] if isinstance(row['prom_genes'], list) else []
        gene_str = ','.join(sorted(set(genes))) if genes else ''

        if peak_key in seen_peaks:
            # Merge genes from multiple promoter anchors
            existing_genes = set(seen_peaks[peak_key]['linked_genes'].split(',')) if seen_peaks[peak_key]['linked_genes'] else set()
            existing_genes.update(genes)
            seen_peaks[peak_key]['linked_genes'] = ','.join(sorted(existing_genes - {''}))
            seen_peaks[peak_key]['max_coaccess'] = max(seen_peaks[peak_key]['max_coaccess'], row['coaccess'])
        else:
            seen_peaks[peak_key] = {
                'chr': row['enh_chr'], 'start': row['enh_start'], 'end': row['enh_end'],
                'CCAN_id': f"pw95_{len(seen_peaks)}",
                'linked_genes': gene_str,
                'max_coaccess': row['coaccess'],
            }

    enhancer_df = pd.DataFrame(list(seen_peaks.values()))
    print(f"\n  Unique enhancer peaks after dedup: {len(enhancer_df)}")

    # Write outputs using shared helper
    _write_outputs(enhancer_df, outdir, suffix, len(conns), len(conns_pos),
                   len(pairs_df), threshold, percentile)


def _write_empty_outputs(outdir, suffix=''):
    """Write empty output files."""
    import gzip as gz
    bed_path = outdir / f'ccan_enhancer_peaks{suffix}.bed.gz'
    with gz.open(bed_path, 'wt') as f:
        pass
    links_path = outdir / f'ccan_enhancer_gene_links{suffix}.tsv'
    pd.DataFrame(columns=['chr', 'start', 'end', 'CCAN_id', 'gene']).to_csv(
        links_path, sep='\t', index=False)
    summary_path = outdir / f'ccan_enhancer_summary{suffix}.txt'
    summary_path.write_text("No enhancer peaks found\n")


def _write_outputs(enhancer_df, outdir, suffix, n_total_conns, n_pos_conns,
                   n_prom_enh_pairs, threshold, percentile):
    """Write BED, gene links TSV, and summary for enhancer peaks."""
    # Write enhancer BED
    if not enhancer_df.empty:
        bed_df = enhancer_df[['chr', 'start', 'end', 'CCAN_id', 'linked_genes']].copy()
        bed_df = bed_df.sort_values(['chr', 'start']).reset_index(drop=True)

        bed_path = outdir / f'ccan_enhancer_peaks{suffix}.bed.gz'
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False, compression='gzip')
        print(f"\nWrote {len(bed_df)} enhancer peaks to {bed_path}")
    else:
        _write_empty_outputs(outdir, suffix)
        return

    # Write gene links TSV
    links_records = []
    for _, row in enhancer_df.iterrows():
        if row['linked_genes']:
            for gene in row['linked_genes'].split(','):
                if gene.strip():
                    links_records.append({
                        'chr': row['chr'], 'start': row['start'],
                        'end': row['end'], 'CCAN_id': row['CCAN_id'],
                        'gene': gene.strip(),
                    })

    links_df = pd.DataFrame(links_records) if links_records else pd.DataFrame(
        columns=['chr', 'start', 'end', 'CCAN_id', 'gene'])
    links_path = outdir / f'ccan_enhancer_gene_links{suffix}.tsv'
    links_df.to_csv(links_path, sep='\t', index=False)
    print(f"Wrote {len(links_df)} enhancer-gene links to {links_path}")

    # Summary
    summary_lines = [
        "CCAN Enhancer Extraction Summary",
        "=" * 40,
        f"Mode: pairwise_95 ({percentile}th percentile)",
        f"Total Cicero connections: {n_total_conns}",
        f"Connections with coaccess > 0: {n_pos_conns}",
        f"Promoter-to-enhancer pairs: {n_prom_enh_pairs}",
        f"Co-accessibility threshold: {threshold:.4f}",
        f"Enhancer peaks above threshold: {len(enhancer_df)}",
        f"Enhancer-gene links: {len(links_df)}",
        f"Unique genes linked: {links_df['gene'].nunique() if not links_df.empty else 0}",
    ]
    if suffix:
        summary_lines.insert(2, f"Condition: {suffix.lstrip('_')}")

    summary_path = outdir / f'ccan_enhancer_summary{suffix}.txt'
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse TSS
    tss_df = parse_gtf_tss(args.gtf, args.promoter_upstream, args.promoter_downstream)

    # Build a per-chromosome index for faster lookup
    tss_by_chr = {chrom: grp for chrom, grp in tss_df.groupby('chr')}

    # Dispatch based on mode
    if args.mode == 'pairwise_95':
        run_pairwise_95_mode(args, tss_df, tss_by_chr, outdir)
        return

    # ---- Original CCAN mode ----
    # Mirror pairwise_95's per-condition filename suffix so CTRL+TRT outputs
    # don't collide when staged together (e.g. BUILD_UNION_ENHANCER_PEAKS).
    cond_label = args.condition_label or ''
    suffix = f"_{cond_label}" if cond_label else ''

    if args.ccan is None:
        print("ERROR: --ccan is required for 'ccan' mode", file=sys.stderr)
        sys.exit(1)

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

        bed_path = outdir / f'ccan_enhancer_peaks{suffix}.bed.gz'
        bed_df.to_csv(bed_path, sep='\t', index=False, header=False, compression='gzip')
        print(f"\nWrote {len(bed_df)} enhancer peaks to {bed_path}")
    else:
        # Write empty BED
        bed_path = outdir / f'ccan_enhancer_peaks{suffix}.bed.gz'
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
    links_path = outdir / f'ccan_enhancer_gene_links{suffix}.tsv'
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
    summary_path = outdir / f'ccan_enhancer_summary{suffix}.txt'
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
