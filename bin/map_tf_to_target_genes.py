#!/usr/bin/env python3
"""
map_tf_to_target_genes.py

Map ChromVAR-enriched TFs to their putative target genes using the
motif-peak match matrix (varm["motif_match"]) and Cicero CCANs.

Logic:
  1. Load chromvar_matrix.h5ad which contains:
     - var_names: peak coordinates (e.g., chr1:10000-10500)
     - varm["motif_match"]: peaks x motifs boolean matrix
     - uns["motif_name"]: motif identifiers (e.g., "MA0489.2\tJun")
  2. Load per_celltype_motifs.json (from EXTRACT_CHROMVAR_MOTIFS):
     - Per cell type: enriched TF names + motif IDs
  3. Load Cicero CCAN assignments (peak → CCAN_id)
  4. Load GTF for TSS → promoter window mapping
  5. For each (cell_type, TF):
     a. Find motif index in uns["motif_name"]
     b. Get peaks where motif_match[:,motif_idx] is True
     c. Map those peaks to CCANs
     d. Find which genes' promoters are in the same CCANs
     e. Those genes are putative TF targets

Outputs:
  - tf_target_genes.json: per-cell-type mapping {ct: {tf: [gene1, gene2, ...]}}
    Plus flat lists for RESOLVE_GENE_COORDINATES and SCPRINTER_FOOTPRINTING
  - tf_target_mapping_report.txt: human-readable summary
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Map TFs to target genes via motif-peak-CCAN linkage")
    p.add_argument("--chromvar-matrix", required=True,
                   help="chromvar_matrix.h5ad with varm['motif_match'] and uns['motif_name']")
    p.add_argument("--motif-json", required=True,
                   help="per_celltype_motifs.json from EXTRACT_CHROMVAR_MOTIFS")
    p.add_argument("--ccan", required=True,
                   help="CCAN assignments TSV (peak → CCAN_id)")
    p.add_argument("--gtf", required=True,
                   help="GTF annotation file for TSS extraction")
    p.add_argument("--promoter-upstream", type=int, default=2000)
    p.add_argument("--promoter-downstream", type=int, default=500)
    p.add_argument("--max-targets-per-tf", type=int, default=20,
                   help="Maximum target genes per TF to keep (ranked by # motif peaks in CCAN)")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def parse_peak_coords(peak_str):
    """Parse peak string like 'chr1_12345_12500' or 'chr1:12345-12500'."""
    if '_' in peak_str and not peak_str.startswith('chr_'):
        parts = peak_str.split('_')
        if len(parts) >= 3:
            try:
                return parts[0], int(parts[1]), int(parts[2])
            except ValueError:
                pass
    if ':' in peak_str and '-' in peak_str:
        chrom, rest = peak_str.split(':', 1)
        start, end = rest.split('-', 1)
        return chrom, int(start), int(end)
    return None, None, None


def parse_gtf_tss(gtf_path, upstream=2000, downstream=500):
    """Extract TSS positions and promoter windows from GTF."""
    print(f"Parsing TSS from GTF: {gtf_path}")
    records = []
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
            records.append({
                'chr': chrom, 'tss': tss,
                'prom_start': prom_start, 'prom_end': prom_end,
                'gene_name': gene_name, 'strand': strand,
            })
    tss_df = pd.DataFrame(records)
    print(f"  Extracted {len(tss_df)} gene TSS positions")
    return tss_df


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Load chromvar matrix with motif_match ----
    print(f"Loading chromvar matrix from {args.chromvar_matrix}")
    adata = ad.read_h5ad(args.chromvar_matrix)
    print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} peaks")

    if 'motif_match' not in adata.varm:
        print("ERROR: chromvar_matrix.h5ad missing varm['motif_match']. "
              "Was chromvar_scan() run?", file=sys.stderr)
        sys.exit(1)

    motif_match = adata.varm['motif_match']  # peaks x motifs
    motif_names_raw = list(adata.uns.get('motif_name', []))
    peak_names = list(adata.var_names)

    print(f"  motif_match: {motif_match.shape} (peaks x motifs)")
    print(f"  motif_names: {len(motif_names_raw)} entries")
    print(f"  Example motif names: {motif_names_raw[:3]}")

    # Parse motif names: "MA0489.2\tJun" → (motif_id, tf_name)
    motif_index = {}  # tf_name → list of motif column indices
    motif_id_to_idx = {}  # motif_id → column index
    for i, raw_name in enumerate(motif_names_raw):
        parts = str(raw_name).split('\t')
        motif_id = parts[0]
        tf_name = parts[-1] if len(parts) >= 2 else parts[0]
        motif_id_to_idx[motif_id] = i
        if tf_name not in motif_index:
            motif_index[tf_name] = []
        motif_index[tf_name].append(i)

    # ---- 2. Load per-cell-type motif JSON ----
    print(f"\nLoading motif JSON from {args.motif_json}")
    with open(args.motif_json) as f:
        motif_data = json.load(f)

    # ---- 3. Load CCAN assignments ----
    print(f"\nLoading CCAN assignments from {args.ccan}")
    ccan_df = pd.read_csv(args.ccan, sep='\t')
    peak_col = ccan_df.columns[0]
    ccan_col = ccan_df.columns[1]
    print(f"  {len(ccan_df)} peak-CCAN assignments, {ccan_df[ccan_col].nunique()} CCANs")

    # Build peak → CCAN lookup (normalize peak format)
    peak_to_ccan = {}
    for _, row in ccan_df.iterrows():
        peak_str = str(row[peak_col])
        peak_to_ccan[peak_str] = row[ccan_col]

    # Also try normalized format (chr:start-end ↔ chr_start_end)
    for _, row in ccan_df.iterrows():
        peak_str = str(row[peak_col])
        # Add both formats
        chrom, start, end = parse_peak_coords(peak_str)
        if chrom:
            peak_to_ccan[f"{chrom}:{start}-{end}"] = row[ccan_col]
            peak_to_ccan[f"{chrom}_{start}_{end}"] = row[ccan_col]

    # ---- 4. Parse GTF for TSS/promoter windows ----
    tss_df = parse_gtf_tss(args.gtf, args.promoter_upstream, args.promoter_downstream)
    tss_by_chr = {chrom: grp for chrom, grp in tss_df.groupby('chr')}

    # ---- 5. Build CCAN → gene mapping (which CCANs contain gene promoters) ----
    print("\nMapping CCANs to genes via promoter overlap...")
    ccan_to_genes = defaultdict(set)

    for peak_str, ccan_id in peak_to_ccan.items():
        chrom, start, end = parse_peak_coords(peak_str)
        if chrom is None:
            continue
        chr_tss = tss_by_chr.get(chrom, pd.DataFrame())
        if chr_tss.empty:
            continue
        # Check if this peak overlaps any promoter
        mask = (
            (chr_tss['prom_start'] < end) &
            (chr_tss['prom_end'] > start)
        )
        if mask.any():
            genes = chr_tss.loc[mask, 'gene_name'].unique()
            for gene in genes:
                ccan_to_genes[ccan_id].add(gene)

    print(f"  {len(ccan_to_genes)} CCANs contain gene promoters")
    total_genes = len(set().union(*ccan_to_genes.values())) if ccan_to_genes else 0
    print(f"  {total_genes} unique genes with CCAN promoter anchors")

    # ---- 6. For each (cell_type, TF): find target genes ----
    print("\nMapping TFs to target genes...")
    result = {
        'cell_types': {},
        'all_target_genes': [],
        'per_celltype_targets': {},  # {ct: [genes]} for SCPRINTER_FOOTPRINTING fan-out
    }
    all_target_genes = set()
    report_lines = [
        "TF → Target Gene Mapping Report",
        "=" * 65,
        f"Motif-peak matrix: {motif_match.shape}",
        f"CCANs with promoter anchors: {len(ccan_to_genes)}",
        f"Max targets per TF: {args.max_targets_per_tf}",
        "",
    ]

    for ct_name, ct_data in motif_data.get('cell_types', {}).items():
        top_tfs = ct_data.get('top_tfs', [])
        motif_details = {d['tf']: d for d in ct_data.get('motif_details', [])}

        ct_result = {}
        ct_all_targets = set()
        report_lines.append(f"--- {ct_name} ({ct_data.get('n_cells', '?')} cells) ---")

        for tf in top_tfs:
            # Find motif column indices for this TF
            if tf not in motif_index:
                report_lines.append(f"  {tf}: no motif match in matrix (skipped)")
                continue

            motif_idxs = motif_index[tf]

            # Get peaks where any of this TF's motifs match
            tf_peak_mask = np.zeros(motif_match.shape[0], dtype=bool)
            for idx in motif_idxs:
                col = motif_match[:, idx]
                if hasattr(col, 'toarray'):
                    col = col.toarray().ravel()
                elif hasattr(col, 'A'):
                    col = col.A.ravel()
                else:
                    col = np.asarray(col).ravel()
                tf_peak_mask |= (col > 0)

            n_motif_peaks = int(tf_peak_mask.sum())

            if n_motif_peaks == 0:
                report_lines.append(f"  {tf}: 0 peaks with motif (skipped)")
                continue

            # Map motif-containing peaks to CCANs, then to genes
            gene_evidence = defaultdict(int)  # gene → count of motif peaks in its CCAN

            for peak_idx in np.where(tf_peak_mask)[0]:
                peak_name = peak_names[peak_idx]
                ccan_id = peak_to_ccan.get(peak_name)
                if ccan_id is None:
                    # Try alternative format
                    chrom, start, end = parse_peak_coords(peak_name)
                    if chrom:
                        ccan_id = peak_to_ccan.get(f"{chrom}:{start}-{end}")
                        if ccan_id is None:
                            ccan_id = peak_to_ccan.get(f"{chrom}_{start}_{end}")
                if ccan_id is None:
                    continue
                for gene in ccan_to_genes.get(ccan_id, set()):
                    gene_evidence[gene] += 1

            if not gene_evidence:
                report_lines.append(
                    f"  {tf}: {n_motif_peaks} motif peaks, but none in CCANs with promoters")
                # Still include TF itself as fallback
                ct_result[tf] = {
                    'target_genes': [tf],
                    'n_motif_peaks': n_motif_peaks,
                    'evidence': {tf: 0},
                    'note': 'fallback_to_self',
                }
                ct_all_targets.add(tf)
                continue

            # Rank genes by evidence (number of motif peaks in their CCAN)
            ranked = sorted(gene_evidence.items(), key=lambda x: -x[1])
            top_targets = [g for g, _ in ranked[:args.max_targets_per_tf]]

            ct_result[tf] = {
                'target_genes': top_targets,
                'n_motif_peaks': n_motif_peaks,
                'evidence': {g: c for g, c in ranked[:args.max_targets_per_tf]},
            }
            ct_all_targets.update(top_targets)

            report_lines.append(
                f"  {tf}: {n_motif_peaks} motif peaks → {len(gene_evidence)} genes "
                f"(top {len(top_targets)}: {', '.join(top_targets[:5])}{'...' if len(top_targets) > 5 else ''})"
            )

        result['cell_types'][ct_name] = ct_result
        result['per_celltype_targets'][ct_name] = sorted(ct_all_targets)
        all_target_genes.update(ct_all_targets)
        report_lines.append(f"  Total unique targets for {ct_name}: {len(ct_all_targets)}")
        report_lines.append("")

    result['all_target_genes'] = sorted(all_target_genes)

    # ---- 7. Write outputs ----
    json_path = outdir / 'tf_target_genes.json'
    with open(json_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nWrote TF→target mapping to {json_path}")
    print(f"  Total unique target genes: {len(all_target_genes)}")

    report_lines.append(f"\nTotal unique target genes across all cell types: {len(all_target_genes)}")
    report_path = outdir / 'tf_target_mapping_report.txt'
    report_path.write_text('\n'.join(report_lines))
    print(f"Wrote report to {report_path}")


if __name__ == '__main__':
    main()
