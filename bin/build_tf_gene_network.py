#!/usr/bin/env python3
"""
build_tf_gene_network.py — Phase 3

Construct a global TF-gene adjacency matrix from scPrinter binding scores
and Cicero co-accessibility, mirroring Shi et al.'s TF_Net methodology
but using continuous binding scores instead of binary TOBIAS classifications.

Inputs:
    - tf_target_genes.json (from MAP_TF_TO_TARGET_GENES)
    - enhancer_fp_summary.csv files (collected from ENHANCER_FOOTPRINTING)
    - ccan_enhancer_gene_links.tsv (from EXTRACT_CCAN_ENHANCERS)
    - cicero_connections.tsv.gz (for co-accessibility edge weights)
    - GTF (for promoter classification)
    - Binding threshold (otsu, percentile_75, or float)

Outputs:
    - tf_gene_adjacency.tsv (full edge list)
    - tf_gene_adjacency_matrix.tsv (pivot: TF x gene)
    - tf_gene_network.graphml (igraph-compatible)
    - network_summary.txt
"""

import argparse
import json
import os
import re
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Build TF-gene regulatory network")
    p.add_argument("--tf-targets", required=True,
                   help="tf_target_genes.json from MAP_TF_TO_TARGET_GENES")
    p.add_argument("--enhancer-gene-links", required=True,
                   help="ccan_enhancer_gene_links.tsv from EXTRACT_CCAN_ENHANCERS")
    p.add_argument("--cicero-connections", required=True,
                   help="Cicero connections TSV (for co-accessibility edge weights)")
    p.add_argument("--summary-dir", required=True,
                   help="Directory containing enhancer_fp_summary.csv files")
    p.add_argument("--gtf", required=True, help="GTF annotation for promoter classification")
    p.add_argument("--binding-threshold", default="otsu",
                   help="'otsu', 'percentile_75', or float value")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def parse_gtf_promoter_genes(gtf_path, upstream=2000, downstream=500):
    """Extract promoter regions and their gene names from GTF."""
    promoter_map = {}
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
            gene_name = attrs.get('gene_name')
            if not gene_name:
                continue
            chrom = parts[0]
            strand = parts[6]
            start = int(parts[3])
            end = int(parts[4])
            tss = start if strand == '+' else end
            prom_start = max(0, tss - upstream)
            prom_end = tss + downstream
            promoter_map[gene_name] = {
                'chr': chrom, 'tss': tss,
                'prom_start': prom_start, 'prom_end': prom_end,
            }
    return promoter_map


def compute_threshold(scores, method='otsu'):
    """Compute bound/unbound threshold."""
    scores = scores[~np.isnan(scores)]
    if len(scores) == 0:
        return 0.5
    try:
        return float(method)
    except (ValueError, TypeError):
        pass
    if method == 'otsu':
        try:
            from skimage.filters import threshold_otsu
            return float(threshold_otsu(scores))
        except (ImportError, ValueError):
            return float(np.percentile(scores, 75))
    return float(np.percentile(scores, 75))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("BUILD TF-GENE NETWORK (Phase 3)")
    print("=" * 60)

    # Load TF-target gene mapping
    print(f"\nLoading TF-target genes from {args.tf_targets}")
    with open(args.tf_targets) as f:
        tf_data = json.load(f)

    # Load enhancer-gene links
    print(f"Loading enhancer-gene links from {args.enhancer_gene_links}")
    enh_links = pd.read_csv(args.enhancer_gene_links, sep='\t')
    print(f"  {len(enh_links)} enhancer-gene links")

    # Load Cicero connections for co-accessibility weights
    print(f"Loading Cicero connections from {args.cicero_connections}")
    conns = pd.read_csv(args.cicero_connections, sep='\t')
    coaccess_col = conns.columns[2] if len(conns.columns) >= 3 else 'coaccess'
    conns[coaccess_col] = pd.to_numeric(conns[coaccess_col], errors='coerce')

    # Build peak -> max coaccess lookup
    peak_coaccess = {}
    for _, row in conns.iterrows():
        for peak in [str(row[conns.columns[0]]), str(row[conns.columns[1]])]:
            score = row[coaccess_col]
            if pd.notna(score):
                peak_coaccess[peak] = max(peak_coaccess.get(peak, 0), score)

    # Load promoter map from GTF
    print(f"Loading GTF promoter annotations from {args.gtf}")
    promoter_map = parse_gtf_promoter_genes(args.gtf)
    print(f"  {len(promoter_map)} genes with promoter annotations")

    # Collect all summary files to understand which (cell_type, TF) pairs exist
    summary_dir = Path(args.summary_dir)
    summary_files = list(summary_dir.glob("**/*.csv"))
    print(f"\nFound {len(summary_files)} enhancer footprint summary files")

    ct_tf_pairs = []
    for sf in summary_files:
        try:
            df = pd.read_csv(sf)
            for _, row in df.iterrows():
                ct_tf_pairs.append((row['cell_type'], row['tf']))
        except Exception:
            continue

    print(f"  {len(ct_tf_pairs)} (cell_type, TF) pairs")

    # Build edge list
    edges = []

    cell_types_data = tf_data.get('cell_types', {})

    for ct, tf in ct_tf_pairs:
        ct_data = cell_types_data.get(ct, {})
        tf_info = ct_data.get(tf, ct_data.get(tf.upper(), {}))

        target_genes = tf_info.get('target_genes', [])
        evidence = tf_info.get('evidence', {})

        for gene in target_genes:
            # Determine edge type: promoter (direct) or enhancer (indirect)
            is_promoter_direct = False

            # Check if TF has motif peaks in the gene's promoter
            if gene in promoter_map:
                is_promoter_direct = True

            # Look up co-accessibility weight from enhancer-gene links
            gene_links = enh_links[enh_links['gene'] == gene] if not enh_links.empty else pd.DataFrame()
            max_coaccess = 0.0

            if not gene_links.empty:
                for _, link in gene_links.iterrows():
                    peak_str = f"{link['chr']}_{link['start']}_{link['end']}"
                    coaccess = peak_coaccess.get(peak_str, 0.0)
                    max_coaccess = max(max_coaccess, coaccess)

            # Evidence count from tf_target_genes.json
            evidence_count = evidence.get(gene, 0)

            edge_type = 'promoter_direct' if is_promoter_direct else 'enhancer_indirect'
            coaccess_weight = 1.0 if is_promoter_direct else max_coaccess

            edges.append({
                'TF': tf,
                'target_gene': gene,
                'cell_type': ct,
                'edge_type': edge_type,
                'evidence_count': evidence_count,
                'coaccessibility': coaccess_weight,
            })

    edges_df = pd.DataFrame(edges)
    print(f"\nTotal edges: {len(edges_df)}")
    if not edges_df.empty:
        print(f"  Promoter direct: {(edges_df['edge_type'] == 'promoter_direct').sum()}")
        print(f"  Enhancer indirect: {(edges_df['edge_type'] == 'enhancer_indirect').sum()}")
        print(f"  Unique TFs: {edges_df['TF'].nunique()}")
        print(f"  Unique target genes: {edges_df['target_gene'].nunique()}")
        print(f"  Unique cell types: {edges_df['cell_type'].nunique()}")

    # Write full edge list
    adj_path = outdir / 'tf_gene_adjacency.tsv'
    edges_df.to_csv(adj_path, sep='\t', index=False)
    print(f"\nWrote edge list: {adj_path}")

    # Write pivot adjacency matrix (TF x gene, values = max evidence count)
    if not edges_df.empty:
        pivot = edges_df.pivot_table(
            index='TF', columns='target_gene',
            values='evidence_count', aggfunc='max', fill_value=0)
        matrix_path = outdir / 'tf_gene_adjacency_matrix.tsv'
        pivot.to_csv(matrix_path, sep='\t')
        print(f"Wrote adjacency matrix: {matrix_path} ({pivot.shape[0]} TFs x {pivot.shape[1]} genes)")
    else:
        pd.DataFrame().to_csv(outdir / 'tf_gene_adjacency_matrix.tsv', sep='\t')

    # Write GraphML for network analysis tools
    try:
        import networkx as nx

        G = nx.DiGraph()
        for _, row in edges_df.iterrows():
            G.add_edge(
                row['TF'], row['target_gene'],
                edge_type=row['edge_type'],
                cell_type=row['cell_type'],
                evidence=int(row['evidence_count']),
                coaccessibility=float(row['coaccessibility']))

        graphml_path = outdir / 'tf_gene_network.graphml'
        nx.write_graphml(G, str(graphml_path))
        print(f"Wrote GraphML network: {graphml_path} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
    except ImportError:
        print("  [WARN] networkx not available — skipping GraphML output")

    # Summary
    summary_lines = [
        "TF-Gene Network Construction Summary",
        "=" * 45,
        f"Total edges: {len(edges_df)}",
    ]
    if not edges_df.empty:
        summary_lines.extend([
            f"Promoter direct edges: {(edges_df['edge_type'] == 'promoter_direct').sum()}",
            f"Enhancer indirect edges: {(edges_df['edge_type'] == 'enhancer_indirect').sum()}",
            f"Unique TFs: {edges_df['TF'].nunique()}",
            f"Unique target genes: {edges_df['target_gene'].nunique()}",
            f"Unique cell types: {edges_df['cell_type'].nunique()}",
            f"",
            f"Per cell type:",
        ])
        for ct_name in sorted(edges_df['cell_type'].unique()):
            ct_edges = edges_df[edges_df['cell_type'] == ct_name]
            summary_lines.append(
                f"  {ct_name}: {len(ct_edges)} edges, "
                f"{ct_edges['TF'].nunique()} TFs, "
                f"{ct_edges['target_gene'].nunique()} targets")

    summary_path = outdir / 'network_summary.txt'
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
