#!/usr/bin/env python3
"""
plot_eregulon_networks.py — GRN visualization with graph-tool

Generates static graph images (PDF + PNG) from SCENIC+ outputs:
  1. Per-cell-type GRN: TF → target gene edges for top eRegulons per cell type
  2. Holistic GRN: top eRegulons across all cell types, TFs colored by specificity
  3. TF-TF co-regulation network: edges between TFs sharing target genes

Inputs:
  --eregulon-tsv      eRegulon_direct.tsv (TF→Region→Gene triplets)
  --aucell-direct     AUCell_direct.h5mu (per-cell eRegulon activity scores)
  --rss-csv           eregulon_rss.csv (precomputed RSS specificity scores)
  --tf-to-gene-adj    tf_to_gene_adj.tsv (TF→gene adjacency weights)
  --cell-type-key     Column name for cell type annotations (default: cell_type)
  --top-n             Number of top eRegulons per cell type (default: 5)
  --outdir            Output directory (default: .)
"""
import argparse
import os
import sys
import warnings

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import graph_tool.all as gt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── Color palette ────────────────────────────────────────────────────────────

def build_color_map(categories):
    """Assign distinct RGBA colors to a list of categories."""
    n = len(categories)
    if n <= 10:
        cmap = plt.get_cmap('tab10')
    elif n <= 20:
        cmap = plt.get_cmap('tab20')
    else:
        cmap = plt.get_cmap('turbo')
    return {cat: cmap(i / max(n - 1, 1)) for i, cat in enumerate(sorted(categories))}


# ── Data loading ─────────────────────────────────────────────────────────────

def load_rss(rss_csv):
    """Load precomputed RSS matrix (regulons × cell types)."""
    return pd.read_csv(rss_csv, index_col=0)


def build_eregulon_edges(eregulon_tsv, regulation_filter='+/+'):
    """Build edge DataFrame from eRegulon triplets."""
    df = pd.read_csv(eregulon_tsv, sep='\t')
    if regulation_filter:
        df = df[df['eRegulon_name'].str.contains(
            regulation_filter.replace('+', r'\+'), na=False)]
    edges = df[['TF', 'Gene', 'Region', 'importance_x_rho',
                'rho_R2G', 'eRegulon_name']].copy()
    edges.rename(columns={'importance_x_rho': 'weight',
                          'rho_R2G': 'correlation'}, inplace=True)
    return edges


def get_top_tfs(rss_df, cell_type=None, top_n=5):
    """Get top TFs from RSS matrix for one or all cell types."""
    if cell_type and cell_type in rss_df.columns:
        top_regs = rss_df[cell_type].nlargest(top_n).index.tolist()
        return set(r.split('_')[0] for r in top_regs)
    all_top = set()
    for ct in rss_df.columns:
        top_regs = rss_df[ct].nlargest(top_n).index.tolist()
        all_top.update(r.split('_')[0] for r in top_regs)
    return all_top


def get_best_celltype(tf, rss_df):
    """Find the cell type where a TF's eRegulon has highest RSS."""
    matching = [r for r in rss_df.index if r.split('_')[0] == tf]
    if matching:
        best_reg = matching[0]
        return rss_df.loc[best_reg].idxmax(), rss_df.loc[best_reg].max()
    return 'unknown', 0.0


def aggregate_edges(edges, top_tfs, max_targets_per_tf=None):
    """Filter edges to top TFs, aggregate TF→Gene, optionally limit targets."""
    sub = edges[edges['TF'].isin(top_tfs)]
    if sub.empty:
        return None
    agg = sub.groupby(['TF', 'Gene']).agg(
        weight=('weight', 'sum'),
        correlation=('correlation', 'mean'),
        n_regions=('Region', 'nunique')
    ).reset_index()

    if max_targets_per_tf:
        # Keep only top N targets per TF by weight
        agg = (agg.groupby('TF', group_keys=False)
               .apply(lambda g: g.nlargest(max_targets_per_tf, 'weight'))
               .reset_index(drop=True))
    return agg


# ── Graph building ───────────────────────────────────────────────────────────

def build_gt_tf_gene_graph(edge_agg, rss_df, cell_type_colors):
    """Build a graph-tool directed graph from TF→Gene edges.

    Returns (graph, pos_layout, property_maps_dict).
    """
    g = gt.Graph(directed=True)

    # Vertex properties
    vlabel = g.new_vp('string')
    vtype = g.new_vp('string')       # 'TF' or 'target'
    vcolor = g.new_vp('vector<double>')
    vsize = g.new_vp('double')
    vfontsize = g.new_vp('double')
    vshape = g.new_vp('string')
    vtextcolor = g.new_vp('vector<double>')
    vhalo = g.new_vp('bool')
    vhalocolor = g.new_vp('vector<double>')

    # Edge properties
    eweight = g.new_ep('double')
    ecolor = g.new_ep('vector<double>')
    ewidth = g.new_ep('double')

    node_map = {}  # name → vertex

    def get_or_add(name):
        if name not in node_map:
            v = g.add_vertex()
            node_map[name] = v
            vlabel[v] = name
        return node_map[name]

    # Add TF nodes first
    tf_nodes = edge_agg['TF'].unique()
    for tf in tf_nodes:
        v = get_or_add(tf)
        vtype[v] = 'TF'
        best_ct, best_rss = get_best_celltype(tf, rss_df)
        ct_color = cell_type_colors.get(best_ct, (0.5, 0.5, 0.5, 1.0))
        vcolor[v] = list(ct_color)
        n_targets = edge_agg[edge_agg['TF'] == tf]['Gene'].nunique()
        vsize[v] = max(20, min(70, 20 + n_targets * 0.4))
        vfontsize[v] = max(9, min(16, 9 + n_targets * 0.12))
        vshape[v] = 'circle'
        vtextcolor[v] = [0.0, 0.0, 0.0, 1.0]
        vhalo[v] = True
        vhalocolor[v] = list(ct_color[:3]) + [0.3]

    # Add target gene nodes and edges
    for _, row in edge_agg.iterrows():
        tf_v = node_map[row['TF']]
        gene_v = get_or_add(row['Gene'])
        if vtype[gene_v] != 'TF':
            vtype[gene_v] = 'target'
            vcolor[gene_v] = [0.65, 0.75, 0.88, 0.9]
            vsize[gene_v] = 8
            vfontsize[gene_v] = 0  # no label for targets in large graphs
            vshape[gene_v] = 'circle'
            vtextcolor[gene_v] = [0.2, 0.2, 0.2, 1.0]
            vhalo[gene_v] = False
            vhalocolor[gene_v] = [0.0, 0.0, 0.0, 0.0]

        e = g.add_edge(tf_v, gene_v)
        eweight[e] = row['weight']
        tf_color = list(vcolor[tf_v])
        ecolor[e] = tf_color[:3] + [0.4]
        ewidth[e] = max(0.5, min(3.0, row['weight'] * 0.8))

    return g, {
        'vlabel': vlabel, 'vtype': vtype, 'vcolor': vcolor, 'vsize': vsize,
        'vfontsize': vfontsize, 'vshape': vshape, 'vtextcolor': vtextcolor,
        'vhalo': vhalo, 'vhalocolor': vhalocolor,
        'eweight': eweight, 'ecolor': ecolor, 'ewidth': ewidth
    }


def build_gt_coregulation_graph(edges, rss_df, cell_type_colors, min_shared=5):
    """Build TF-TF co-regulation graph from shared target genes."""
    tf_targets = edges.groupby('TF')['Gene'].apply(set).to_dict()
    tfs = list(tf_targets.keys())

    coedges = []
    for i in range(len(tfs)):
        for j in range(i + 1, len(tfs)):
            shared = tf_targets[tfs[i]] & tf_targets[tfs[j]]
            if len(shared) >= min_shared:
                union = tf_targets[tfs[i]] | tf_targets[tfs[j]]
                jaccard = len(shared) / len(union)
                coedges.append({
                    'source': tfs[i], 'target': tfs[j],
                    'shared': len(shared), 'jaccard': jaccard
                })

    if not coedges:
        return None, None

    coedge_df = pd.DataFrame(coedges)
    g = gt.Graph(directed=False)

    vlabel = g.new_vp('string')
    vcolor = g.new_vp('vector<double>')
    vsize = g.new_vp('double')
    vfontsize = g.new_vp('double')

    eweight = g.new_ep('double')
    ecolor = g.new_ep('vector<double>')
    ewidth = g.new_ep('double')
    elabel = g.new_ep('string')

    node_map = {}
    all_tfs_in_net = set(coedge_df['source']) | set(coedge_df['target'])
    for tf in sorted(all_tfs_in_net):
        v = g.add_vertex()
        node_map[tf] = v
        vlabel[v] = tf
        best_ct, _ = get_best_celltype(tf, rss_df)
        ct_color = cell_type_colors.get(best_ct, (0.5, 0.5, 0.5, 1.0))
        vcolor[v] = list(ct_color)
        n_tgt = len(tf_targets.get(tf, []))
        vsize[v] = max(20, min(70, 20 + n_tgt * 0.05))
        vfontsize[v] = 10

    for _, row in coedge_df.iterrows():
        e = g.add_edge(node_map[row['source']], node_map[row['target']])
        eweight[e] = row['jaccard']
        ewidth[e] = max(0.5, row['jaccard'] * 8)
        ecolor[e] = [0.4, 0.4, 0.4, max(0.2, row['jaccard'])]
        elabel[e] = str(row['shared'])

    return g, {
        'vlabel': vlabel, 'vcolor': vcolor, 'vsize': vsize,
        'vfontsize': vfontsize,
        'eweight': eweight, 'ecolor': ecolor, 'ewidth': ewidth,
        'elabel': elabel
    }


# ── Rendering ────────────────────────────────────────────────────────────────

def render_graph(g, props, title, outpath_stem, output_size=(2400, 2400),
                 show_edge_labels=False, label_all_vertices=False):
    """Render graph to PDF and PNG using graph-tool's graph_draw."""
    if g.num_vertices() == 0:
        print(f"  [WARN] Empty graph for {title} — skipping")
        return

    # Choose layout
    n = g.num_vertices()
    if n <= 80:
        pos = gt.fruchterman_reingold_layout(g, n_iter=1000)
    else:
        # sfdp with edge weights as spring constants for better clustering
        pos = gt.sfdp_layout(g, eweight=props.get('eweight'),
                             cooling_step=0.9, epsilon=0.005,
                             K=3.0)

    # Build vertex text: filter labels based on context
    vtext = g.new_vp('string')
    for v in g.vertices():
        if label_all_vertices:
            vtext[v] = props['vlabel'][v]
        elif props.get('vfontsize') and props['vfontsize'][v] > 0:
            vtext[v] = props['vlabel'][v]
        else:
            vtext[v] = ''

    # Text position: -1 puts label outside (above); 0 = automatic
    text_pos = -1 if n < 150 else 0

    draw_kwargs = dict(
        pos=pos,
        vertex_text=vtext,
        vertex_fill_color=props['vcolor'],
        vertex_size=props['vsize'],
        vertex_font_size=props.get('vfontsize', 8),
        vertex_text_color=props.get('vtextcolor'),
        vertex_text_position=text_pos,
        vertex_halo=props.get('vhalo'),
        vertex_halo_color=props.get('vhalocolor'),
        vertex_halo_size=1.5,
        edge_color=props['ecolor'],
        edge_pen_width=props['ewidth'],
        edge_marker_size=10 if g.is_directed() else 0,
        output_size=output_size,
        fit_view=0.9,
    )

    if show_edge_labels and 'elabel' in props:
        draw_kwargs['edge_text'] = props['elabel']
        draw_kwargs['edge_font_size'] = 8
        draw_kwargs['edge_text_color'] = [0.2, 0.2, 0.2, 0.9]

    for fmt in ['pdf', 'png']:
        outpath = f"{outpath_stem}.{fmt}"
        gt.graph_draw(g, output=outpath, **draw_kwargs)
        print(f"  Saved: {outpath}")


def render_legend(cell_type_colors, outpath, title="Cell Type Legend"):
    """Render a standalone legend mapping colors to cell types."""
    fig, ax = plt.subplots(figsize=(4, max(3, len(cell_type_colors) * 0.3)))
    for i, (ct, color) in enumerate(sorted(cell_type_colors.items())):
        ax.barh(i, 1, color=color, height=0.8)
        ax.text(1.05, i, ct, va='center', fontsize=8)
    ax.set_xlim(0, 3)
    ax.set_ylim(-0.5, len(cell_type_colors) - 0.5)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='SCENIC+ GRN network visualization with graph-tool')
    parser.add_argument('--eregulon-tsv', required=True)
    parser.add_argument('--aucell-direct', required=True)
    parser.add_argument('--rss-csv', required=True)
    parser.add_argument('--tf-to-gene-adj', required=True)
    parser.add_argument('--cell-type-key', default='cell_type')
    parser.add_argument('--top-n', type=int, default=5)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=== SCENIC+ GRN Network Visualization (graph-tool) ===")

    # Load data
    print("Loading eRegulon data...")
    rss_df = load_rss(args.rss_csv)
    edges = build_eregulon_edges(args.eregulon_tsv, regulation_filter='+/+')
    print(f"  RSS matrix: {rss_df.shape[0]} eRegulons × {rss_df.shape[1]} cell types")
    print(f"  Activating edges (TF→Region→Gene): {len(edges)}")

    cell_types = list(rss_df.columns)
    cell_type_colors = build_color_map(cell_types)

    # ── 1. Per-cell-type GRNs ─────────────────────────────────────────
    print(f"\n[1/3] Building per-cell-type GRNs ({len(cell_types)} cell types)...")
    ct_outdir = os.path.join(args.outdir, 'per_celltype')
    os.makedirs(ct_outdir, exist_ok=True)

    for ct in cell_types:
        ct_safe = (ct.replace('/', '_').replace(' ', '_')
                   .replace('+', 'pos').replace('-', 'neg'))
        top_tfs = get_top_tfs(rss_df, cell_type=ct, top_n=args.top_n)
        edge_agg = aggregate_edges(edges, top_tfs)
        if edge_agg is None or edge_agg.empty:
            print(f"  {ct}: no edges — skipping")
            continue

        print(f"  {ct}: {len(top_tfs)} TFs, {len(edge_agg)} edges")
        g, props = build_gt_tf_gene_graph(edge_agg, rss_df, cell_type_colors)

        # For per-cell-type: label target genes too (smaller graphs)
        for v in g.vertices():
            if props['vtype'][v] == 'target':
                props['vfontsize'][v] = 5
                props['vsize'][v] = 8

        render_graph(g, props, f'GRN — {ct}',
                     os.path.join(ct_outdir, f'grn_{ct_safe}'),
                     output_size=(1600, 1600), label_all_vertices=True)

    # ── 2. Holistic GRN ──────────────────────────────────────────────
    print("\n[2/3] Building holistic GRN (top TFs across all cell types)...")
    all_top_tfs = get_top_tfs(rss_df, cell_type=None, top_n=args.top_n)
    # Limit targets per TF to keep the graph readable
    edge_agg = aggregate_edges(edges, all_top_tfs, max_targets_per_tf=15)
    if edge_agg is not None and not edge_agg.empty:
        n_tfs = len(all_top_tfs)
        n_targets = edge_agg['Gene'].nunique()
        print(f"  Nodes: {n_tfs} TFs + {n_targets} targets")
        print(f"  Edges: {len(edge_agg)}")

        g, props = build_gt_tf_gene_graph(edge_agg, rss_df, cell_type_colors)
        render_graph(g, props,
                     'Holistic GRN — Top eRegulons',
                     os.path.join(args.outdir, 'grn_holistic'),
                     output_size=(3000, 3000))

        # Save legend
        # Filter to only cell types that appear as best_ct for at least one TF
        used_cts = set()
        for tf in all_top_tfs:
            ct, _ = get_best_celltype(tf, rss_df)
            used_cts.add(ct)
        legend_colors = {ct: cell_type_colors[ct] for ct in used_cts
                         if ct in cell_type_colors}
        render_legend(legend_colors,
                      os.path.join(args.outdir, 'grn_holistic_legend.png'),
                      title='TF Cell Type Specificity')
    else:
        print("  [WARN] No edges for holistic GRN — skipping")

    # ── 3. TF co-regulation network ──────────────────────────────────
    print("\n[3/3] Building TF co-regulation network...")
    coreg_g, coreg_props = build_gt_coregulation_graph(
        edges, rss_df, cell_type_colors, min_shared=5)
    if coreg_g is not None:
        print(f"  TFs: {coreg_g.num_vertices()}, edges: {coreg_g.num_edges()}")
        render_graph(coreg_g, coreg_props,
                     'TF Co-Regulation (≥5 shared targets)',
                     os.path.join(args.outdir, 'grn_tf_coregulation'),
                     output_size=(1600, 1600),
                     show_edge_labels=True, label_all_vertices=True)
    else:
        print("  [WARN] No TF pairs with ≥5 shared targets — skipping")

    print(f"\n=== GRN visualization complete. Outputs in {args.outdir} ===")


if __name__ == '__main__':
    main()
