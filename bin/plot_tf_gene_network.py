#!/usr/bin/env python3
"""
plot_tf_gene_network.py — graph-tool visualization of SDas TF-gene network.

Consumes tf_gene_adjacency.tsv from BUILD_TF_GENE_NETWORK (Shi et al. TF_Net
equivalent, built from scPrinter binding + Cicero co-accessibility) and
renders per-cell-type GRNs, a holistic GRN, and a TF co-regulation network.

Input schema (tab-separated):
    TF, target_gene, cell_type, edge_type, evidence_count, coaccessibility
Where edge_type ∈ {promoter_direct, enhancer_indirect}.

Outputs (PDF + PNG):
    per_celltype/grn_<celltype>.{pdf,png}
    grn_holistic.{pdf,png}
    grn_holistic_legend.png
    grn_tf_coregulation.{pdf,png}
"""
import argparse
import os
import warnings

warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import graph_tool.all as gt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def build_color_map(categories):
    n = len(categories)
    if n <= 10:
        cmap = plt.get_cmap('tab10')
    elif n <= 20:
        cmap = plt.get_cmap('tab20')
    else:
        cmap = plt.get_cmap('turbo')
    return {cat: cmap(i / max(n - 1, 1)) for i, cat in enumerate(sorted(categories))}


def load_edges(path):
    df = pd.read_csv(path, sep='\t')
    required = {'TF', 'target_gene', 'cell_type', 'edge_type',
                'evidence_count', 'coaccessibility'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    df['evidence_count'] = pd.to_numeric(df['evidence_count'], errors='coerce').fillna(0)
    df['coaccessibility'] = pd.to_numeric(df['coaccessibility'], errors='coerce').fillna(0)
    df['weight'] = df['evidence_count'] * df['coaccessibility'].clip(lower=0.1)
    return df


def rank_tfs_per_celltype(edges, cell_type, top_n):
    sub = edges[edges['cell_type'] == cell_type]
    if sub.empty:
        return set()
    rank = (sub.groupby('TF')
              .agg(n_targets=('target_gene', 'nunique'),
                   mean_weight=('weight', 'mean'))
              .assign(score=lambda x: x['n_targets'] * x['mean_weight'].clip(lower=0.1))
              .sort_values('score', ascending=False))
    return set(rank.head(top_n).index)


def get_dominant_celltype(tf, edges):
    sub = edges[edges['TF'] == tf]
    if sub.empty:
        return 'unknown'
    return (sub.groupby('cell_type')['target_gene']
               .nunique()
               .idxmax())


def aggregate_tf_gene(edges, tfs, max_targets_per_tf=None):
    sub = edges[edges['TF'].isin(tfs)]
    if sub.empty:
        return None
    agg = (sub.groupby(['TF', 'target_gene'])
              .agg(weight=('weight', 'sum'),
                   evidence=('evidence_count', 'sum'),
                   coaccess=('coaccessibility', 'max'),
                   n_cts=('cell_type', 'nunique'),
                   edge_type=('edge_type', lambda s: 'promoter_direct'
                              if (s == 'promoter_direct').any() else 'enhancer_indirect'))
              .reset_index())
    if max_targets_per_tf:
        agg = (agg.groupby('TF', group_keys=False)
                  .apply(lambda g: g.nlargest(max_targets_per_tf, 'weight'))
                  .reset_index(drop=True))
    return agg


def build_tf_gene_graph(edge_agg, edges_all, cell_type_colors):
    g = gt.Graph(directed=True)

    vlabel = g.new_vp('string')
    vtype = g.new_vp('string')
    vcolor = g.new_vp('vector<double>')
    vsize = g.new_vp('double')
    vfontsize = g.new_vp('double')
    vshape = g.new_vp('string')
    vtextcolor = g.new_vp('vector<double>')
    vhalo = g.new_vp('bool')
    vhalocolor = g.new_vp('vector<double>')

    eweight = g.new_ep('double')
    ecolor = g.new_ep('vector<double>')
    ewidth = g.new_ep('double')

    node_map = {}

    def get_or_add(name):
        if name not in node_map:
            v = g.add_vertex()
            node_map[name] = v
            vlabel[v] = name
        return node_map[name]

    tf_nodes = edge_agg['TF'].unique()
    for tf in tf_nodes:
        v = get_or_add(tf)
        vtype[v] = 'TF'
        dom_ct = get_dominant_celltype(tf, edges_all)
        ct_color = cell_type_colors.get(dom_ct, (0.5, 0.5, 0.5, 1.0))
        vcolor[v] = list(ct_color)
        n_targets = edge_agg[edge_agg['TF'] == tf]['target_gene'].nunique()
        vsize[v] = max(20, min(70, 20 + n_targets * 0.4))
        vfontsize[v] = max(9, min(16, 9 + n_targets * 0.12))
        vshape[v] = 'circle'
        vtextcolor[v] = [0.0, 0.0, 0.0, 1.0]
        vhalo[v] = True
        vhalocolor[v] = list(ct_color[:3]) + [0.3]

    for _, row in edge_agg.iterrows():
        tf_v = node_map[row['TF']]
        gene_v = get_or_add(row['target_gene'])
        if vtype[gene_v] != 'TF':
            vtype[gene_v] = 'target'
            vcolor[gene_v] = [0.65, 0.75, 0.88, 0.9]
            vsize[gene_v] = 8
            vfontsize[gene_v] = 0
            vshape[gene_v] = 'circle'
            vtextcolor[gene_v] = [0.2, 0.2, 0.2, 1.0]
            vhalo[gene_v] = False
            vhalocolor[gene_v] = [0.0, 0.0, 0.0, 0.0]

        e = g.add_edge(tf_v, gene_v)
        eweight[e] = float(row['weight'])
        tf_color = list(vcolor[tf_v])
        # promoter_direct edges opaque; enhancer_indirect translucent
        alpha = 0.55 if row['edge_type'] == 'promoter_direct' else 0.28
        ecolor[e] = tf_color[:3] + [alpha]
        ewidth[e] = max(0.5, min(3.0, np.log1p(row['weight']) * 0.8))

    return g, {
        'vlabel': vlabel, 'vtype': vtype, 'vcolor': vcolor, 'vsize': vsize,
        'vfontsize': vfontsize, 'vshape': vshape, 'vtextcolor': vtextcolor,
        'vhalo': vhalo, 'vhalocolor': vhalocolor,
        'eweight': eweight, 'ecolor': ecolor, 'ewidth': ewidth
    }


def build_coregulation_graph(edges_all, tfs_of_interest, cell_type_colors, min_shared=5):
    sub = edges_all[edges_all['TF'].isin(tfs_of_interest)]
    tf_targets = sub.groupby('TF')['target_gene'].apply(set).to_dict()
    tfs = list(tf_targets.keys())

    coedges = []
    for i in range(len(tfs)):
        ti = tf_targets[tfs[i]]
        for j in range(i + 1, len(tfs)):
            tj = tf_targets[tfs[j]]
            shared = ti & tj
            if len(shared) >= min_shared:
                jaccard = len(shared) / len(ti | tj)
                coedges.append({'source': tfs[i], 'target': tfs[j],
                                'shared': len(shared), 'jaccard': jaccard})

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
    tfs_in_net = set(coedge_df['source']) | set(coedge_df['target'])
    for tf in sorted(tfs_in_net):
        v = g.add_vertex()
        node_map[tf] = v
        vlabel[v] = tf
        dom_ct = get_dominant_celltype(tf, edges_all)
        ct_color = cell_type_colors.get(dom_ct, (0.5, 0.5, 0.5, 1.0))
        vcolor[v] = list(ct_color)
        n_tgt = len(tf_targets.get(tf, []))
        vsize[v] = max(20, min(70, 20 + n_tgt * 0.05))
        vfontsize[v] = 10

    for _, row in coedge_df.iterrows():
        e = g.add_edge(node_map[row['source']], node_map[row['target']])
        eweight[e] = float(row['jaccard'])
        ewidth[e] = max(0.5, row['jaccard'] * 8)
        ecolor[e] = [0.4, 0.4, 0.4, max(0.2, row['jaccard'])]
        elabel[e] = str(row['shared'])

    return g, {
        'vlabel': vlabel, 'vcolor': vcolor, 'vsize': vsize,
        'vfontsize': vfontsize,
        'eweight': eweight, 'ecolor': ecolor, 'ewidth': ewidth,
        'elabel': elabel
    }


def _layout_with_component_packing(g, eweight=None):
    """Layout each connected component separately, then pack in a grid.

    graph-tool's force-directed layouts spread disconnected components to
    infinity; packing brings them into a compact grid so canvas space isn't
    wasted on sparse graphs (common for TF co-regulation / per-cell-type
    GRNs where TFs don't share targets).
    """
    comp, _ = gt.label_components(g, directed=False)
    n_comp = int(comp.a.max()) + 1 if g.num_vertices() > 0 else 0
    pos = g.new_vp('vector<double>')

    if n_comp <= 1:
        base = (gt.fruchterman_reingold_layout(g, n_iter=1000)
                if g.num_vertices() <= 80
                else gt.sfdp_layout(g, eweight=eweight, cooling_step=0.9,
                                    epsilon=0.005, K=3.0))
        for v in g.vertices():
            pos[v] = list(base[v])
        return pos

    # Per-component layout + bounding box
    comp_sizes = np.bincount(comp.a)
    order = np.argsort(-comp_sizes)  # biggest first
    boxes = []
    for c in order:
        vfilt = g.new_vp('bool')
        for v in g.vertices():
            vfilt[v] = (comp[v] == c)
        sub = gt.GraphView(g, vfilt=vfilt)
        n_sub = sub.num_vertices()
        if n_sub == 0:
            continue
        sub_pos = (gt.fruchterman_reingold_layout(sub, n_iter=500)
                   if n_sub <= 50
                   else gt.sfdp_layout(sub, cooling_step=0.9, K=2.0))
        coords = np.array([list(sub_pos[v]) for v in sub.vertices()])
        mn, mx = coords.min(axis=0), coords.max(axis=0)
        extent = np.maximum(mx - mn, 1.0)
        # Normalize each component to a unit box, then scale by sqrt(n_sub)
        scale = np.sqrt(max(n_sub, 1))
        norm = (coords - mn) / extent * scale
        boxes.append((c, list(sub.vertices()), norm, scale))

    # Grid arrangement: square-ish layout
    n_boxes = len(boxes)
    n_cols = int(np.ceil(np.sqrt(n_boxes)))
    max_scale = max((b[3] for b in boxes), default=1.0)
    pad = max_scale * 0.15

    for i, (_, verts, norm, scale) in enumerate(boxes):
        row, col = divmod(i, n_cols)
        dx = col * (max_scale + pad)
        dy = -row * (max_scale + pad)  # top-down
        for v_idx, v in enumerate(verts):
            pos[v] = [float(norm[v_idx, 0] + dx), float(norm[v_idx, 1] + dy)]
    return pos


def render_graph(g, props, outpath_stem, output_size=(2400, 2400),
                 show_edge_labels=False, label_all_vertices=False):
    if g.num_vertices() == 0:
        print(f"  [WARN] Empty graph — skipping {outpath_stem}")
        return

    pos = _layout_with_component_packing(g, eweight=props.get('eweight'))
    n = g.num_vertices()

    vtext = g.new_vp('string')
    for v in g.vertices():
        if label_all_vertices:
            vtext[v] = props['vlabel'][v]
        elif props.get('vfontsize') and props['vfontsize'][v] > 0:
            vtext[v] = props['vlabel'][v]
        else:
            vtext[v] = ''

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


def sanitize(name):
    return (str(name).replace('/', '_').replace(' ', '_')
                     .replace('+', 'pos').replace('-', 'neg'))


def main():
    ap = argparse.ArgumentParser(description='SDas TF-gene network visualization (graph-tool)')
    ap.add_argument('--adjacency', required=True,
                    help='tf_gene_adjacency.tsv from BUILD_TF_GENE_NETWORK')
    ap.add_argument('--top-n', type=int, default=5,
                    help='Top N TFs per cell type (default: 5)')
    ap.add_argument('--max-targets-per-tf', type=int, default=15,
                    help='Max targets per TF in holistic GRN (default: 15)')
    ap.add_argument('--min-shared-coreg', type=int, default=5,
                    help='Min shared targets for TF co-reg edges (default: 5)')
    ap.add_argument('--outdir', default='.')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=== SDas TF-Gene Network Visualization (graph-tool) ===")

    edges_all = load_edges(args.adjacency)
    if edges_all.empty:
        print("[WARN] Adjacency is empty — nothing to plot.")
        return

    cell_types = sorted(edges_all['cell_type'].unique())
    cell_type_colors = build_color_map(cell_types)
    print(f"  Loaded {len(edges_all)} edges: "
          f"{edges_all['TF'].nunique()} TFs × {edges_all['target_gene'].nunique()} genes "
          f"× {len(cell_types)} cell types")

    # 1. Per-cell-type GRNs
    print(f"\n[1/3] Per-cell-type GRNs (top-{args.top_n} TFs each)")
    ct_outdir = os.path.join(args.outdir, 'per_celltype')
    os.makedirs(ct_outdir, exist_ok=True)

    for ct in cell_types:
        top_tfs = rank_tfs_per_celltype(edges_all, ct, args.top_n)
        if not top_tfs:
            print(f"  {ct}: no TFs — skipping")
            continue
        ct_edges = edges_all[edges_all['cell_type'] == ct]
        edge_agg = aggregate_tf_gene(ct_edges, top_tfs)
        if edge_agg is None or edge_agg.empty:
            print(f"  {ct}: no edges — skipping")
            continue

        print(f"  {ct}: {len(top_tfs)} TFs, {len(edge_agg)} edges")
        g, props = build_tf_gene_graph(edge_agg, edges_all, cell_type_colors)

        for v in g.vertices():
            if props['vtype'][v] == 'target':
                props['vfontsize'][v] = 5
                props['vsize'][v] = 8

        stem = os.path.join(ct_outdir, f'grn_{sanitize(ct)}')
        render_graph(g, props, stem, output_size=(1600, 1600),
                     label_all_vertices=True)

    # 2. Holistic GRN
    print("\n[2/3] Holistic GRN (union of per-celltype top TFs)")
    all_top_tfs = set()
    for ct in cell_types:
        all_top_tfs |= rank_tfs_per_celltype(edges_all, ct, args.top_n)

    edge_agg = aggregate_tf_gene(edges_all, all_top_tfs,
                                 max_targets_per_tf=args.max_targets_per_tf)
    if edge_agg is not None and not edge_agg.empty:
        print(f"  Nodes: {len(all_top_tfs)} TFs + "
              f"{edge_agg['target_gene'].nunique()} targets; "
              f"Edges: {len(edge_agg)}")
        g, props = build_tf_gene_graph(edge_agg, edges_all, cell_type_colors)
        render_graph(g, props, os.path.join(args.outdir, 'grn_holistic'),
                     output_size=(3000, 3000))

        used_cts = {get_dominant_celltype(tf, edges_all) for tf in all_top_tfs}
        legend_colors = {ct: cell_type_colors[ct] for ct in used_cts
                         if ct in cell_type_colors}
        render_legend(legend_colors,
                      os.path.join(args.outdir, 'grn_holistic_legend.png'),
                      title='TF Dominant Cell Type')
    else:
        print("  [WARN] No edges for holistic GRN")

    # 3. TF co-regulation network
    print(f"\n[3/3] TF co-regulation (≥{args.min_shared_coreg} shared targets)")
    coreg_g, coreg_props = build_coregulation_graph(
        edges_all, all_top_tfs, cell_type_colors,
        min_shared=args.min_shared_coreg)
    if coreg_g is not None:
        print(f"  TFs: {coreg_g.num_vertices()}, edges: {coreg_g.num_edges()}")
        render_graph(coreg_g, coreg_props,
                     os.path.join(args.outdir, 'grn_tf_coregulation'),
                     output_size=(1600, 1600),
                     show_edge_labels=True, label_all_vertices=True)
    else:
        print(f"  [WARN] No TF pairs with ≥{args.min_shared_coreg} shared targets")

    print(f"\n=== Done. Outputs in {args.outdir} ===")


if __name__ == '__main__':
    main()
