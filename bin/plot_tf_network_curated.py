#!/usr/bin/env python3
"""
plot_tf_network_curated.py — Shi Fig 4C / 5B / 5D / 5F equivalent.

For each broad cell type, render a small focused TF→target subgraph using
the top-3 TFs from select_shi_candidates.py, drawn from the global
tf_gene_adjacency.tsv produced by build_tf_gene_network.py. Node shape
encodes the differential-TF-accessibility direction (square = up in
treatment, circle = up in control). Without snRNA, target nodes are not
DE-coloured — flagged in the legend.

Inputs:
  --adjacency      tf_gene_adjacency.tsv (TF, target_gene, cell_type, edge_type, ...)
  --candidates     candidates.json from select_shi_candidates.py
  --diff-tf-dir    results/differential_tf  (for direction symbol on TFs)
  --treatment, --control
"""

import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adjacency", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--diff-tf-dir", required=True)
    p.add_argument("--treatment", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--max-targets-per-tf", type=int, default=12)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def _label_text(motif):
    parts = str(motif).replace("\t", " ").split()
    if len(parts) >= 2 and parts[0].startswith("MA"):
        return parts[1]
    return parts[-1] if parts else str(motif)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adj = pd.read_csv(args.adjacency, sep="\t")
    cands = json.loads(Path(args.candidates).read_text())
    suffix = f"_{args.treatment}_vs_{args.control}.csv"

    rendered = []
    for entry in cands.get("top_tfs_per_celltype", []):
        ct = entry["cell_type"]
        tf_records = entry.get("top_tfs", [])
        if not tf_records:
            continue
        tf_labels = [_label_text(r["motif"]) for r in tf_records]
        tf_dir = {_label_text(r["motif"]): "up_trt" if r["logfoldchange"] > 0 else "up_ctrl"
                  for r in tf_records}

        # Find adjacency rows for this cell type whose TF (case-insensitive) matches
        # one of the curated TFs. The adjacency uses TF gene symbols; motif names like
        # "MA1483.2 ELF2" already collapsed to "ELF2".
        # Adjacency uses spaces in cell-type names ("CNU-LGE GABA"), candidates
        # may use underscores ("CNU-LGE_GABA"). Try both forms.
        ct_alt = ct.replace("_", " ")
        adj_ct = adj[adj["cell_type"].isin([ct, ct_alt])].copy()
        adj_ct["TF_upper"] = adj_ct["TF"].astype(str).str.upper()
        wanted_upper = {t.upper() for t in tf_labels}
        sub = adj_ct[adj_ct["TF_upper"].isin(wanted_upper)].copy()
        if sub.empty:
            print(f"[curated_net] {ct} — none of {tf_labels} found in adjacency; skip",
                  flush=True)
            continue

        # Cap targets per TF for legibility, ranked by evidence_count then coaccess
        sort_cols = [c for c in ("evidence_count", "coaccessibility")
                    if c in sub.columns]
        sub = sub.sort_values(sort_cols, ascending=False) if sort_cols else sub
        kept = sub.groupby("TF").head(args.max_targets_per_tf)

        G = nx.DiGraph()
        for tf in tf_labels:
            G.add_node(tf, kind="tf", direction=tf_dir.get(tf, "neutral"))
        for _, row in kept.iterrows():
            tf = row["TF"]
            tg = row["target_gene"]
            if tg in [t for t in tf_labels]:
                continue
            G.add_node(tg, kind="target")
            etype = row.get("edge_type", "enhancer_linked")
            G.add_edge(tf, tg, edge_type=etype)

        if G.number_of_nodes() < 2:
            continue

        pos = nx.spring_layout(G, k=1.4 / max(1, G.number_of_nodes() ** 0.5),
                              seed=42, iterations=80)

        fig, ax = plt.subplots(figsize=(8, 7))
        # Edges
        prom_edges = [(u, v) for u, v, d in G.edges(data=True)
                     if d.get("edge_type") == "promoter_direct"]
        enh_edges = [(u, v) for u, v, d in G.edges(data=True)
                    if d.get("edge_type") != "promoter_direct"]
        nx.draw_networkx_edges(G, pos, edgelist=prom_edges, alpha=0.7,
                              edge_color="#444444", width=1.0,
                              arrowsize=10, ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=enh_edges, alpha=0.4,
                              edge_color="#aaaaaa", width=0.7, style="dashed",
                              arrowsize=10, ax=ax)

        # Target nodes (small light blue circles)
        target_nodes = [n for n, d in G.nodes(data=True) if d["kind"] == "target"]
        nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_color="#cfe6ff",
                              node_shape="o", node_size=320, edgecolors="#4477aa",
                              linewidths=0.8, ax=ax)

        # TF nodes — square if up in treatment, circle if up in control
        for tf in tf_labels:
            if tf not in G:
                continue
            shape = "s" if tf_dir[tf] == "up_trt" else "o"
            colour = "#ffd966" if tf_dir[tf] == "up_trt" else "#fcec88"
            nx.draw_networkx_nodes(G, pos, nodelist=[tf], node_color=colour,
                                  node_shape=shape, node_size=900,
                                  edgecolors="black", linewidths=1.2, ax=ax)

        nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)
        ax.set_title(
            f"{ct} — top {len(tf_labels)} TFs ({args.treatment} vs {args.control})\n"
            f"square = up in {args.treatment}, circle = up in {args.control}; "
            f"dashed = enhancer-linked",
            fontsize=9)
        ax.set_axis_off()
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{ct}_curated_network.{ext}", dpi=200,
                       bbox_inches="tight")
        plt.close(fig)
        rendered.append({"cell_type": ct, "tfs": tf_labels,
                        "n_targets": int(G.number_of_nodes()) - len(tf_labels)})
        print(f"[curated_net] {ct} — TFs={tf_labels} targets={G.number_of_nodes() - len(tf_labels)}",
              flush=True)

    if not rendered:
        raise SystemExit("[curated_net] no networks rendered")
    (outdir / "curated_networks_manifest.json").write_text(
        json.dumps({"rendered": rendered}, indent=2))
    print(f"[curated_net] wrote {len(rendered)} networks", flush=True)


if __name__ == "__main__":
    main()
