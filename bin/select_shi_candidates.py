#!/usr/bin/env python3
"""
select_shi_candidates.py — pick our analogues for the Shi et al. curated-example panels.

Outputs a single JSON manifest at <outdir>/candidates.json with:

  marker_genes              hard-coded mouse markers for Fig 1E
  top_tfs_per_celltype      top-N TFs per broad cell type by |logfoldchange| at FDR<cutoff
  top_loci_per_celltype     top-N gene loci per broad cell type by |delta_coaccess| on
                            promoter-linked CCANs (uses stratified Cicero output)

Inputs:
  --diff-tf-dir         results/differential_tf/  (tf_differential_<ct>_<trt>_vs_<ctrl>.csv)
  --connections-ctrl    Cicero connections TSV for control condition
  --connections-trt     Cicero connections TSV for treatment condition
  --gene-links          ccan_enhancer_gene_links.tsv (from extract_ccan_enhancers, control)
  --gene-links-trt      ccan_enhancer_gene_links_<trt>.tsv (optional)
  --seed-genes          comma-separated AD-relevant seed genes (always retained)
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


MOUSE_MARKERS = {
    "Astrocyte":      "Gfap",
    "Excitatory":     "Slc17a7",
    "Inhibitory":     "Gad2",
    "Microglia":      "Csf1r",
    "Oligodendrocyte": "Mobp",
    "OPC":            "Pdgfra",
    "Vascular":       "Cldn5",
}

DEFAULT_SEED = "Adam10,App,Apoe,Mapt,Trem2,Cx3cr1,Bin1,Picalm,Spi1,Srebf1,Hes1,Maf"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--diff-tf-dir", required=True)
    p.add_argument("--connections-ctrl", required=True)
    p.add_argument("--connections-trt", required=True)
    p.add_argument("--gene-links", required=True,
                   help="ccan_enhancer_gene_links.tsv (peak->gene mapping)")
    p.add_argument("--gene-links-trt", default=None,
                   help="treatment-condition gene links (optional; falls back to --gene-links)")
    p.add_argument("--treatment", required=True, help="Treatment label, e.g. SREBF1_OE")
    p.add_argument("--control", required=True, help="Control label, e.g. 5xFAD")
    p.add_argument("--seed-genes", default=DEFAULT_SEED)
    p.add_argument("--tf-allowlist", default=None,
                   help="Optional file (TSV/CSV) whose first column lists TF symbols. "
                        "When given, only motifs whose label intersects this list are "
                        "considered for top_tfs_per_celltype. Useful to restrict to "
                        "TFs present in the tf_gene_adjacency.tsv.")
    p.add_argument("--adjacency", default=None,
                   help="Optional tf_gene_adjacency.tsv. When given, the top loci list "
                        "is augmented with the most-targeted genes from the adjacency, "
                        "so that 4E (TF binding bars) has data for every selected gene.")
    p.add_argument("--top-tfs", type=int, default=3)
    p.add_argument("--top-loci", type=int, default=5)
    p.add_argument("--fdr-cutoff", type=float, default=0.05)
    p.add_argument("--out", default="candidates.json")
    return p.parse_args()


def _motif_label(motif):
    parts = str(motif).replace("\t", " ").split()
    if len(parts) >= 2 and parts[0].startswith("MA"):
        return parts[1]
    return parts[-1] if parts else str(motif)


def collect_top_tfs(diff_tf_dir, top_n, fdr, treatment, control, allowlist=None):
    rows = []
    suffix = f"_{treatment}_vs_{control}.csv"
    for csv in sorted(Path(diff_tf_dir).glob(f"tf_differential_*{suffix}")):
        ct = csv.name[len("tf_differential_"):-len(suffix)]
        df = pd.read_csv(csv)
        if df.empty or "pval_adj" not in df.columns:
            continue
        df["motif"] = df["motif"].astype(str).str.replace(r"\s+", " ", regex=True)
        df["label"] = df["motif"].map(_motif_label)
        if allowlist is not None:
            df = df[df["label"].str.upper().isin(allowlist)]
            if df.empty:
                continue
        sig = df[df["pval_adj"] < fdr].copy()
        if sig.empty:
            sig = df.head(top_n).copy()
        sig["abs_lfc"] = sig["logfoldchange"].abs()
        top = sig.sort_values("abs_lfc", ascending=False).head(top_n)
        rows.append({
            "cell_type": ct,
            "treatment": treatment,
            "control":   control,
            "top_tfs":   top[["motif", "logfoldchange", "pval_adj"]].to_dict("records"),
        })
    return rows


def load_connections(path):
    df = pd.read_csv(path, sep="\t")
    cols = list(df.columns)
    df = df.rename(columns={cols[0]: "Peak1", cols[1]: "Peak2", cols[2]: "coaccess"})
    df["coaccess"] = pd.to_numeric(df["coaccess"], errors="coerce")
    df = df.dropna(subset=["coaccess"])
    return df


def _links_to_map(path):
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        return {}
    df["peak_id"] = (df["chr"].astype(str) + ":" +
                    df["start"].astype(str) + "-" + df["end"].astype(str))
    return dict(zip(df["peak_id"], df["gene"]))


def collect_top_loci(conn_ctrl_path, conn_trt_path, gene_links_ctrl, gene_links_trt,
                    top_n, seed_genes):
    ctrl = load_connections(conn_ctrl_path)
    trt = load_connections(conn_trt_path)
    map_ctrl = _links_to_map(gene_links_ctrl)
    map_trt = _links_to_map(gene_links_trt) if gene_links_trt else map_ctrl
    if not map_ctrl and not map_trt:
        return {"note": "no gene mappings available; skipping"}

    def join_delta(conn_df, mapping):
        conn = conn_df.copy()
        conn["g1"] = conn["Peak1"].map(mapping)
        conn["g2"] = conn["Peak2"].map(mapping)
        conn["gene"] = conn["g1"].fillna(conn["g2"])
        return conn[["Peak1", "Peak2", "coaccess", "gene"]].dropna(subset=["gene"])

    c = join_delta(ctrl, map_ctrl).rename(columns={"coaccess": "coaccess_ctrl"})
    t = join_delta(trt, map_trt).rename(columns={"coaccess": "coaccess_trt"})
    merged = c.merge(t, on=["Peak1", "Peak2", "gene"], how="outer")
    merged["coaccess_ctrl"] = merged["coaccess_ctrl"].fillna(0)
    merged["coaccess_trt"] = merged["coaccess_trt"].fillna(0)
    merged["delta"] = (merged["coaccess_trt"] - merged["coaccess_ctrl"]).abs()

    per_gene = merged.groupby("gene")["delta"].sum().reset_index()
    per_gene = per_gene.sort_values("delta", ascending=False)

    seed = [g.strip() for g in seed_genes.split(",") if g.strip()]
    seed_present = [g for g in seed if g in set(per_gene["gene"])]
    seed_df = per_gene[per_gene["gene"].isin(seed_present)]
    top_df = per_gene.head(top_n)
    selected = pd.concat([seed_df, top_df]).drop_duplicates("gene").head(top_n + len(seed_present))
    return {
        "global_top": top_df.head(top_n).to_dict("records"),
        "seed_present": seed_df.to_dict("records"),
        "selected_genes": selected["gene"].tolist(),
    }


def _load_tf_allowlist(path):
    if not path:
        return None
    out = set()
    with open(path) as f:
        for ln in f:
            tok = ln.strip().split("\t")[0].split(",")[0]
            if tok and not tok.lower().startswith("tf"):
                out.add(tok.upper())
    return out


def main():
    args = parse_args()
    allowlist = _load_tf_allowlist(args.tf_allowlist)
    top_loci = collect_top_loci(
        args.connections_ctrl, args.connections_trt,
        args.gene_links, args.gene_links_trt,
        args.top_loci, args.seed_genes)

    # Augment selected_genes with top-targeted adjacency genes (drives 4E coverage)
    if args.adjacency and isinstance(top_loci, dict):
        try:
            adj = pd.read_csv(args.adjacency, sep="\t")
            top_adj = (adj.groupby("target_gene")["evidence_count"].sum()
                       .sort_values(ascending=False).head(args.top_loci).index.tolist())
            existing = set(top_loci.get("selected_genes", []))
            top_loci.setdefault("adjacency_top_genes", top_adj)
            top_loci["selected_genes"] = list(dict.fromkeys(
                top_loci.get("selected_genes", []) + top_adj))
        except Exception as e:
            print(f"[select_candidates] adjacency augmentation failed: {e}", flush=True)

    out = {
        "marker_genes": MOUSE_MARKERS,
        "top_tfs_per_celltype": collect_top_tfs(
            args.diff_tf_dir, args.top_tfs, args.fdr_cutoff,
            args.treatment, args.control, allowlist),
        "top_loci": top_loci,
        "params": {
            "top_tfs": args.top_tfs,
            "top_loci": args.top_loci,
            "fdr_cutoff": args.fdr_cutoff,
            "seed_genes": args.seed_genes,
        },
    }
    Path(args.out).write_text(json.dumps(out, indent=2))
    print(f"[select_candidates] wrote {args.out}")
    print(f"[select_candidates]   top_tfs cell types: "
          f"{[r['cell_type'] for r in out['top_tfs_per_celltype']]}")
    if isinstance(out["top_loci"], dict):
        print(f"[select_candidates]   selected loci: {out['top_loci']['selected_genes']}")


if __name__ == "__main__":
    main()
