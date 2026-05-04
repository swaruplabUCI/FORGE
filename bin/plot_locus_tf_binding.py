#!/usr/bin/env python3
"""
plot_locus_tf_binding.py — Shi Fig 4E approximation.

For each curated gene from candidates.json, render a bar chart of differential
TF activity (ChromVAR-based) for the TFs predicted to regulate that gene
according to the TF→target adjacency. Without GWAS SuSiE windows we substitute
the Cicero CCAN window of the gene as the locus restriction (the adjacency
already captures TF→gene edges anchored to that window).

Inputs:
  --adjacency    tf_gene_adjacency.tsv  (TF × target_gene × cell_type)
  --diff-tf-dir  results/differential_tf/  (tf_differential_<ct>_*.csv)
  --candidates   candidates.json from select_shi_candidates.py
  --treatment, --control   condition labels for the differential CSV suffix

Outputs:
  <gene>_locus_tf_binding.{pdf,png}    — one figure per curated gene
  locus_tf_binding_summary.tsv
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--adjacency", required=True)
    p.add_argument("--diff-tf-dir", required=True)
    p.add_argument("--candidates", required=True)
    p.add_argument("--treatment", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--top-tfs-per-gene", type=int, default=15)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def _label(motif):
    parts = str(motif).replace("\t", " ").split()
    if len(parts) >= 2 and parts[0].startswith("MA"):
        return parts[1]
    return parts[-1] if parts else str(motif)


def load_diff_tables(diff_dir, treatment, control):
    """Returns dict ct -> DataFrame(motif, label, logfoldchange, pval_adj)."""
    suffix = f"_{treatment}_vs_{control}.csv"
    tables = {}
    for csv in Path(diff_dir).glob(f"tf_differential_*{suffix}"):
        ct = csv.name[len("tf_differential_"):-len(suffix)]
        df = pd.read_csv(csv)
        if df.empty:
            continue
        df["label"] = df["motif"].map(_label).str.upper()
        tables[ct] = df[["label", "logfoldchange", "pval_adj"]]
    return tables


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cands = json.loads(Path(args.candidates).read_text())
    genes = []
    if isinstance(cands.get("top_loci"), dict):
        genes = cands["top_loci"].get("selected_genes", [])
    if not genes:
        raise SystemExit("[4E] no curated genes in candidates.json:top_loci.selected_genes")

    adj = pd.read_csv(args.adjacency, sep="\t")
    diff_tables = load_diff_tables(args.diff_tf_dir, args.treatment, args.control)
    if not diff_tables:
        raise SystemExit("[4E] no diff TF CSVs found")

    summary_rows = []
    for gene in genes:
        sub = adj[adj["target_gene"].astype(str) == str(gene)].copy()
        if sub.empty:
            print(f"[4E] {gene}: no TF predicted to target it; skip", flush=True)
            continue

        # Join differential TF activity per (TF, cell_type) using the adjacency cell_type
        sub["TF_upper"] = sub["TF"].astype(str).str.upper()
        rows = []
        for _, r in sub.iterrows():
            ct = r["cell_type"]
            ct_alt = ct.replace(" ", "_")
            for ct_try in (ct, ct_alt):
                tab = diff_tables.get(ct_try)
                if tab is None:
                    continue
                hit = tab[tab["label"] == r["TF_upper"]]
                if not hit.empty:
                    h = hit.iloc[0]
                    rows.append({
                        "gene": gene,
                        "tf": r["TF"],
                        "cell_type": ct_try,
                        "logfoldchange": float(h["logfoldchange"]),
                        "pval_adj": float(h["pval_adj"]),
                        "edge_type": r.get("edge_type", ""),
                        "evidence_count": r.get("evidence_count", np.nan),
                    })
                    break
        if not rows:
            print(f"[4E] {gene}: TFs in adjacency but none matched in diff-TF; skip",
                  flush=True)
            continue

        df = pd.DataFrame(rows)
        df["abs_lfc"] = df["logfoldchange"].abs()
        df = df.sort_values("abs_lfc", ascending=False).head(args.top_tfs_per_gene)
        df = df.sort_values("logfoldchange")

        labels = df["tf"].astype(str) + " — " + df["cell_type"].astype(str)
        fig, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(df))))
        colors = ["#d62728" if v >= 0 else "#1f77b4" for v in df["logfoldchange"]]
        ax.barh(labels, df["logfoldchange"], color=colors, edgecolor="white",
                linewidth=0.4)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel(f"log2 FC ({args.treatment} / {args.control})  — ChromVAR")
        ax.set_title(f"TF activity at {gene} CCAN window\n"
                     f"top {len(df)} TFs by |log2FC|", fontsize=10)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"{gene}_locus_tf_binding.{ext}", dpi=200,
                       bbox_inches="tight")
        plt.close(fig)
        summary_rows.extend(df.to_dict("records"))
        print(f"[4E] {gene}: {len(df)} TFs plotted", flush=True)

    if not summary_rows:
        raise SystemExit("[4E] no curated gene produced a chart")
    pd.DataFrame(summary_rows).to_csv(outdir / "locus_tf_binding_summary.tsv",
                                      sep="\t", index=False)
    print(f"[4E] wrote summary for {len({r['gene'] for r in summary_rows})} genes",
          flush=True)


if __name__ == "__main__":
    main()
