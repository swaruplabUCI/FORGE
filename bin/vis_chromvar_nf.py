#!/usr/bin/env python

import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import scprinter as scp
from statsmodels.stats.multitest import multipletests
from itertools import combinations

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--dev-h5ad", required=True,
        help="Path to chromVAR deviations h5ad (from GPU_CHROMVAR)"
    )
    p.add_argument(
        "--out-dir", required=True,
        help="Output directory for PDFs and TSVs"
    )
    p.add_argument(
        "--cell-type-key", default="cell_type",
        help="obs column to use as cell type (default: cell_type)"
    )
    return p.parse_args()

def select_motifs_by_tf(dev, tf_list):
    """Select one motif per TF (highest variance across cells)"""
    X = np.asarray(dev.X)
    hits = {}
    for tf in tf_list:
        ids = dev.var_names[dev.var["tf"].str.contains(tf, case=False, regex=False)]
        if len(ids) > 0:
            hits[tf] = ids.tolist()
    selected = {}
    for tf, ids in hits.items():
        if len(ids) == 1:
            selected[tf] = ids[0]
        else:
            idxs = [dev.var_names.get_loc(i) for i in ids]
            variances = np.nanvar(X[:, idxs], axis=0)
            best_id = ids[int(np.nanargmax(variances))]
            selected[tf] = best_id
    return selected

def find_cell_type(dev, search_terms):
    """
    Find cell type in dev.obs['cell_type'] using fuzzy matching.
    
    Args:
        dev: AnnData with .obs['cell_type']
        search_terms: list of strings to search for (case-insensitive)
    
    Returns:
        Exact cell type name found in data, or None
    """
    available_types = dev.obs['cell_type'].unique()
    
    for term in search_terms:
        for ct in available_types:
            if term.lower() in str(ct).lower():
                return str(ct)
    
    return None

def diff_within_celltype(dev, cell_type, group_col="Condition",
                         group1="Group1", group2="Group2",
                         motif_ids=None, n_iter=1000):
    """Run chromVAR permutation test for one cell type"""
    sub = dev[(dev.obs["cell_type"] == cell_type) & (~dev.obs[group_col].isna())].copy()
    if sub.n_obs == 0:
        raise ValueError(f"No cells for cell_type={cell_type}")
    if motif_ids is not None:
        sub = sub[:, motif_ids]
    n1 = (sub.obs[group_col].astype(str) == group1).sum()
    n2 = (sub.obs[group_col].astype(str) == group2).sum()
    print(f"{cell_type}: {group1}={n1}, {group2}={n2}, motifs={sub.n_vars}")
    if n1 < 10 or n2 < 10:
        print(f"Warning: small group sizes in {cell_type}.")
    res = scp.chromvar.permutation_test(
        adata=sub, groupby=group_col, group1=group1, group2=group2,
        n_iterations=n_iter
    )
    res = res.copy()
    res["TF"] = pd.Index(res.index).str.split("\t").str[-1]
    res["qval"] = multipletests(res["pval"].values, method="fdr_bh")[1]
    return res, sub

def main():
    args = parse_args()
    DEV_H5AD = args.dev_h5ad
    OUT_DIR = args.out_dir
    os.makedirs(OUT_DIR, exist_ok=True)

    # -------------------------------------------------------------------
    # Load chromVAR deviations
    # -------------------------------------------------------------------
    dev = ad.read_h5ad(DEV_H5AD)
    print(dev)

    # Remap cell type column if needed (so downstream code can use 'cell_type')
    ct_key = args.cell_type_key
    if ct_key != "cell_type" and ct_key in dev.obs.columns:
        dev.obs["cell_type"] = dev.obs[ct_key]
        print(f"✓ Mapped '{ct_key}' → 'cell_type'")

    # -------------------------------------------------------------------
    # Require condition column
    # -------------------------------------------------------------------
    condition_col_candidates = ["condition_group", "Condition", "condition", "group", "diagnosis"]
    cond_col = None
    
    for c in condition_col_candidates:
        if c in dev.obs.columns:
            cond_col = c
            break

    if cond_col is None:
        raise ValueError(
            f"No condition column found in dev.obs. Expected one of: {condition_col_candidates}\n"
            f"Available columns: {list(dev.obs.columns)}\n"
            f"Ensure your peak matrix includes condition metadata before running chromVAR."
        )

    print(f"✓ Using condition column: dev.obs['{cond_col}']")
    dev.obs["Condition"] = dev.obs[cond_col].astype("category")

    # Verify we have at least 2 groups
    n_groups = dev.obs["Condition"].nunique()
    if n_groups < 2:
        raise ValueError(
            f"Need at least 2 conditions for differential analysis, found {n_groups}\n"
            f"Condition values: {list(dev.obs['Condition'].unique())}"
        )

    print(f"✓ Condition counts: {dev.obs['Condition'].value_counts().to_dict()}")

    # -------------------------------------------------------------------
    # Require cell_type column
    # -------------------------------------------------------------------
    if "cell_type" not in dev.obs.columns:
        raise ValueError("dev.obs['cell_type'] not found; needed for per-cell-type analyses.")

    # -------------------------------------------------------------------
    # Parse TF names from motif var_names
    # -------------------------------------------------------------------
    dev.var["tf"] = pd.Index(dev.var_names).str.split("\t").str[-1]

    # -------------------------------------------------------------------
    # Glial TF query lists (unchanged - these are biological constants)
    # -------------------------------------------------------------------
    tf_queries_astro = ["SOX9","NFIA","NFIB","STAT3","CEBPB","CEBPD","RELA","NFKB1","JUN","FOS"]
    tf_queries_micro = ["SPI1","IRF8","CEBPA","CEBPB","CEBPD","STAT1","STAT2","RELA","NFKB1","JUN","FOS"]
    tf_queries_oligo = ["SOX10","OLIG2","MYRF","SOX8","NKX2-2"]

    sel_astro = select_motifs_by_tf(dev, tf_queries_astro)
    sel_micro = select_motifs_by_tf(dev, tf_queries_micro)
    sel_oligo = select_motifs_by_tf(dev, tf_queries_oligo)

    pd.DataFrame({"TF": list(sel_astro.keys()), "motif_id": list(sel_astro.values())}).to_csv(
        os.path.join(OUT_DIR, "selected_motifs_astro.tsv"), sep="\t", index=False)
    pd.DataFrame({"TF": list(sel_micro.keys()), "motif_id": list(sel_micro.values())}).to_csv(
        os.path.join(OUT_DIR, "selected_motifs_micro.tsv"), sep="\t", index=False)
    pd.DataFrame({"TF": list(sel_oligo.keys()), "motif_id": list(sel_oligo.values())}).to_csv(
        os.path.join(OUT_DIR, "selected_motifs_oligo.tsv"), sep="\t", index=False)

    # -------------------------------------------------------------------
    #  Flexible cell type matching (no hardcoding)
    # -------------------------------------------------------------------
    print("\n" + "="*60)
    print("Identifying cell types in dataset...")
    print("="*60)

    available_cell_types = list(dev.obs['cell_type'].unique())
    print(f"Available cell types: {available_cell_types}")

    # Define search terms for each glia type
    glia_search_config = {
        'Astrocytes': {
            'search_terms': ['astro', 'astrocyte'],
            'tf_selection': sel_astro
        },
        'Microglia': {
            'search_terms': ['micro', 'microglia', 'pvm'],
            'tf_selection': sel_micro
        },
        'Oligodendrocytes': {
            'search_terms': ['oligo', 'oligodendrocyte'],
            'tf_selection': sel_oligo
        }
    }

    # Build dynamic glia_targets list
    glia_targets = []
    for label, config in glia_search_config.items():
        matched_ct = find_cell_type(dev, config['search_terms'])
        
        if matched_ct:
            print(f"  ✓ Found {label}: '{matched_ct}'")
            glia_targets.append((matched_ct, config['tf_selection']))
        else:
            print(f"  ⚠ No match for {label} (searched: {config['search_terms']})")

    if not glia_targets:
        print("\n  WARNING: No glial cell types matched. Available types:")
        for ct in available_cell_types:
            print(f"    - {ct}")
        print("\n  Skipping glial differential analysis.")
        return
    else:
        print(f"\n  ✓ Proceeding with {len(glia_targets)} matched cell types")

    # -------------------------------------------------------------------
    #  Auto-detect conditions (no hardcoded group names)
    # -------------------------------------------------------------------
    print("\n" + "="*60)
    print("Setting up pairwise comparisons...")
    print("="*60)

    available_conditions = sorted(dev.obs["Condition"].dropna().unique())
    print(f"Available conditions: {available_conditions}")

    if len(available_conditions) < 2:
        raise ValueError(
            f"Need at least 2 conditions for differential analysis.\n"
            f"Found: {available_conditions}"
        )

    # Generate all pairwise comparisons
    comparisons = list(combinations(available_conditions, 2))
    print(f"\n✓ Running {len(comparisons)} pairwise comparisons:")
    for cond1, cond2 in comparisons:
        print(f"  - {cond2} vs {cond1}")

    # -------------------------------------------------------------------
    # Run differential tests for each cell type × comparison
    # -------------------------------------------------------------------
    all_results = []

    for ct, sel in glia_targets:
        if len(sel) == 0:
            print(f"\nNo selected motifs found for {ct}, skipping.")
            continue
        
        motif_ids = list(sel.values())
        
        # Run comparison for each condition pair
        for group2, group1 in comparisons:
            print(f"\n{'='*60}")
            print(f"Cell type: {ct} | Comparison: {group2} vs {group1}")
            print(f"{'='*60}")
            
            try:
                res_ct, sub_ct = diff_within_celltype(
                    dev, ct, "Condition", group1, group2, motif_ids, n_iter=1000
                )
            except Exception as e:
                print(f"  ✗ Error: {e}")
                continue
            
            # Save raw results
            out_tsv_raw = os.path.join(OUT_DIR, f"{ct}_diffmotifs_{group2}_vs_{group1}.tsv")
            res_ct.sort_values("pval").to_csv(out_tsv_raw, sep="\t")
            print(f"  ✓ Saved: {out_tsv_raw}")

            # Process results for plotting
            res_trim = res_ct.loc[res_ct.index.intersection(motif_ids)].copy()
            if res_trim.index.duplicated().any():
                n_dup = int(res_trim.index.duplicated().sum())
                print(f"  Warning: {n_dup} duplicate motif IDs; keeping first occurrence.")
                res_trim = res_trim[~res_trim.index.duplicated(keep="first")]
            
            if "qval" not in res_trim.columns:
                res_trim["qval"] = multipletests(res_trim["pval"].values, method="fdr_bh")[1]

            order = res_trim.sort_values("pval").index.tolist()
            labels = [mid.split("\t")[-1] for mid in order]
            delta  = np.asarray(res_trim.loc[order, "stat"].values).reshape(-1)
            qvals  = np.asarray(res_trim.loc[order, "qval"].values).reshape(-1)

            assert len(labels) == len(delta) == len(qvals)

            plot_df = pd.DataFrame({"TF": labels, "delta": delta, "qval": qvals})

            # Bar plot
            plt.figure(figsize=(8, max(3, 0.4 * len(plot_df))))
            sns.barplot(data=plot_df, x="delta", y="TF", color="steelblue")
            plt.axvline(0, color="k", linewidth=1)
            plt.xlabel(f"Mean chromVAR z ({group2} − {group1})")
            plt.ylabel("Motif (TF)")
            plt.title(f"{ct}: differential motif activity")
            for i, row in plot_df.iterrows():
                if row["qval"] < 0.05:
                    x = row["delta"]
                    plt.text(x + (0.02 if x >= 0 else -0.02), i, "*", va="center",
                             ha="left" if x >= 0 else "right", fontsize=12)
            bar_pdf = os.path.join(OUT_DIR, f"{ct}_diffmotifs_bar_{group2}_vs_{group1}.pdf")
            plt.tight_layout(); plt.savefig(bar_pdf); plt.close()
            print(f"  ✓ Saved: {bar_pdf}")

            # Heatmap of per-group means
            df_ct = pd.DataFrame(np.asarray(sub_ct.X), index=sub_ct.obs_names, columns=sub_ct.var_names)
            df_ct["Condition"] = sub_ct.obs["Condition"].astype(str).values
            mean_g1 = df_ct[df_ct["Condition"] == group1].drop(columns=["Condition"]).mean(axis=0)
            mean_g2 = df_ct[df_ct["Condition"] == group2].drop(columns=["Condition"]).mean(axis=0)
            mean_mat = pd.DataFrame({group1: mean_g1, group2: mean_g2})
            mean_mat = mean_mat.loc[order]

            plt.figure(figsize=(4.8, max(3, 0.3 * len(order))))
            sns.heatmap(mean_mat, cmap="vlag", center=0, cbar_kws={"label": "mean z"})
            plt.title(f"{ct}: mean chromVAR z by group")
            plt.ylabel("Motif (TF)")
            heat_pdf = os.path.join(OUT_DIR, f"{ct}_group_means_heatmap_{group2}_vs_{group1}.pdf")
            plt.tight_layout(); plt.savefig(heat_pdf); plt.close()
            print(f"  ✓ Saved: {heat_pdf}")

            # Save trimmed results
            out_tsv = os.path.join(OUT_DIR, f"{ct}_diffmotifs_{group2}_vs_{group1}.trimmed.tsv")
            res_trim.sort_values("pval").to_csv(out_tsv, sep="\t")
            print(f"  ✓ Saved: {out_tsv}")

            res_trim["cell_type"] = ct
            res_trim["comparison"] = f"{group2}_vs_{group1}"
            all_results.append(res_trim)

    # Summary heatmap across glia
    if all_results:
        df_all = pd.concat(all_results, axis=0)
        
        # Create one heatmap per comparison
        for group2, group1 in comparisons:
            df_comp = df_all[df_all["comparison"] == f"{group2}_vs_{group1}"]
            if df_comp.empty:
                continue
            
            mat = df_comp.pivot_table(index="cell_type", columns="TF", values="stat", aggfunc="first")
            
            plt.figure(figsize=(0.5*mat.shape[1] + 3, 0.7*mat.shape[0] + 2))
            sns.heatmap(mat, cmap="vlag", center=0, cbar_kws={"label": f"{group2} − {group1} mean z"})
            plt.title(f"Differential motif activity across glia ({group2} vs {group1})")
            glia_pdf = os.path.join(OUT_DIR, f"Glia_delta_heatmap_{group2}_vs_{group1}.pdf")
            plt.tight_layout(); plt.savefig(glia_pdf); plt.close()
            print(f"\n✓ Saved summary: {glia_pdf}")
    else:
        print("\nNo glial differential results produced.")

    print("\n" + "="*60)
    print("vis_chromvar_nf.py finished successfully.")
    print("="*60)

if __name__ == "__main__":
    main()
