#!/usr/bin/env python3
"""
tf_chip_corroboration.py  —  TF_CHIP_CORROBORATION (Gap B closer vs MAESTRO)

Closes the one genuine methodological gap between FORGE and MAESTRO: FORGE's TF
calls are motif- + footprint-based (chromVAR, SCENIC+ motifs, scPrinter), whereas
MAESTRO grounds its TF ranking in *empirical* ChIP-seq occupancy from CistromeDB
via GIGGLE/LISA.

This module adds that empirical layer as an orthogonal CORROBORATION input: for
each FORGE-nominated (cell_type, TF), it tests whether the TF's external ChIP-seq
peaks (ChIP-Atlas / Cistrome, one BED per TF) are enriched within that cell type's
characteristic ATAC peaks — a GIGGLE/LOLA-style Fisher's-exact overlap test over
the consensus-peak universe.

The result is a *triangulation* table: motif (chromVAR z) x empirical (ChIP OR/p)
[x footprint, if provided]. This is a strictly stronger TF claim than MAESTRO's
ChIP-only ranking, because FORGE keeps the motif/footprint evidence too.

Nominations are consumed from EXTRACT_CHROMVAR_MOTIFS's per_celltype_motifs.json
(canonical per-CT TF lists); optionally augmented with SCENIC+ eRegulon TFs.
Query peaks are derived from the annotated peak matrix by one-vs-rest accessibility
(no dependency on differential-peak calls being present).

Dependencies: anndata, numpy, pandas, scipy  (all in the `scgpu` container).
No pyranges / bedtools / statsmodels — interval overlap is a sorted-sweep, BH is
implemented inline — to keep the container contract minimal.
"""
import argparse
import gzip
import json
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------------
# Peak-coordinate parsing
# ---------------------------------------------------------------------------
def parse_peak_coords(adata):
    """Return (chrom[], start[], end[]) arrays aligned to adata.var order.

    Tries structured .var columns first (chrom/seqnames + start + end), then
    falls back to parsing var_names of the form chr:start-end or chr-start-end.
    """
    var = adata.var
    cols = {c.lower(): c for c in var.columns}

    def _find(*cands):
        for cand in cands:
            if cand in cols:
                return cols[cand]
        return None

    c_chrom = _find("chrom", "chr", "seqnames", "seqname", "chromosome")
    c_start = _find("start", "chromstart", "peak_start")
    c_end = _find("end", "chromend", "peak_end")

    if c_chrom and c_start and c_end:
        chrom = var[c_chrom].astype(str).values
        start = var[c_start].astype(np.int64).values
        end = var[c_end].astype(np.int64).values
        return chrom, start, end

    # Fall back to var_names.  Accept "chr1:100-200" and "chr1-100-200".
    chrom, start, end = [], [], []
    for name in adata.var_names:
        s = str(name)
        if ":" in s:
            c, rng = s.split(":", 1)
            a, b = rng.split("-")
        else:
            parts = s.rsplit("-", 2)
            if len(parts) != 3:
                raise ValueError(f"Cannot parse peak coordinate from var_name: {name!r}")
            c, a, b = parts
        chrom.append(c)
        start.append(int(a))
        end.append(int(b))
    return np.array(chrom, dtype=object), np.array(start, dtype=np.int64), np.array(end, dtype=np.int64)


def normalize_chrom(c):
    """Normalize chromosome names so 'chr1' and '1' compare equal."""
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


# ---------------------------------------------------------------------------
# ChIP-reference loading + overlap
# ---------------------------------------------------------------------------
def _open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def index_chip_ref_dir(chip_ref_dir):
    """Map lowercase TF symbol -> BED path.  Filenames like SPI1.bed(.gz)."""
    idx = {}
    if not os.path.isdir(chip_ref_dir):
        raise NotADirectoryError(f"--chip-ref-dir not a directory: {chip_ref_dir}")
    for fn in os.listdir(chip_ref_dir):
        low = fn.lower()
        if low.endswith(".bed") or low.endswith(".bed.gz"):
            tf = fn.split(".")[0]
            idx[tf.lower()] = os.path.join(chip_ref_dir, fn)
    return idx


def load_chip_intervals(bed_path):
    """Load a ChIP BED into {chrom_norm: (starts_sorted, ends_sorted)} arrays."""
    by_chrom = defaultdict(list)
    with _open_maybe_gzip(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            by_chrom[normalize_chrom(f[0])].append((int(f[1]), int(f[2])))
    out = {}
    for c, ivals in by_chrom.items():
        ivals.sort()
        starts = np.fromiter((s for s, _ in ivals), dtype=np.int64, count=len(ivals))
        ends = np.fromiter((e for _, e in ivals), dtype=np.int64, count=len(ivals))
        out[c] = (starts, ends)
    return out


def peaks_overlapping_chip(peak_chrom_norm, peak_start, peak_end, chip):
    """Boolean mask over peaks: True if peak overlaps >=1 ChIP interval.

    Sorted-sweep: for each peak, binary-search ChIP starts on its chromosome and
    scan the small candidate window. O(n_peaks * log n_chip) — no pyranges needed.
    """
    mask = np.zeros(len(peak_start), dtype=bool)
    # Group peak indices by chromosome to reuse each chrom's sorted arrays.
    order = defaultdict(list)
    for i, c in enumerate(peak_chrom_norm):
        order[c].append(i)
    for c, idxs in order.items():
        if c not in chip:
            continue
        cstarts, cends = chip[c]
        # A peak [ps,pe) overlaps some chip [cs,ce) iff there exists an interval
        # with cs < pe and ce > ps.  Find chip intervals whose start < pe via
        # searchsorted, then confirm ce > ps within a bounded back-scan.
        cmax_end = np.maximum.accumulate(cends)  # running max end, monotone
        for i in idxs:
            ps, pe = peak_start[i], peak_end[i]
            hi = np.searchsorted(cstarts, pe, side="left")  # chip starts < pe
            if hi == 0:
                continue
            # Among chip[:hi], overlap iff any end > ps. Running-max end lets us
            # test the whole prefix in O(1); if the prefix max end <= ps, none.
            if cmax_end[hi - 1] > ps:
                mask[i] = True
    return mask


# ---------------------------------------------------------------------------
# Per-cell-type characteristic peaks (one-vs-rest accessibility)
# ---------------------------------------------------------------------------
def characteristic_peaks(adata, ct_mask, top_k):
    """Indices of the top_k peaks most specific to the CT (in-vs-rest fraction)."""
    X = adata.X
    binX = (X > 0)
    if hasattr(binX, "toarray"):
        in_frac = np.asarray(binX[ct_mask].mean(axis=0)).ravel()
        rest_frac = np.asarray(binX[~ct_mask].mean(axis=0)).ravel()
    else:
        binX = np.asarray(binX)
        in_frac = binX[ct_mask].mean(axis=0)
        rest_frac = binX[~ct_mask].mean(axis=0)
    spec = in_frac - rest_frac
    k = min(top_k, spec.shape[0])
    return np.argpartition(-spec, k - 1)[:k]


# ---------------------------------------------------------------------------
# Benjamini-Hochberg
# ---------------------------------------------------------------------------
def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


# ---------------------------------------------------------------------------
def load_nominations(motif_json, eregulon_json):
    """Return {cell_type: {tf: chromvar_zscore_or_nan}} and a global TF set."""
    with open(motif_json) as f:
        mj = json.load(f)
    noms = defaultdict(dict)
    for ct, d in mj.get("cell_types", {}).items():
        zmap = {det["tf"]: det.get("mean_abs_zscore", np.nan)
                for det in d.get("motif_details", [])}
        for tf in d.get("top_tfs", []):
            noms[ct][tf] = zmap.get(tf, np.nan)

    if eregulon_json and os.path.exists(eregulon_json):
        with open(eregulon_json) as f:
            ej = json.load(f)
        # Flexible: accept {cell_type: [tf,...]} or {tf: {...}} or {"eregulons":[...]}
        if isinstance(ej, dict) and any(isinstance(v, list) for v in ej.values()):
            for ct, tfs in ej.items():
                if isinstance(tfs, list):
                    for tf in tfs:
                        noms.setdefault(ct, {}).setdefault(str(tf), np.nan)
        else:
            # No CT structure: add TFs to every nominated CT as candidates.
            tf_keys = list(ej.keys()) if isinstance(ej, dict) else list(ej)
            for ct in list(noms.keys()):
                for tf in tf_keys:
                    noms[ct].setdefault(str(tf), np.nan)
    return noms


def resolve_ct_col(adata, requested):
    if requested and requested in adata.obs.columns:
        return requested
    for cand in ("cell_type", "broad_cell_type", "celltypist_prediction",
                 "cell_type_broad", "predicted_labels"):
        if cand in adata.obs.columns:
            return cand
    raise ValueError(
        f"No usable cell-type column (requested={requested!r}); "
        f"available: {list(adata.obs.columns)}")


def main():
    ap = argparse.ArgumentParser(description="Empirical ChIP corroboration of FORGE TF calls")
    ap.add_argument("--motif-json", required=True, help="per_celltype_motifs.json (nominations)")
    ap.add_argument("--peak-matrix", required=True, help="annotated peak matrix h5ad (query source)")
    ap.add_argument("--chip-ref-dir", required=True, help="dir of per-TF ChIP BEDs (<TF>.bed[.gz])")
    ap.add_argument("--eregulon-json", default=None, help="optional SCENIC+ eRegulon TF list")
    ap.add_argument("--cell-type-col", default=None, help="obs column for cell type")
    ap.add_argument("--top-peaks-per-ct", type=int, default=2000)
    ap.add_argument("--min-cells", type=int, default=50, help="resolution floor: min cells")
    ap.add_argument("--min-pct", type=float, default=0.01, help="resolution floor: min fraction")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-odds-ratio", type=float, default=1.5)
    ap.add_argument("--min-overlap", type=int, default=3, help="min query peaks overlapping ChIP")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("TF_CHIP_CORROBORATION — empirical ChIP corroboration of FORGE TF calls", flush=True)
    print(f"  peak matrix: {args.peak_matrix}", flush=True)
    print(f"  chip ref:    {args.chip_ref_dir}", flush=True)

    adata = ad.read_h5ad(args.peak_matrix)
    ct_col = resolve_ct_col(adata, args.cell_type_col)
    print(f"  cell-type column: {ct_col}", flush=True)

    chrom, start, end = parse_peak_coords(adata)
    chrom_norm = np.array([normalize_chrom(c) for c in chrom], dtype=object)
    N = len(start)
    print(f"  peaks: {N}  cells: {adata.n_obs}", flush=True)

    noms = load_nominations(args.motif_json, args.eregulon_json)
    chip_idx = index_chip_ref_dir(args.chip_ref_dir)
    print(f"  ChIP refs available for {len(chip_idx)} TFs", flush=True)

    total_cells = adata.n_obs
    res_floor = max(args.min_cells, int(args.min_pct * total_cells))

    # Cache per-TF ChIP-overlap masks over the universe (reused across CTs).
    tf_universe_mask = {}
    missing_refs = set()

    def universe_mask_for(tf):
        key = tf.lower()
        if key in tf_universe_mask:
            return tf_universe_mask[key]
        if key not in chip_idx:
            missing_refs.add(tf)
            tf_universe_mask[key] = None
            return None
        chip = load_chip_intervals(chip_idx[key])
        m = peaks_overlapping_chip(chrom_norm, start, end, chip)
        tf_universe_mask[key] = m
        return m

    rows = []
    skipped_cts = []
    for ct in sorted(noms.keys()):
        ct_mask = (adata.obs[ct_col].astype(str) == str(ct)).values
        n_ct = int(ct_mask.sum())
        if n_ct < res_floor:
            skipped_cts.append((ct, n_ct, res_floor))
            print(f"  SKIP {ct}: n_cells={n_ct} < floor={res_floor}", flush=True)
            continue

        q_idx = characteristic_peaks(adata, ct_mask, args.top_peaks_per_ct)
        q_mask = np.zeros(N, dtype=bool)
        q_mask[q_idx] = True
        n_q = int(q_mask.sum())

        for tf, zscore in noms[ct].items():
            umask = universe_mask_for(tf)
            if umask is None:
                continue  # no ChIP ref — recorded in missing_refs
            a = int(np.sum(q_mask & umask))          # query & chip
            m = int(np.sum(umask))                   # universe & chip
            b = n_q - a                              # query, no chip
            c = m - a                                # non-query & chip
            d = (N - n_q) - c                        # non-query, no chip
            if a + b == 0 or a + c == 0:
                continue
            odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            expected = n_q * (m / N) if N else 0.0
            log2_or = float(np.log2(odds)) if (odds and np.isfinite(odds) and odds > 0) else np.nan
            rows.append({
                "cell_type": ct,
                "tf": tf,
                "n_ct_cells": n_ct,
                "chromvar_zscore": zscore,
                "n_query_peaks": n_q,
                "query_overlap_chip": a,
                "expected_overlap": round(expected, 2),
                "universe_overlap_chip": m,
                "odds_ratio": round(float(odds), 4) if np.isfinite(odds) else np.inf,
                "log2_odds_ratio": round(log2_or, 4) if np.isfinite(log2_or) else np.nan,
                "p_value": p,
            })

    df = pd.DataFrame(rows)
    if len(df):
        df["q_value"] = bh_fdr(df["p_value"].values)
        df["empirical_support"] = (
            (df["q_value"] < args.fdr)
            & (df["odds_ratio"] >= args.min_odds_ratio)
            & (df["query_overlap_chip"] >= args.min_overlap)
        )
        # Evidence-layer triangulation label.
        def layers(r):
            L = []
            if pd.notna(r["chromvar_zscore"]):
                L.append("motif")
            if r["empirical_support"]:
                L.append("chip")
            return "+".join(L) if L else "none"
        df["evidence_layers"] = df.apply(layers, axis=1)
        df = df.sort_values(["cell_type", "empirical_support", "odds_ratio"],
                            ascending=[True, False, False])
    else:
        df = pd.DataFrame(columns=[
            "cell_type", "tf", "n_ct_cells", "chromvar_zscore", "n_query_peaks",
            "query_overlap_chip", "expected_overlap", "universe_overlap_chip",
            "odds_ratio", "log2_odds_ratio", "p_value", "q_value",
            "empirical_support", "evidence_layers"])

    corr_path = os.path.join(args.outdir, "tf_chip_corroboration.tsv")
    df.to_csv(corr_path, sep="\t", index=False)

    # Triangulation view: one row per (CT, TF) with the compact evidence summary.
    tri = df[["cell_type", "tf", "chromvar_zscore", "odds_ratio", "q_value",
              "empirical_support", "evidence_layers"]].copy()
    tri_path = os.path.join(args.outdir, "tf_triangulation.tsv")
    tri.to_csv(tri_path, sep="\t", index=False)

    # Missing-reference log — SILENT CAPS ARE FORBIDDEN (feedback_no-silent-caps).
    miss_path = os.path.join(args.outdir, "missing_chip_refs.txt")
    with open(miss_path, "w") as f:
        f.write("TFs nominated by FORGE with NO ChIP reference BED available\n")
        f.write("(these were NOT corroborated — absence of evidence, not evidence of absence)\n")
        f.write("=" * 60 + "\n")
        for tf in sorted(missing_refs):
            f.write(tf + "\n")

    n_support = int(df["empirical_support"].sum()) if len(df) else 0
    summ_path = os.path.join(args.outdir, "tf_chip_corroboration_summary.txt")
    with open(summ_path, "w") as f:
        f.write("TF_CHIP_CORROBORATION summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Cell types tested:            {df['cell_type'].nunique() if len(df) else 0}\n")
        f.write(f"(CT,TF) tests run:            {len(df)}\n")
        f.write(f"Empirically corroborated:     {n_support}\n")
        f.write(f"TFs with no ChIP reference:   {len(missing_refs)}\n")
        f.write(f"Cell types skipped (floor):   {len(skipped_cts)}\n")
        f.write(f"Resolution floor (cells):     {res_floor}\n")
        if skipped_cts:
            f.write("\nSkipped CTs:\n")
            for ct, n, fl in skipped_cts:
                f.write(f"  {ct}: {n} < {fl}\n")

    print(f"  wrote {corr_path}  ({len(df)} tests, {n_support} corroborated)", flush=True)
    print(f"  wrote {tri_path}", flush=True)
    print(f"  wrote {miss_path}  ({len(missing_refs)} TFs without refs)", flush=True)
    print("TF_CHIP_CORROBORATION complete", flush=True)


if __name__ == "__main__":
    sys.exit(main())
