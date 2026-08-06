"""
rank_enhancer_strip_genes.py

Ranks per-(CT, TF) target genes for MSFP enhancer strip visualization using:
  1. scPrinter binding scores from enhancer_tfbs_*.h5ad obsm profiles
  2. Cicero peak-to-peak co-accessibility connections

For each (CT, TF):
  - Select top-N regions by max binding score (position-wise max in obsm profile)
  - Find co-accessible Cicero peaks for those regions (overlap-based matching)
  - Map connected peaks to gene promoters (TSS ±promoter_window bp from GTF)
  - Score genes: sum(binding_score × coaccess_score) across all linked regions
  - Take top-K genes per (CT, TF)

Outputs (all in working directory):
  per_ct_genes.csv          — cell_type, strip_target_genes (comma-sep union per CT)
  per_ct_tf_genes.json      — {ct: {tf: [gene, ...]}} for RENDER wiring
  strip_gene_ranking.json   — full scored ranking for diagnostics
"""
import argparse, gzip, json, os, re, sys
from bisect import bisect_left, bisect_right
from collections import defaultdict

import anndata as ad
import numpy as np
import pandas as pd


# ── helpers ───────────────────────────────────────────────────────────────────

def parse_locus(s):
    m = re.match(r'(chr[^:]+):(\d+)-(\d+)', s)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), int(m.group(3))


def build_cicero_index(connections_path, min_coaccess=0.05):
    """Return {peak_id: [(connected_peak_id, coaccess), ...]} and sorted interval lists."""
    print(f"[rank] Loading Cicero connections from {connections_path}...", flush=True)
    df = pd.read_csv(connections_path, sep='\t',
                     compression='gzip' if connections_path.endswith('.gz') else None)
    df = df[df['coaccess'] >= min_coaccess]

    conn_dict = defaultdict(list)
    for _, row in df.iterrows():
        p1, p2, ca = row['Peak1'], row['Peak2'], float(row['coaccess'])
        conn_dict[p1].append((p2, ca))
        conn_dict[p2].append((p1, ca))

    # Build sorted interval structure for overlap queries: chrom -> [(start, end, peak_id)]
    peaks_by_chrom = defaultdict(list)
    all_peaks = set(df['Peak1'].tolist() + df['Peak2'].tolist())
    for pk in all_peaks:
        chrom, s, e = parse_locus(pk)
        if chrom:
            peaks_by_chrom[chrom].append((s, e, pk))
    for chrom in peaks_by_chrom:
        peaks_by_chrom[chrom].sort()

    print(f"[rank]   {len(all_peaks):,} unique Cicero peaks, {len(df):,} connections (coaccess≥{min_coaccess})", flush=True)
    return conn_dict, peaks_by_chrom


def find_overlapping_cicero_peaks(chrom, q_start, q_end, peaks_by_chrom):
    """Find Cicero peaks that overlap [q_start, q_end]."""
    entries = peaks_by_chrom.get(chrom, [])
    if not entries:
        return []
    starts = [e[0] for e in entries]
    # All peaks with start <= q_end
    hi = bisect_right(starts, q_end)
    # Filter those where end >= q_start
    return [(e[1], e[2]) for e in entries[:hi] if e[1] >= q_start]


def build_gene_tss_index(gtf_path, promoter_window=2000):
    """Return {chrom: sorted [(tss_start, tss_end, gene_name)]} for promoter windows."""
    print(f"[rank] Loading GTF from {gtf_path}...", flush=True)
    tss_by_chrom = defaultdict(list)
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            chrom = parts[0]
            start, end, strand = int(parts[3]), int(parts[4]), parts[6]
            m = re.search(r'gene_name "([^"]+)"', parts[8])
            if not m:
                continue
            gene = m.group(1)
            tss = start if strand == '+' else end
            tss_by_chrom[chrom].append((tss - promoter_window, tss + promoter_window, gene))

    for chrom in tss_by_chrom:
        tss_by_chrom[chrom].sort()

    n_genes = sum(len(v) for v in tss_by_chrom.values())
    print(f"[rank]   {n_genes:,} genes loaded", flush=True)
    return tss_by_chrom


def find_genes_for_peak(peak_id, tss_by_chrom):
    """Return list of gene names whose promoter window overlaps this peak."""
    chrom, s, e = parse_locus(peak_id)
    if not chrom:
        return []
    entries = tss_by_chrom.get(chrom, [])
    if not entries:
        return []
    starts = [x[0] for x in entries]
    # promoter windows with start <= peak_end
    hi = bisect_right(starts, e)
    # filter: window_end >= peak_start
    return [x[2] for x in entries[:hi] if x[1] >= s]


# ── per-(CT, TF) ranking ──────────────────────────────────────────────────────

def rank_genes_for_ct_tf(h5ad_path, ct_name, conn_dict, peaks_by_chrom,
                          tss_by_chrom, top_n_regions, top_k_genes):
    """
    Returns list of (gene, score) tuples sorted descending for this (CT, TF).
    """
    a = ad.read_h5ad(h5ad_path)
    n_obsm = len(a.obsm)
    if n_obsm == 0:
        return []

    # Identify focal-CT row: use obs_names, fall back to row 0
    if ct_name in a.obs_names:
        ct_idx = list(a.obs_names).index(ct_name)
    else:
        ct_idx = 0

    # Score each region: max binding probability for focal CT across position bins
    region_scores = {}
    for region_key in a.obsm.keys():
        arr = np.asarray(a.obsm[region_key])
        if arr.ndim == 1:
            score = float(arr.max())
        elif arr.ndim == 2:
            # shape (n_ct, n_bins) or (n_conditions, n_bins)
            if arr.shape[0] == 1:
                score = float(arr[0].max())
            elif ct_idx < arr.shape[0]:
                score = float(arr[ct_idx].max())
            else:
                score = float(arr.max())
        else:
            score = float(arr.max())
        region_scores[region_key] = score

    # Select top-N by binding score
    top_regions = sorted(region_scores.items(), key=lambda x: -x[1])[:top_n_regions]

    # For each top region: find Cicero connections → gene promoters
    gene_scores = defaultdict(float)
    for region_key, binding_score in top_regions:
        chrom, rs, re_ = parse_locus(region_key)
        if not chrom:
            continue

        # Find Cicero peaks that overlap this enhancer region
        overlapping = find_overlapping_cicero_peaks(chrom, rs, re_, peaks_by_chrom)

        for _cic_end, cic_peak_id in overlapping:
            # Find all peaks co-accessible with this Cicero peak
            for connected_peak, coaccess in conn_dict.get(cic_peak_id, []):
                # Map connected peak to gene promoters
                for gene in find_genes_for_peak(connected_peak, tss_by_chrom):
                    gene_scores[gene] += binding_score * coaccess

    ranked = sorted(gene_scores.items(), key=lambda x: -x[1])[:top_k_genes]
    return ranked


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tfbs-h5ads',      nargs='+', required=True,
                   help='enhancer_tfbs_*.h5ad files (all TFs for all CTs)')
    p.add_argument('--cicero-connections', required=True,
                   help='cicero_connections.tsv.gz')
    p.add_argument('--gtf',             required=True)
    p.add_argument('--top-n-regions',   type=int, default=100)
    p.add_argument('--top-k-genes',     type=int, default=5)
    p.add_argument('--min-coaccess',    type=float, default=0.05)
    p.add_argument('--promoter-window', type=int, default=2000,
                   help='bp upstream/downstream of TSS for promoter definition')
    args = p.parse_args()

    # Parse CT and TF from h5ad filenames and obs_names
    # Filename convention: enhancer_tfbs_{safe_ct}_{tf}.h5ad
    # safe_ct = ct.replace('/', '_') — spaces preserved
    h5ad_meta = []
    for path in args.tfbs_h5ads:
        fname = os.path.basename(path)
        stem = fname.replace('enhancer_tfbs_', '').replace('.h5ad', '')
        # Try uns first; fall back to obs_names for CT, filename split for TF
        try:
            a_meta = ad.read_h5ad(path, backed='r')
            ct_obs = list(a_meta.obs_names)[0] if len(a_meta.obs_names) > 0 else None
            tf_uns = a_meta.uns.get('tf', None)
            a_meta.file.close()
        except Exception:
            ct_obs, tf_uns = None, None

        ct = ct_obs or stem  # fallback
        if tf_uns:
            tf = tf_uns
        else:
            # Use CT obs name to strip it from stem
            ct_safe = (ct or '').replace('/', '_')
            if stem.startswith(ct_safe + '_'):
                tf = stem[len(ct_safe) + 1:]
            else:
                tf = stem  # last resort
        h5ad_meta.append({'ct': ct, 'tf': tf, 'path': path})
        print(f"[rank] h5ad: CT={ct!r}  TF={tf!r}  file={fname}", flush=True)

    print(f"[rank] Total h5ads to process: {len(h5ad_meta)}", flush=True)

    # Build shared indexes (done once)
    conn_dict, peaks_by_chrom = build_cicero_index(args.cicero_connections, args.min_coaccess)
    tss_by_chrom = build_gene_tss_index(args.gtf, args.promoter_window)

    # Rank genes per (CT, TF)
    full_ranking = {}   # {ct: {tf: [(gene, score), ...]}}
    for meta in h5ad_meta:
        ct, tf, path = meta['ct'], meta['tf'], meta['path']
        print(f"[rank] Processing CT={ct!r} TF={tf!r}...", flush=True)
        ranked = rank_genes_for_ct_tf(
            path, ct, conn_dict, peaks_by_chrom, tss_by_chrom,
            args.top_n_regions, args.top_k_genes
        )
        if ct not in full_ranking:
            full_ranking[ct] = {}
        full_ranking[ct][tf] = ranked
        print(f"[rank]   top genes: {[(g, round(s, 4)) for g, s in ranked]}", flush=True)

    # Build per-CT gene union and per-(CT, TF) gene list
    per_ct_genes = {}     # {ct: set of genes}
    per_ct_tf_genes = {}  # {ct: {tf: [genes]}}

    for ct, tf_map in full_ranking.items():
        per_ct_tf_genes[ct] = {}
        per_ct_genes[ct] = set()
        for tf, ranked in tf_map.items():
            genes = [g for g, _ in ranked]
            per_ct_tf_genes[ct][tf] = genes
            per_ct_genes[ct].update(genes)

    # Write outputs
    # 1. per_ct_genes.csv — for ENHANCER_FOOTPRINTING_PER_CT_STRIP --strip-target-genes
    rows = []
    for ct, genes in per_ct_genes.items():
        rows.append({'cell_type': ct, 'strip_target_genes': ','.join(sorted(genes))})
    pd.DataFrame(rows).to_csv('per_ct_genes.csv', index=False)
    print(f"[rank] Wrote per_ct_genes.csv ({len(rows)} CTs)", flush=True)

    # 2. per_ct_tf_genes.json — for RENDER_MSFP_ENHANCER_STRIP channel wiring
    with open('per_ct_tf_genes.json', 'w') as f:
        json.dump(per_ct_tf_genes, f, indent=2)
    print(f"[rank] Wrote per_ct_tf_genes.json", flush=True)

    # 3. strip_gene_ranking.json — full scored ranking for diagnostics
    ranking_serializable = {
        ct: {tf: [{'gene': g, 'score': round(s, 6)} for g, s in ranked]
             for tf, ranked in tf_map.items()}
        for ct, tf_map in full_ranking.items()
    }
    with open('strip_gene_ranking.json', 'w') as f:
        json.dump(ranking_serializable, f, indent=2)
    print(f"[rank] Wrote strip_gene_ranking.json", flush=True)

    # Summary
    total_ct_tf = sum(len(tf_map) for tf_map in full_ranking.values())
    total_genes = sum(len(g) for g in per_ct_genes.values())
    print(f"\n[rank] Done. {len(full_ranking)} CTs, {total_ct_tf} CT×TF pairs, "
          f"{total_genes} total per-CT gene assignments", flush=True)


if __name__ == '__main__':
    main()
