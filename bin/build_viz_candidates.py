#!/usr/bin/env python3
"""
build_viz_candidates.py — Swarup-chain (gene, TF, cell_type) candidate selector.

Replaces the 17,490-panel cartesian fan-out (genes × TFs) with a data-driven
filter mirroring Swarup PiD/AD lab's chain (Science Adv. 2024-2025; repo
swaruplabUCI/scMultiomics_identifies_shared_and_distinct_pathways_in_PiDandAD):

  gene-window enhancer ∩ TF motif site (in that ct)
  ∧ TF differential FDR < cutoff (in that ct)
  ∧ |log2FC| > cutoff
  → top-N TFs per (gene, ct) by |log2FC|

Inputs:
  --track-manifest       JSON from PREPARE_ENHANCER_VIZ_TRACKS (gene windows +
                         per-gene .ini)
  --motif-scan-manifest  JSON from MOTIF_SCAN_ENHANCERS (region_sets per ct)
  --motif-bed-dir        Directory holding *.bed files referenced by region_sets
  --diff-csv-dir         Directory holding tf_differential_<ct>_<trt>_vs_<ctrl>.csv
                         (only the .csv files; .summary stubs from skipped ct's
                         are ignored)
  --top-n-per-gene-per-ct Top-N per (gene, ct) by |log2FC| (default 3)
  --fdr-cutoff           Max FDR per TF (default 0.05)
  --lfc-cutoff           Min |log2FC| per TF (default 0.5)
  --min-motif-sites      Min motif sites in gene-window enhancers (default 1)
  --out-csv              Output CSV path (default candidates.csv)

Output schema (candidates.csv):
  gene,tf,cell_type,n_motif_sites,log2fc,pval_adj,score

Where `score = |log2fc| * log10(n_motif_sites + 1)` provides a single ranking
column the consumer can use without re-deriving it. Consumer is main.nf via
`splitCsv(header: true)` → tuple(gene, tf, cell_type) fan-out.

Falls back to all-genes × top-N-TFs-by-overlap if no diff CSVs are provided
(e.g. condition not configured) — prints a [WARN] and logs which mode was used.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--track-manifest', required=True)
    p.add_argument('--motif-scan-manifest', required=True)
    p.add_argument('--motif-bed-dir', required=True,
                   help='Directory containing the .bed files referenced from '
                        'region_set_manifest (e.g. tf_enhancer_region_sets/)')
    p.add_argument('--diff-csv-dir', default='',
                   help='Directory with tf_differential_<ct>_<trt>_vs_<ctrl>.csv. '
                        'When empty, falls back to motif-overlap-only ranking.')
    p.add_argument('--treatment', default='',
                   help='Treatment condition label (used to robustly parse the '
                        'ct from filenames where trt contains underscores).')
    p.add_argument('--control', default='',
                   help='Control condition label.')
    p.add_argument('--top-n-per-gene-per-ct', type=int, default=3)
    p.add_argument('--fdr-cutoff', type=float, default=0.05)
    p.add_argument('--lfc-cutoff', type=float, default=0.5)
    p.add_argument('--min-motif-sites', type=int, default=1)
    p.add_argument('--out-csv', default='candidates.csv')
    p.add_argument('--debug-log', default='candidates.log')
    return p.parse_args()


def load_gene_windows(track_manifest_path):
    """Return {gene: (chrom, start, end)} from track_manifest.gene_regions."""
    t = json.load(open(track_manifest_path))
    out = {}
    for g, info in (t.get('gene_regions') or {}).items():
        if g == '_global':
            continue
        region = info.get('region', '')
        if not region or ':' not in region:
            continue
        c, se = region.split(':')
        s, e = se.split('-')
        out[g] = (c, int(s), int(e))
    return out


def load_motif_beds(region_set_manifest_path, motif_bed_dir):
    """Return {(cell_type, tf): [(chrom, start, end), ...]} indexed by chrom for fast overlap."""
    m = json.load(open(region_set_manifest_path))
    by_chrom = defaultdict(list)  # (ct, tf, chrom) -> [(start, end)]
    pair_total = 0
    pair_loaded = 0
    for r in m.get('region_sets', []):
        pair_total += 1
        bed_path = os.path.join(motif_bed_dir, r['bed_file'])
        if not os.path.exists(bed_path):
            continue
        ct = r['cell_type']
        tf = r['tf']
        with open(bed_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    try:
                        by_chrom[(ct, tf, parts[0])].append((int(parts[1]), int(parts[2])))
                    except ValueError:
                        continue
        pair_loaded += 1
    return by_chrom, pair_total, pair_loaded


def parse_ct_from_diff_filename(name, trt='', ctrl=''):
    """tf_differential_<ct>_<sanitize(trt)>_vs_<sanitize(ctrl)>.csv → (ct, trt, ctrl).

    The differential_tf_accessibility.py output sanitizes each component via
    re.sub(r"[^A-Za-z0-9_.-]", "_", x). Both `trt` and `ctrl` may contain
    underscores (e.g. trt='SREBF1_OE'), so a naive rsplit('_', 1) on the
    middle won't recover ct correctly. When trt/ctrl are passed in, this
    function strips the known suffix `_<sanitize(trt)>_vs_<sanitize(ctrl)>`
    explicitly. Without trt/ctrl, falls back to the legacy heuristic and may
    misparse names where trt has `_`.
    """
    def san(s):
        return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s)) if s else ''

    base = os.path.basename(name)
    if not base.startswith('tf_differential_') or not base.endswith('.csv'):
        return ('', '', '')
    middle = base[len('tf_differential_'):-len('.csv')]

    if trt and ctrl:
        suffix = f'_{san(trt)}_vs_{san(ctrl)}'
        if middle.endswith(suffix):
            return (middle[:-len(suffix)], trt, ctrl)
        return ('', '', '')

    # Legacy heuristic: assume single-token trt/ctrl
    if '_vs_' not in middle:
        return ('', '', '')
    pre_vs, ctrl_part = middle.rsplit('_vs_', 1)
    if '_' not in pre_vs:
        return ('', '', '')
    ct, trt_part = pre_vs.rsplit('_', 1)
    return (ct, trt_part, ctrl_part)


def load_diff_csvs(diff_csv_dir, treatment='', control=''):
    """Return {ct: {tf: (log2fc, pval_adj)}}; skips empty/.summary stubs.

    The ct key is the sanitized form (re.sub(r"[^A-Za-z0-9_.-]", "_", ct)),
    matching what differential_tf_accessibility.py wrote into the filename.
    Pass `treatment`/`control` to robustly parse filenames where trt/ctrl
    contain underscores (e.g. trt='SREBF1_OE'). Without them, parsing falls
    back to a legacy heuristic that may misclaim ct.
    """
    out = {}
    if not diff_csv_dir or not os.path.isdir(diff_csv_dir):
        return out
    for fname in os.listdir(diff_csv_dir):
        if not fname.startswith('tf_differential_') or not fname.endswith('.csv'):
            continue
        path = os.path.join(diff_csv_dir, fname)
        if os.path.getsize(path) < 200:
            continue  # likely a stub summary mistakenly named .csv
        ct_san, _trt, _ctrl = parse_ct_from_diff_filename(
            fname, trt=treatment, ctrl=control,
        )
        if not ct_san:
            continue
        per_tf = {}
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # 'motif' column is "<MA-id>\t<TF_name>" emitted by scanpy via a
                # tab-joined index; split to recover the TF name. If no tab,
                # treat the full string as the TF name.
                m = row.get('motif', '')
                if '\t' in m:
                    tf_name = m.split('\t', 1)[1].strip()
                else:
                    tf_name = m.strip()
                if not tf_name:
                    continue
                try:
                    log2fc = float(row.get('logfoldchange', 'nan'))
                    pval_adj = float(row.get('pval_adj', 'nan'))
                except ValueError:
                    continue
                if math.isnan(log2fc) or math.isnan(pval_adj):
                    continue
                # Keep most-significant entry per TF (in case of duplicate motif IDs).
                prev = per_tf.get(tf_name)
                if prev is None or pval_adj < prev[1]:
                    per_tf[tf_name] = (log2fc, pval_adj)
        if per_tf:
            out[ct_san] = per_tf
    return out


def count_overlaps(by_chrom, ct, tf, chrom, gs, ge):
    ivs = by_chrom.get((ct, tf, chrom), [])
    return sum(1 for s, e in ivs if e > gs and s < ge)


def main():
    args = parse_args()

    log_lines = []
    def log(msg):
        log_lines.append(msg)
        print(msg, flush=True)

    log(f"[build_viz_candidates] track_manifest={args.track_manifest}")
    log(f"[build_viz_candidates] motif_scan_manifest={args.motif_scan_manifest}")
    log(f"[build_viz_candidates] motif_bed_dir={args.motif_bed_dir}")
    log(f"[build_viz_candidates] diff_csv_dir={args.diff_csv_dir or '<none>'}")
    log(f"[build_viz_candidates] top_n={args.top_n_per_gene_per_ct}, "
        f"fdr<{args.fdr_cutoff}, |log2fc|>{args.lfc_cutoff}, "
        f"min_sites>={args.min_motif_sites}")

    # 1. Gene windows
    gene_window = load_gene_windows(args.track_manifest)
    log(f"[1] {len(gene_window)} gene windows loaded")

    # 2. TF motif sites per (ct, tf)
    by_chrom, n_total, n_loaded = load_motif_beds(
        args.motif_scan_manifest, args.motif_bed_dir
    )
    cell_types = sorted({k[0] for k in by_chrom.keys()})
    log(f"[2] {n_loaded}/{n_total} (ct, tf) bed files loaded across {len(cell_types)} cell types")

    # 3. Differential TF data per ct (optional)
    diff = (load_diff_csvs(args.diff_csv_dir, args.treatment, args.control)
            if args.diff_csv_dir else {})
    if diff:
        ct_with_diff = sorted(diff.keys())
        log(f"[3] differential CSVs available for {len(diff)} cell types: {ct_with_diff}")
        # Sanity: how many TFs per ct
        for ct in ct_with_diff[:5]:
            log(f"     {ct}: {len(diff[ct])} TFs")
    else:
        log("[3] No differential CSVs — falling back to motif-overlap-only ranking.")

    # 4. Build candidates (gene, tf, ct)
    candidates = []  # (gene, tf, ct, n_sites, log2fc, pval_adj)
    n_skip_no_diff = 0
    n_skip_thresh = 0
    n_skip_no_sites = 0
    for g, (gc, gs, ge) in gene_window.items():
        for ct in cell_types:
            # Sanitized form for diff dict lookup
            ct_san = re.sub(r"[^A-Za-z0-9_.-]", "_", ct)
            tfs_in_ct = sorted({tf for (cct, tf, _c) in by_chrom.keys() if cct == ct})
            for tf in tfs_in_ct:
                n_sites = count_overlaps(by_chrom, ct, tf, gc, gs, ge)
                if n_sites < args.min_motif_sites:
                    n_skip_no_sites += 1
                    continue
                if diff:
                    diff_ct = diff.get(ct_san) or diff.get(ct)
                    if diff_ct is None:
                        n_skip_no_diff += 1
                        # No diff data for this ct — keep with NaN to allow ranking by motif if user opts in
                        log2fc = float('nan')
                        pval_adj = float('nan')
                    else:
                        hit = diff_ct.get(tf)
                        if hit is None:
                            n_skip_no_diff += 1
                            continue
                        log2fc, pval_adj = hit
                        if pval_adj >= args.fdr_cutoff or abs(log2fc) <= args.lfc_cutoff:
                            n_skip_thresh += 1
                            continue
                else:
                    log2fc = float('nan')
                    pval_adj = float('nan')
                candidates.append((g, tf, ct, n_sites, log2fc, pval_adj))

    log(f"[4] pre-rank candidates: {len(candidates)} "
        f"(skipped: no_sites={n_skip_no_sites}, no_diff={n_skip_no_diff}, "
        f"below_threshold={n_skip_thresh})")

    # 5. Top-N per (gene, ct) by |log2fc| (fallback: by n_sites if log2fc is NaN)
    grouped = defaultdict(list)
    for c in candidates:
        grouped[(c[0], c[2])].append(c)

    final = []
    for key, group in grouped.items():
        def rank_key(t):
            _, _, _, n_sites, log2fc, _ = t
            if math.isnan(log2fc):
                # No-diff fallback: rank by motif site count
                return (0.0, n_sites)
            return (abs(log2fc), n_sites)
        group_sorted = sorted(group, key=rank_key, reverse=True)
        final.extend(group_sorted[:args.top_n_per_gene_per_ct])

    log(f"[5] top-{args.top_n_per_gene_per_ct} per (gene, ct): {len(final)} triples kept")
    if final:
        unique_pairs = {(t[0], t[2]) for t in final}
        unique_genes = {t[0] for t in final}
        unique_cts = {t[2] for t in final}
        log(f"     covering {len(unique_pairs)} (gene, ct) pairs, "
            f"{len(unique_genes)} genes, {len(unique_cts)} cell types")

    # 6. Write CSV
    with open(args.out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['gene', 'tf', 'cell_type', 'n_motif_sites', 'log2fc', 'pval_adj', 'score'])
        for g, tf, ct, n_sites, log2fc, pval_adj in final:
            score = (abs(log2fc) if not math.isnan(log2fc) else 0.0) * math.log10(n_sites + 1)
            w.writerow([g, tf, ct, n_sites, log2fc, pval_adj, f'{score:.6f}'])
    log(f"[6] wrote {args.out_csv}")

    # 7. Debug log
    Path(args.debug_log).write_text('\n'.join(log_lines) + '\n')


if __name__ == '__main__':
    main()
