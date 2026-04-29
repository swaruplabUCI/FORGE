#!/usr/bin/env python3
"""
aggregate_fp_stats.py — per-(cell_type, TF) footprint scalar metrics.

Walks the per-(ct, TF) h5ads emitted by ENHANCER_FOOTPRINTING and rolls them
up into one CSV row per pair, then optionally joins by gene-window overlap to
produce a per-(gene, TF, ct) triple CSV.

The output columns are RAW measurements only — no composite/fused scores. Each
column means one specific thing. Users compose their own filter cascades from
these atoms; see `fp_stats_README.md` (also written by this script) for a
recommended workflow following the standard mechanistic-evidence hierarchy
(direct binding > co-accessibility linkage > differential activity).

OUTPUT COLUMNS (per-pair CSV; triple CSV adds gene + n_sites_in_gene_window):

Identifiers
  cell_type                          broad scATAnno label (33 ABCA classes)
  tf                                 TF / motif name (JASPAR2022_core)

Enhancer-TF binding evidence (mechanistic — strongest single signal that the
TF physically binds in this cell type. NOTE: 'enhancer' here means CCAN-
annotated peaks from MOTIF_SCAN_ENHANCERS. Promoter regions are NOT scanned
by this pipeline — see README §"Note on TF-promoter binding".)
  n_motif_sites_in_enhancers         total motif occurrences across all enhancer peaks
  n_sites_in_h5ad                    sites surviving scPRINTER's per-(ct) footprinting
  bind_dip_depth_mean                BINDetect-equivalent: mean(flank) − min(center)
                                     averaged across sites. Higher = stronger footprint.
  bind_dip_depth_p95                 95th percentile of per-site dip depth.
                                     Robust to outlier sites.
  n_sites_high_quality               count of sites with dip > --dip-threshold
                                     (default 0.05). Use as a "support count".

Cell-type specificity (does this TF preferentially bind in THIS ct vs others)
  fp_depth_mean                      scPRINTER footprint depth scalar, mean
  fp_depth_p95, fp_depth_max         distribution of per-site fp_depth
  fp_depth_rank_among_cts            rank of this ct among all 33 by fp_depth_mean
                                     for this TF. 1 = most-bound ct. Use when you
                                     want to find ct-specific binders.

Enhancer-gene linkage (regulatory plausibility — does this TF reach genes
through co-accessible enhancers?)
  n_genes_in_network                 number of target genes Cicero-linked to
                                     this TF in this ct via enhancer coaccessibility
  mean_coaccess_to_genes             mean Cicero coaccess score across linked genes

Differential accessibility (perturbation signal — only condition-aware metric)
  chromvar_log2fc                    log2 fold-change of chromVAR z-score
                                     (treatment vs control). Sign = direction.
  chromvar_p_adj                     BH-adjusted p-value from rank_genes_groups
  chromvar_signed_score              sign(log2fc) × −log10(p_adj). Use for sorting
                                     by directional significance.

Provenance
  evidence_tier                      'effect_curated' if chromvar_log2fc not NaN
                                     (i.e., the ct passed min_cells in both
                                     conditions and ran differential testing),
                                     else 'motif_only' (no diff data, but TF
                                     motif sites + footprint were measured)

Inputs:
  --fp-h5ads   FILE [FILE ...]   enhancer_footprints_*.h5ad files (collected)
  --tfbs-h5ads FILE [FILE ...]   enhancer_tfbs_*.h5ad files (collected, optional)
  --diff-csv-dir DIR             tf_differential_<ct>_<trt>_vs_<ctrl>.csv
  --treatment STR, --control STR for ct-name parsing of diff filenames
  --tf-gene-adjacency FILE       tf_gene_adjacency.tsv from BUILD_TF_GENE_NETWORK
  --motif-bed-dir DIR            tf_enhancer_region_sets/ for n_motif_sites count
  --region-set-manifest FILE     region_set_manifest.json
  --track-manifest FILE          track_manifest.json (for triple-form gene joining)
  --promoter-fp-dir DIR          results/scprinter/footprints/ root containing
                                 per-(ct, gene) {footprints,tfbs}_*.h5ad files.
                                 Optional; if provided, triple CSV gets Tier-1
                                 promoter-level evidence columns (see README).
  --dip-threshold FLOAT          per-site dip cutoff (default 0.05)
  --out-pair-csv FILE            (default fp_stats_per_pair.csv)
  --out-triple-csv FILE          (default fp_stats_per_triple.csv; written iff --track-manifest)
  --out-readme FILE              (default fp_stats_README.md; filter-workflow guide)
  --debug-log FILE               (default fp_stats.log)
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

import anndata as ad
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--fp-h5ads', nargs='+', required=True)
    p.add_argument('--tfbs-h5ads', nargs='*', default=[])
    p.add_argument('--diff-csv-dir', default='')
    p.add_argument('--treatment', default='')
    p.add_argument('--control', default='')
    p.add_argument('--tf-gene-adjacency', default='')
    p.add_argument('--motif-bed-dir', default='')
    p.add_argument('--region-set-manifest', default='')
    p.add_argument('--track-manifest', default='')
    p.add_argument('--promoter-fp-dir', default='')
    p.add_argument('--dip-threshold', type=float, default=0.05)
    p.add_argument('--out-pair-csv', default='fp_stats_per_pair.csv')
    p.add_argument('--out-triple-csv', default='fp_stats_per_triple.csv')
    p.add_argument('--out-readme', default='fp_stats_README.md')
    p.add_argument('--debug-log', default='fp_stats.log')
    return p.parse_args()


def _sanitize(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(s)) if s else ''


def _parse_ct_from_diff_filename(name, trt='', ctrl=''):
    """Same logic as build_viz_candidates.parse_ct_from_diff_filename. Recovers ct
    sanitized form from `tf_differential_<ct>_<sanitize(trt)>_vs_<sanitize(ctrl)>.csv`.
    """
    base = os.path.basename(name)
    if not base.startswith('tf_differential_') or not base.endswith('.csv'):
        return ''
    middle = base[len('tf_differential_'):-len('.csv')]
    if trt and ctrl:
        suffix = f'_{_sanitize(trt)}_vs_{_sanitize(ctrl)}'
        if middle.endswith(suffix):
            return middle[:-len(suffix)]
        return ''
    # Legacy fallback (single-token trt/ctrl)
    if '_vs_' not in middle:
        return ''
    pre_vs, _ = middle.rsplit('_vs_', 1)
    if '_' not in pre_vs:
        return ''
    ct, _ = pre_vs.rsplit('_', 1)
    return ct


def load_diff_csvs(diff_csv_dir, treatment, control):
    """Return {ct_sanitized: {tf: (log2fc, pval_adj)}}."""
    out = {}
    if not diff_csv_dir or not os.path.isdir(diff_csv_dir):
        return out
    for fname in os.listdir(diff_csv_dir):
        if not (fname.startswith('tf_differential_') and fname.endswith('.csv')):
            continue
        path = os.path.join(diff_csv_dir, fname)
        if os.path.getsize(path) < 200:
            continue
        ct_san = _parse_ct_from_diff_filename(fname, treatment, control)
        if not ct_san:
            continue
        per_tf = {}
        with open(path) as f:
            for row in csv.DictReader(f):
                m = row.get('motif', '')
                tf_name = m.split('\t', 1)[1].strip() if '\t' in m else m.strip()
                if not tf_name:
                    continue
                try:
                    log2fc = float(row.get('logfoldchange', 'nan'))
                    pval_adj = float(row.get('pval_adj', 'nan'))
                except ValueError:
                    continue
                if math.isnan(log2fc) or math.isnan(pval_adj):
                    continue
                prev = per_tf.get(tf_name)
                if prev is None or pval_adj < prev[1]:
                    per_tf[tf_name] = (log2fc, pval_adj)
        if per_tf:
            out[ct_san] = per_tf
    return out


def load_adjacency_counts(tsv_path):
    """Return {(ct, tf): (n_genes, mean_coaccess)} from tf_gene_adjacency.tsv."""
    if not tsv_path or not os.path.exists(tsv_path):
        return {}
    df = pd.read_csv(tsv_path, sep='\t')
    if df.empty:
        return {}
    grouped = df.groupby(['cell_type', 'TF']).agg(
        n_genes=('target_gene', 'count'),
        mean_coaccess=('coaccessibility', 'mean'),
    )
    return {(ct, tf): (int(row.n_genes), float(row.mean_coaccess))
            for (ct, tf), row in grouped.iterrows()}


def load_motif_bed_counts(motif_bed_dir, manifest_path):
    """{(ct, tf): n_motif_sites_in_enhancers} from per-(ct, tf) bed files."""
    if not manifest_path or not os.path.exists(manifest_path):
        return {}
    m = json.load(open(manifest_path))
    out = {}
    for r in m.get('region_sets', []):
        ct, tf = r['cell_type'], r['tf']
        bed = os.path.join(motif_bed_dir, r['bed_file']) if motif_bed_dir else None
        if bed and os.path.exists(bed):
            with open(bed) as f:
                out[(ct, tf)] = sum(1 for _ in f)
        else:
            out[(ct, tf)] = int(r.get('n_regions', 0))
    return out


def compute_dip_depth(tfbs_h5ad, primary_ct_idx):
    """Per-site dip depth = mean(flank) − min(center). Returns array (n_sites,)."""
    sites = list(tfbs_h5ad.obsm.keys())
    if not sites:
        return np.array([])
    out = []
    for s in sites:
        arr = tfbs_h5ad.obsm[s]
        if arr.ndim != 2 or arr.shape[0] <= primary_ct_idx:
            continue
        row = arr[primary_ct_idx, :]
        n = row.shape[0]
        flank_n = n // 4
        if flank_n <= 0 or n - 2 * flank_n <= 0:
            continue
        flank = np.concatenate([row[:flank_n], row[-flank_n:]])
        center = row[flank_n:-flank_n]
        out.append(float(flank.mean()) - float(center.min()))
    return np.array(out, dtype=np.float64)


def process_pair(fp_path, tfbs_path, dip_threshold):
    """Extract per-(primary_ct, TF) metrics from one (fp h5ad, tfbs h5ad) pair."""
    fp = ad.read_h5ad(fp_path)
    primary_ct = str(fp.uns.get('cell_type', ''))
    tf = str(fp.uns.get('tf', ''))
    if not primary_ct or not tf:
        return None
    obs_ix = list(fp.obs.index)
    if primary_ct not in obs_ix:
        return None
    ct_idx = obs_ix.index(primary_ct)
    n_sites = int(fp.shape[1])

    row = np.asarray(fp.X[ct_idx, :]).ravel()
    fp_mean = float(row.mean()) if len(row) else 0.0
    fp_p95 = float(np.quantile(row, 0.95)) if len(row) else 0.0
    fp_max = float(row.max()) if len(row) else 0.0

    all_means = np.asarray(fp.X.mean(axis=1)).ravel()
    rank = int(np.argsort(-all_means).tolist().index(ct_idx)) + 1 if len(all_means) else 0

    dip_mean = dip_p95 = 0.0
    n_high_q = 0
    if tfbs_path and os.path.exists(tfbs_path):
        try:
            tfbs = ad.read_h5ad(tfbs_path)
            tfbs_obs_ix = list(tfbs.obs.index)
            tfbs_idx = tfbs_obs_ix.index(primary_ct) if primary_ct in tfbs_obs_ix else ct_idx
            dip = compute_dip_depth(tfbs, tfbs_idx)
            if len(dip):
                dip_mean = float(dip.mean())
                dip_p95 = float(np.quantile(dip, 0.95))
                n_high_q = int((dip > dip_threshold).sum())
        except Exception as e:
            print(f"[WARN] tfbs load failed for {tfbs_path}: {e}", flush=True)

    return {
        'cell_type': primary_ct,
        'tf': tf,
        'n_sites_in_h5ad': n_sites,
        'fp_depth_mean': fp_mean,
        'fp_depth_p95': fp_p95,
        'fp_depth_max': fp_max,
        'fp_depth_rank_among_cts': rank,
        'bind_dip_depth_mean': dip_mean,
        'bind_dip_depth_p95': dip_p95,
        'n_sites_high_quality': n_high_q,
    }


def pair_tfbs_to_fp(fp_path, tfbs_paths_by_basename):
    """Match tfbs h5ad to fp h5ad by basename suffix swap."""
    base = os.path.basename(fp_path)
    swapped = base.replace('enhancer_footprints_', 'enhancer_tfbs_', 1)
    return tfbs_paths_by_basename.get(swapped)


def index_promoter_fp_dir(promoter_fp_dir):
    """Walk results/scprinter/footprints/{ct_safe}/{footprints,tfbs}_{ct}_{gene}.h5ad.

    Returns dict {(ct_label, gene): {'fp': fp_path, 'tfbs': tfbs_path}}. The ct
    LABEL inside the filename uses the unsanitized form ('OPC-Oligo', 'ABC NN'),
    matching the cell_type_broad column values used elsewhere. Filename pattern:
        footprints_{ct_label}_{gene}.h5ad
        tfbs_{ct_label}_{gene}.h5ad
    """
    out = {}
    if not promoter_fp_dir or not os.path.isdir(promoter_fp_dir):
        return out
    # Iterate per-ct subdir to disambiguate ct from gene (both can have spaces/_)
    for ct_safe in os.listdir(promoter_fp_dir):
        ct_dir = os.path.join(promoter_fp_dir, ct_safe)
        if not os.path.isdir(ct_dir):
            continue
        for fname in os.listdir(ct_dir):
            if not fname.endswith('.h5ad'):
                continue
            kind = None
            stem = None
            if fname.startswith('footprints_'):
                kind = 'fp'
                stem = fname[len('footprints_'):-len('.h5ad')]
            elif fname.startswith('tfbs_'):
                kind = 'tfbs'
                stem = fname[len('tfbs_'):-len('.h5ad')]
            else:
                continue
            # stem = "{ct_label}_{gene}" — split on LAST underscore so ct labels
            # with embedded underscores (e.g. 'IT-ET Glut') still parse cleanly.
            # Allen cell_type_broad labels use spaces in the unsanitized form,
            # but on disk the file basename keeps the space (no extra escaping).
            if '_' not in stem:
                continue
            ct_label, gene = stem.rsplit('_', 1)
            key = (ct_label, gene)
            entry = out.setdefault(key, {})
            entry[kind] = os.path.join(ct_dir, fname)
    return out


def compute_promoter_metrics(promoter_index, ct_label, gene):
    """Extract per-(ct, gene) promoter-level Tier-1 evidence from scPRINTER h5ads.

    Returns a dict with the four promoter columns (or NaN values if no h5ads
    exist for this key). Loads each h5ad lazily; called per-triple.

    Metrics (all per-(ct, gene), NOT per-TF):
      promoter_max_tfbs_score  Max value in tfbs obsm matrix (1 × ~380 positions)
                               across the 4kb promoter window. Range typically
                               0.002–0.9. >0.5 = strong TF binding signal
                               somewhere in the promoter for this ct.
      promoter_fp_depth_max    Max value in footprints obsm (1 × 99 scales ×
                               4000 bp). Strongest footprint anywhere in window.
      has_promoter_fp_data     True if both h5ads exist for this (ct, gene).
      promoter_fp_path         File path to the per-(ct, gene) footprint h5ad
                               for click-through inspection by the user.

    Why these are not TF-specific: the scPRINTER Bind model integrates over
    its motif catalog at each genomic position. To get per-TF dip metrics at
    promoters we'd need to run motif scanning over the 4kb window, which is
    not currently in the pipeline. Combined with TF-specific enhancer evidence
    (n_sites_high_quality, bind_dip_depth_mean), the user can build the
    composite logic: "TF has confirmed enhancer binding ∧ this gene's promoter
    is bound by something in this ct" — strong Tier-1+2 evidence.
    """
    out = {
        'promoter_max_tfbs_score': float('nan'),
        'promoter_fp_depth_max': float('nan'),
        'has_promoter_fp_data': False,
        'promoter_fp_path': '',
    }
    entry = promoter_index.get((ct_label, gene))
    if not entry:
        return out
    tfbs_path = entry.get('tfbs')
    fp_path = entry.get('fp')

    if tfbs_path and os.path.exists(tfbs_path):
        try:
            t = ad.read_h5ad(tfbs_path)
            if t.obsm:
                m = t.obsm[list(t.obsm.keys())[0]]
                out['promoter_max_tfbs_score'] = float(np.asarray(m).max())
        except Exception as e:
            print(f"[WARN] tfbs promoter load failed for {tfbs_path}: {e}", flush=True)

    if fp_path and os.path.exists(fp_path):
        try:
            f = ad.read_h5ad(fp_path)
            if f.obsm:
                m = f.obsm[list(f.obsm.keys())[0]]
                out['promoter_fp_depth_max'] = float(np.asarray(m).max())
            out['promoter_fp_path'] = fp_path
            out['has_promoter_fp_data'] = bool(tfbs_path)
        except Exception as e:
            print(f"[WARN] fp promoter load failed for {fp_path}: {e}", flush=True)

    return out


def main():
    args = parse_args()
    log = []

    def log_print(msg):
        log.append(msg)
        print(msg, flush=True)

    log_print(f"[1] fp_h5ads: {len(args.fp_h5ads)}, tfbs_h5ads: {len(args.tfbs_h5ads)}")
    log_print(f"[1] dip threshold: {args.dip_threshold}")

    # Index tfbs h5ads by basename for fast lookup
    tfbs_by_base = {os.path.basename(p): p for p in args.tfbs_h5ads}

    # Load joinable tables
    diff = load_diff_csvs(args.diff_csv_dir, args.treatment, args.control) if args.diff_csv_dir else {}
    log_print(f"[2] diff CSVs loaded for {len(diff)} cell types")

    adj_counts = load_adjacency_counts(args.tf_gene_adjacency) if args.tf_gene_adjacency else {}
    log_print(f"[3] tf_gene_adjacency entries: {len(adj_counts)} (ct, tf) pairs")

    motif_counts = load_motif_bed_counts(args.motif_bed_dir, args.region_set_manifest)
    log_print(f"[4] motif bed counts: {len(motif_counts)} (ct, tf) pairs")

    # Tier-1 promoter index — per-(ct, gene) scPRINTER outputs from
    # SCPRINTER_FOOTPRINTING. Used only by the triple CSV (which has a gene
    # column). pair CSV stays per-(ct, TF) and can't join on gene.
    promoter_index = index_promoter_fp_dir(args.promoter_fp_dir)
    log_print(f"[4b] promoter h5ad index: {len(promoter_index)} (ct, gene) pairs")

    # Process each pair
    rows = []
    for i, fp_path in enumerate(args.fp_h5ads):
        if not os.path.exists(fp_path):
            continue
        tfbs_path = pair_tfbs_to_fp(fp_path, tfbs_by_base)
        try:
            r = process_pair(fp_path, tfbs_path, args.dip_threshold)
        except Exception as e:
            log_print(f"[WARN] failed {fp_path}: {e}")
            continue
        if r is None:
            continue
        # Joins
        ct_san = _sanitize(r['cell_type'])
        diff_hit = (diff.get(ct_san) or diff.get(r['cell_type']) or {}).get(r['tf'])
        r['chromvar_log2fc'] = diff_hit[0] if diff_hit else float('nan')
        r['chromvar_p_adj'] = diff_hit[1] if diff_hit else float('nan')
        r['evidence_tier'] = 'effect_curated' if diff_hit else 'motif_only'

        adj_hit = adj_counts.get((r['cell_type'], r['tf']), (0, 0.0))
        r['n_genes_in_network'] = adj_hit[0]
        r['mean_coaccess_to_genes'] = adj_hit[1]

        r['n_motif_sites_in_enhancers'] = motif_counts.get((r['cell_type'], r['tf']), r['n_sites_in_h5ad'])

        # Single derived metric: chromvar_signed_score = sign(log2fc) × −log10(p_adj).
        # Useful as a directional sort key when filtering on differential signal.
        # Not a fused composite — purely a re-expression of the two chromvar columns.
        if not math.isnan(r['chromvar_log2fc']):
            r['chromvar_signed_score'] = (
                (1 if r['chromvar_log2fc'] >= 0 else -1)
                * -math.log10(max(r['chromvar_p_adj'], 1e-300))
            )
        else:
            r['chromvar_signed_score'] = float('nan')

        rows.append(r)
        if (i + 1) % 100 == 0:
            log_print(f"  processed {i+1}/{len(args.fp_h5ads)}")

    log_print(f"[5] {len(rows)} pair rows assembled")

    # Write per-pair CSV. Column ordering: identifiers → binding evidence →
    # ct-specificity → gene linkage → differential → provenance. Default sort
    # is by n_sites_high_quality desc (binding evidence is the priority axis
    # per the recommended filter workflow in fp_stats_README.md). Users
    # composing their own filters can sort by any column.
    if rows:
        df_pair = pd.DataFrame(rows)
        col_order = [
            # Identifiers
            'cell_type', 'tf',
            # Enhancer-TF binding evidence (mechanistic)
            'n_motif_sites_in_enhancers', 'n_sites_in_h5ad',
            'bind_dip_depth_mean', 'bind_dip_depth_p95', 'n_sites_high_quality',
            # Cell-type specificity
            'fp_depth_mean', 'fp_depth_p95', 'fp_depth_max', 'fp_depth_rank_among_cts',
            # Enhancer-gene linkage
            'n_genes_in_network', 'mean_coaccess_to_genes',
            # Differential accessibility (perturbation signal)
            'chromvar_log2fc', 'chromvar_p_adj', 'chromvar_signed_score',
            # Provenance
            'evidence_tier',
        ]
        df_pair = df_pair.reindex(columns=col_order)
        df_pair = df_pair.sort_values(
            ['n_sites_high_quality', 'bind_dip_depth_mean'],
            ascending=[False, False], kind='mergesort')
        df_pair.to_csv(args.out_pair_csv, index=False)
        log_print(f"[6] wrote {args.out_pair_csv} ({len(df_pair)} rows)")
    else:
        log_print(f"[6] no pair rows; skipping {args.out_pair_csv}")

    # Per-triple CSV — joins each pair to gene windows by motif-site overlap
    if args.track_manifest and os.path.exists(args.track_manifest) and rows:
        track = json.load(open(args.track_manifest))
        gene_windows = {}
        for g, info in (track.get('gene_regions') or {}).items():
            if g == '_global':
                continue
            region = info.get('region', '')
            if region and ':' in region:
                c, se = region.split(':')
                s, e = se.split('-')
                gene_windows[g] = (c, int(s), int(e))

        log_print(f"[7] gene_windows: {len(gene_windows)}")

        # For each (ct, tf), iterate the bed file and intersect with gene windows
        triple_rows = []
        manifest = json.load(open(args.region_set_manifest)) if args.region_set_manifest else {}
        bed_by_pair = {}
        for r in manifest.get('region_sets', []):
            bed_path = os.path.join(args.motif_bed_dir, r['bed_file']) if args.motif_bed_dir else None
            if bed_path and os.path.exists(bed_path):
                bed_by_pair[(r['cell_type'], r['tf'])] = bed_path

        # Index pairs for fast metric lookup
        pair_by_id = {(p['cell_type'], p['tf']): p for p in rows}

        # For each pair, count motif sites overlapping each gene window
        for (ct, tf), bed_path in bed_by_pair.items():
            pair_metrics = pair_by_id.get((ct, tf))
            if pair_metrics is None:
                continue
            with open(bed_path) as f:
                intervals = []
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            intervals.append((parts[0], int(parts[1]), int(parts[2])))
                        except ValueError:
                            continue
            for g, (gc, gs, ge) in gene_windows.items():
                n = sum(1 for c, s, e in intervals if c == gc and e > gs and s < ge)
                if n == 0:
                    continue
                triple = dict(pair_metrics)  # copy
                triple['gene'] = g
                triple['n_sites_in_gene_window'] = n
                # Tier-1 promoter join (per-(ct, gene), not per-TF). Adds
                # promoter_max_tfbs_score, promoter_fp_depth_max,
                # has_promoter_fp_data, promoter_fp_path. NaN if scPRINTER
                # didn't run for this (ct, gene).
                triple.update(compute_promoter_metrics(promoter_index, ct, g))
                triple_rows.append(triple)

        if triple_rows:
            df_triple = pd.DataFrame(triple_rows)
            # Identifiers + window-overlap + Tier-1 promoter (per-(ct, gene))
            # leading the row, then the per-pair columns in their established
            # order. Promoter columns sit between identifiers and the per-TF
            # enhancer columns to mirror the filter cascade in the README.
            promoter_cols = [
                'promoter_max_tfbs_score', 'promoter_fp_depth_max',
                'has_promoter_fp_data', 'promoter_fp_path',
            ]
            front = ['gene', 'cell_type', 'tf', 'n_sites_in_gene_window'] + promoter_cols
            tail = [c for c in df_pair.columns
                    if c not in ('cell_type', 'tf') and c not in promoter_cols]
            col_order = front + tail
            df_triple = df_triple.reindex(columns=col_order)
            # Sort key: by promoter binding strength first (Tier 1 — strongest
            # mechanistic signal of TF activity at this gene), then by enhancer
            # support. NaN promoter scores (no scPRINTER data) sink to bottom.
            df_triple = df_triple.sort_values(
                ['promoter_max_tfbs_score', 'n_sites_in_gene_window',
                 'n_sites_high_quality', 'bind_dip_depth_mean'],
                ascending=[False, False, False, False],
                na_position='last', kind='mergesort')
            df_triple.to_csv(args.out_triple_csv, index=False)
            log_print(f"[8] wrote {args.out_triple_csv} ({len(df_triple)} rows)")

    # Filter-workflow guide: explains each column and gives a tiered cascade
    # the user can copy-paste into pandas. No fused composite score is shipped
    # in the CSVs by design — see Option B / "no magic weights" rationale.
    Path(args.out_readme).write_text(_readme_text(args.dip_threshold))
    log_print(f"[9] wrote {args.out_readme}")

    Path(args.debug_log).write_text('\n'.join(log) + '\n')


def _readme_text(dip_threshold):
    """Filter workflow guide written alongside the CSVs.

    Tiered cascade matching the standard mechanistic-evidence hierarchy in
    regulatory genomics: physical binding > coaccessibility linkage >
    differential activity. Filter thresholds below are starting points;
    users should inspect the per-column distributions in their dataset
    before locking values down.
    """
    return f"""# fp_stats — filter workflow guide

This directory contains per-(cell_type, TF) and per-(gene, cell_type, TF) summary
statistics rolled up from ENHANCER_FOOTPRINTING. Use them to find interesting
regulators outside the top-3-per-(gene, ct) panels rendered by COMPOSITE_ENHANCER_VIZ.

## Files

- `fp_stats_per_pair.csv` — one row per (cell_type, TF). ~657 rows for this run.
- `fp_stats_per_triple.csv` — joined to gene windows: rows where TF motif sites
  overlap a target-gene window (gene ± window_kb). Thousands of rows.

## No fused composite score by design

These CSVs export only **raw measurements**. There is no `score_combined` column
because any single-number fusion requires arbitrary weights that don't port
across tissues / species / coverage levels. Compose your own filter cascade
from the columns below — see the recommended workflow at the bottom.

## Column reference

### Identifiers
| column | meaning |
|---|---|
| `cell_type` | broad scATAnno label (33 ABCA classes for mouse brain) |
| `tf` | TF / motif name from JASPAR2022_core |
| `gene` *(triple only)* | target gene window the TF motif sites overlap |
| `n_sites_in_gene_window` *(triple only)* | motif sites overlapping this gene's window |

### Tier 1 — Promoter-level binding evidence *(triple CSV only)*
**Source**: `SCPRINTER_FOOTPRINTING` outputs at
`results/scprinter/footprints/{{ct}}/{{footprints,tfbs}}_{{ct}}_{{gene}}.h5ad`.
Each h5ad covers one cell type × one ~4kb window centered on the gene's TSS
(strand-aware, from `gene_coordinates.json`).

These columns are **per-(ct, gene), NOT per-TF**: the scPRINTER Bind model
integrates over its motif catalog at each genomic position, so the per-position
score reflects "any TF" binding strength. To make them per-TF, combine with
the TF-specific enhancer columns below — e.g. "TF X has confirmed enhancer
binding (Tier-2 columns) ∧ this gene's promoter shows strong binding signal
in this ct (Tier-1 columns)" is the strongest practical Tier-1+2 filter.

| column | meaning |
|---|---|
| `promoter_max_tfbs_score` | Max value of the scPRINTER TFBS binding score across the 4kb promoter window (1 ct × ~380 positions). Range typically 0.002–0.9. **>0.5 = strong TF binding signal somewhere in this gene's promoter for this ct.** |
| `promoter_fp_depth_max` | Max value of the multi-scale footprint depth across the promoter window (99 scales × 4000 bp). Strongest footprint anywhere in the window. |
| `has_promoter_fp_data` | True if scPRINTER ran for this (ct, gene). False = the gene wasn't in `target_genes` or the ct didn't pass min_cells; treat as missing data, not as evidence of absence. |
| `promoter_fp_path` | Click-through path to the per-(ct, gene) footprint h5ad. Load with `anndata.read_h5ad()` to inspect the full per-position binding profile. |

### Tier 2 — Enhancer-TF binding evidence (per-TF, mechanistic)
**Source**: `ENHANCER_FOOTPRINTING` outputs at
`results/enhancer_footprinting/footprints/{{ct}}/{{tf}}/`. 'Enhancer' here means
CCAN-annotated peaks from MOTIF_SCAN_ENHANCERS — direct TF-binding evidence
in this ct's regulatory regions. These columns ARE per-TF (one row per TF).

| column | meaning |
|---|---|
| `n_motif_sites_in_enhancers` | total motif occurrences across all enhancer peaks (this ct) |
| `n_sites_in_h5ad` | sites surviving scPRINTER's per-(ct) footprinting |
| `bind_dip_depth_mean` | BINDetect-equivalent footprint metric: mean(flank counts) − min(center counts), averaged across sites. Larger = stronger footprint. Typical 0.01–0.5. |
| `bind_dip_depth_p95` | 95th percentile of per-site dip depth. Robust to outlier sites. |
| `n_sites_high_quality` | count of sites passing dip > {dip_threshold:.3f}. **Use as a "support count" — high count + high mean = confident footprint** |

### Cell-type specificity (find ct-specific binders)
| column | meaning |
|---|---|
| `fp_depth_mean` | scPRINTER footprint depth scalar, mean across sites |
| `fp_depth_p95` / `fp_depth_max` | distribution tail of fp_depth |
| `fp_depth_rank_among_cts` | rank of this ct among all 33 by `fp_depth_mean` for this TF. **`=1` means this is the most-bound ct for this TF** — strong signal of ct-specific binding |

### Tier 2 — Enhancer-gene linkage (regulatory plausibility — does this TF reach this gene?)
| column | meaning |
|---|---|
| `n_genes_in_network` | count of target genes Cicero-linked to this TF in this ct via co-accessible enhancers |
| `mean_coaccess_to_genes` | mean Cicero coaccess score across linked genes |

A TF with high `n_genes_in_network` is connected to many genes through enhancer
co-accessibility. Combined with the Tier-1 promoter columns, this lets you
distinguish "TF binds enhancers that loop to gene's promoter" from "TF binds
the promoter directly."

### Tier 3 — Differential accessibility (perturbation signal — only condition-aware metric)
**Note**: this is the only column set that uses condition information. The
binding and linkage metrics above are computed across all cells, stratified
by cell type — they tell you "what does this TF do here," not "what changed
in the perturbation." For perturbation-specific questions, weight differential
metrics higher in your filter than the standard mechanistic hierarchy implies.

| column | meaning |
|---|---|
| `chromvar_log2fc` | log2 fold-change of chromVAR z-score (treatment / control). Sign = direction. NaN if this ct didn't pass min_cells in differential testing. |
| `chromvar_p_adj` | Benjamini-Hochberg adjusted p-value from `rank_genes_groups` |
| `chromvar_signed_score` | sign(log2fc) × −log10(p_adj). Sort key for directional significance. |

### Provenance
| column | meaning |
|---|---|
| `evidence_tier` | `effect_curated` if differential testing ran for this ct (≥50 cells per group per condition); `motif_only` if it did not. Use to flag pairs without perturbation signal. |

## Recommended filter workflow

Follow the standard mechanistic-evidence hierarchy: **promoter binding
(strongest) > enhancer-TF binding > enhancer-gene linkage > perturbation
correlation (weakest, but the only condition-aware signal)**. Run as a cascade
on the **triple CSV** (the per-pair CSV doesn't have promoter columns since
those are per-(ct, gene), not per-(ct, TF)).

```python
import pandas as pd
df = pd.read_csv('fp_stats_per_triple.csv')

# --- Tier 1: TF-promoter binding evidence ---
# Strongest single signal. The 4kb promoter window of this gene shows strong
# TF binding signal in this cell type. Per-(ct, gene), not per-TF — combined
# with Tier 2 below it becomes "TF binds enhancers ∧ this gene's promoter is
# bound" which is the practical per-(ct, gene, TF) Tier-1+2 evidence.
df_t1 = df[
    df['has_promoter_fp_data'] &              # scPRINTER ran for this (ct, gene)
    (df['promoter_max_tfbs_score'] > 0.5)     # strong binding somewhere in window
]

# --- Tier 2: enhancer-TF binding (per-TF mechanism) ---
# This specific TF physically binds enhancers in this ct. Combined with
# Tier 1, you have promoter binding (any TF) + enhancer binding (this TF) —
# the strongest practical evidence the pipeline can give for "TF X regulates
# gene Y in cell type Z."
df_t2 = df_t1[
    (df_t1['n_sites_high_quality'] >= 5) &    # ≥5 confident enhancer footprints
    (df_t1['bind_dip_depth_mean'] > 0.02)     # depth: footprint visible
]

# --- Tier 2b: enhancer-gene linkage (regulatory reach) ---
# The TF reaches multiple genes via Cicero coaccess. Optional — high
# n_genes_in_network adds confidence the TF is broadly regulatory in this ct,
# but a TF can be a real regulator of a specific gene without linking widely.
df_t2b = df_t2[
    df_t2['n_genes_in_network'] >= 3
]

# --- Tier 3: differential accessibility (perturbation-aware) ---
# Optional. Use when prioritizing TFs that RESPOND to your perturbation.
# Drops 'motif_only' rows automatically since they have NaN log2fc.
df_t3 = df_t2b[
    (df_t2b['chromvar_p_adj'].fillna(1) < 0.05) &
    (df_t2b['chromvar_log2fc'].abs() > 0.5)
]

print(f"  starting:    {{len(df)}}")
print(f"  after T1:    {{len(df_t1)}}   (promoter binding)")
print(f"  after T2:    {{len(df_t2)}}   (+ enhancer binding for this TF)")
print(f"  after T2b:   {{len(df_t2b)}}  (+ enhancer-gene linkage ≥ 3 genes)")
print(f"  after T3:    {{len(df_t3)}}   (+ perturbation differential)")

# Inspect the survivors, sorted by promoter binding strength then ct-specificity
df_t3.sort_values(
    ['promoter_max_tfbs_score', 'fp_depth_rank_among_cts'],
    ascending=[False, True]
).head(20)
```

### Caveat: gene-list coverage gap

Tier 1 (promoter binding) requires that scPRINTER ran a per-promoter
footprinting task for the (ct, gene) you're inspecting. The `target_genes`
list passed to SCPRINTER_FOOTPRINTING is set independently of the visualization
target genes. In the current run, 446 genes were footprinted (from chromVAR-
discovery mode); for the ~17 visualization target genes you may have **zero
overlap** with that 446 list, so all `has_promoter_fp_data` will be False on
those rows. Two ways to fix:

1. *Quick*: set `params.scprinter.target_genes` explicitly to your viz target
   gene list (and re-launch with `-resume`). Only the new (ct × new_gene)
   tasks run; the existing 39k stay cached. Cost: ~17 genes × 33 cts = ~561
   new SCPRINTER_FOOTPRINTING tasks.
2. *Comprehensive*: union the chromVAR-discovery list with the viz list at
   runtime. Better long-term but a bigger change to the pipeline plumbing.

### When to deviate from this cascade

- **For perturbation discovery** (you want anything moving under your perturbation,
  even weak binders): start with Tier 3 first, then apply Tier 1 as a confidence
  filter.
- **For ct-specific marker TFs** (you want TFs that bind in *one* cell type): filter
  by `fp_depth_rank_among_cts == 1` AND high `n_sites_high_quality`, ignore tiers.
- **For very rare cell types** (n_cells < 50 in either condition): differential
  testing was skipped → all pairs in those cts are `motif_only`. Drop Tier 3 for
  those rows or you'll filter them all out.

### Adapting thresholds for other tissues / species

The default thresholds above (`n_sites_high_quality >= 5`, `bind_dip_depth_mean > 0.02`,
`n_genes_in_network >= 3`) are starting points calibrated for mouse brain scATAC-seq.
Before locking thresholds for a different dataset, plot the per-column distributions:

```python
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(12, 3))
df['n_sites_high_quality'].hist(ax=axes[0], bins=50); axes[0].set_xlabel('n_sites_high_quality')
df['bind_dip_depth_mean'].hist(ax=axes[1], bins=50);  axes[1].set_xlabel('bind_dip_depth_mean')
df['n_genes_in_network'].hist(ax=axes[2], bins=50);   axes[2].set_xlabel('n_genes_in_network')
```

Pick thresholds at natural distribution breakpoints (often around the 50th–75th
percentile for the "binding" metrics, lower for "linkage" because the network
is sparser).
"""


if __name__ == '__main__':
    main()
