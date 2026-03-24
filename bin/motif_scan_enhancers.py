#!/usr/bin/env python3
"""
motif_scan_enhancers.py — Recipe Steps A5+A6

Scan enhancer regions for chromVAR-nominated TF motifs using scPRINTER's
MOODS-based motif scanner.  For each cell type's top TFs, extract enhancer
regions containing those motifs and write per-TF BED files suitable for
scPRINTER footprinting.

Inputs:
    - ccan_enhancer_peaks.bed.gz (from extract_ccan_enhancers.py)
    - per_celltype_motifs.json (from EXTRACT_CHROMVAR_MOTIFS)
    - JASPAR PFM file

Outputs:
    - enhancer_motif_scan.tsv.gz (region x TF presence matrix)
    - tf_enhancer_region_sets/ (per-TF BED files, 4-col: chr, start, end, linked_genes)
    - region_set_manifest.json
"""

import argparse
import json
from pathlib import Path

import scprinter as scp
import pandas as pd
import numpy as np


def init_scprinter(cache_dir):
    """Set scprinter dataset cache path."""
    dataset_obj = scp.datasets.datasets()
    dataset_obj.path = cache_dir
    return scp


def parse_jaspar_headers(pfms_path):
    """Parse JASPAR file headers to build motif_id -> TF name map.

    JASPAR format: >MA0489.2\tJun
    Returns dict: {'MA0489.2': 'Jun', ...}
    """
    headers = {}
    with open(pfms_path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split('\t')
                if len(parts) >= 2:
                    headers[parts[0]] = parts[1]
    return headers


def build_tf_name_map(all_tfs, motif_data, pfms_path):
    """Build mapping between chromVAR TF names and JASPAR motif names.

    Uses motif_id from chromVAR motif_details as bridge to JASPAR,
    with case-insensitive name fallback.

    Returns:
        jaspar_to_cv: dict mapping JASPAR motif name -> chromVAR TF name
        cv_to_jaspar: dict mapping chromVAR TF name -> JASPAR motif name
    """
    # Extract motif_id -> chromVAR TF from motif_data
    cv_motif_ids = {}  # motif_id -> cv_tf_name
    for ct_data in motif_data.get('cell_types', {}).values():
        for detail in ct_data.get('motif_details', []):
            cv_motif_ids[detail['motif_id']] = detail['tf']

    # Parse JASPAR headers
    jaspar_headers = parse_jaspar_headers(pfms_path)

    # Build maps
    jaspar_to_cv = {}
    cv_to_jaspar = {}
    matched = []
    unmatched = []

    for tf in all_tfs:
        # Strategy 1: motif_id bridge (most reliable)
        found = False
        for mid, cv_name in cv_motif_ids.items():
            if cv_name == tf and mid in jaspar_headers:
                jname = jaspar_headers[mid]
                jaspar_to_cv[jname] = tf
                cv_to_jaspar[tf] = jname
                matched.append(f"  {tf} -> {mid} -> JASPAR '{jname}'")
                found = True
                break

        if found:
            continue

        # Strategy 2: case-insensitive name match
        for jid, jname in jaspar_headers.items():
            if jname.lower() == tf.lower():
                jaspar_to_cv[jname] = tf
                cv_to_jaspar[tf] = jname
                matched.append(f"  {tf} -> name match -> JASPAR '{jname}' ({jid})")
                found = True
                break

        if not found:
            mid_for_tf = next((m for m, t in cv_motif_ids.items() if t == tf), 'N/A')
            unmatched.append(f"  {tf} (motif_id={mid_for_tf}) -- no JASPAR match")

    print(f"\nTF -> JASPAR mapping ({len(matched)} matched, {len(unmatched)} unmatched):")
    for m in matched:
        print(m)
    for u in unmatched:
        print(u)

    return jaspar_to_cv, cv_to_jaspar


def parse_args():
    p = argparse.ArgumentParser(description="Motif scan on enhancer regions")
    p.add_argument("--enhancer-peaks", required=True, help="Enhancer peaks BED.gz")
    p.add_argument("--motif-list", required=True, help="Per-celltype motifs JSON")
    p.add_argument("--pfms", required=True, help="JASPAR PFM file path")
    p.add_argument("--cache-dir", required=True, help="scPrinter cache dir")
    p.add_argument("--genome", default="hg38", help="Genome key")
    p.add_argument("--cpus", type=int, default=4, help="Number of parallel jobs for MOODS")
    p.add_argument("--min-regions", type=int, default=10,
                   help="Minimum regions per TF to write a region set")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    region_sets_dir = outdir / "tf_enhancer_region_sets"
    region_sets_dir.mkdir(parents=True, exist_ok=True)

    # Load enhancer peaks (5-col BED: chr, start, end, CCAN_id, linked_genes)
    print(f"Loading enhancer peaks from {args.enhancer_peaks}")
    enh_df = pd.read_csv(
        args.enhancer_peaks, sep='\t', header=None,
        names=['chr', 'start', 'end', 'CCAN_id', 'linked_genes'],
        compression='gzip'
    )
    print(f"  Loaded {len(enh_df)} enhancer peaks")

    if enh_df.empty:
        print("WARNING: No enhancer peaks to scan. Writing empty outputs.")
        pd.DataFrame().to_csv(outdir / "enhancer_motif_scan.tsv.gz", sep='\t',
                              index=False, compression='gzip')
        with open(outdir / "region_set_manifest.json", 'w') as f:
            json.dump({"region_sets": []}, f, indent=2)
        (outdir / "motif_scan_summary.txt").write_text("No enhancer peaks to scan\n")
        return

    # Load motif list
    print(f"Loading motif list from {args.motif_list}")
    with open(args.motif_list) as f:
        motif_data = json.load(f)

    all_tfs = motif_data.get('all_unique_tfs', [])
    ct_tfs = motif_data.get('cell_types', {})
    print(f"  {len(all_tfs)} unique TFs across {len(ct_tfs)} cell types")

    # Build TF name map (chromVAR <-> JASPAR)
    jaspar_to_cv, cv_to_jaspar = build_tf_name_map(all_tfs, motif_data, args.pfms)

    # Initialize scPRINTER for genome access
    global scp
    scp = init_scprinter(args.cache_dir)
    genome_obj = getattr(scp.genome, args.genome)

    # Build motif scanner from JASPAR PFMs
    jaspar_names_needed = list(cv_to_jaspar.values())
    motif_hits = {}  # cv_tf_name -> set of peak indices

    if not jaspar_names_needed:
        print("WARNING: No TFs matched JASPAR motifs. No region sets will be written.")
    else:
        print(f"\nBuilding motif scanner from {args.pfms}")
        motif_scanner = scp.motifs.Motifs(
            ref_path_motif=args.pfms,
            ref_path_fa=genome_obj.fetch_fa(),
            bg="even",
            n_jobs=args.cpus,
            motif_name_func=lambda x: x.split('\t')[-1],
        )

        # Filter scanner to only our TFs of interest
        motif_scanner.prep_scanner(tf_genes=jaspar_names_needed)
        print(f"  Scanner prepared for {len(jaspar_names_needed)} motifs: {jaspar_names_needed}")

        # Scan enhancer peaks for motif occurrences
        print(f"\nScanning {len(enh_df)} enhancer regions for motifs...")
        peaks_list = enh_df[['chr', 'start', 'end']].values.tolist()

        try:
            hits = motif_scanner.scan_motif(peaks_list, verbose=True, split_tfs=True)

            # Group hits by chromVAR TF name
            for hit in hits:
                peak_idx = hit[3]
                jaspar_tf = hit[4]
                cv_tf = jaspar_to_cv.get(jaspar_tf)
                if cv_tf:
                    motif_hits.setdefault(cv_tf, set()).add(peak_idx)

            print(f"\nMotif scan results:")
            for tf in all_tfs:
                n = len(motif_hits.get(tf, set()))
                if n >= args.min_regions:
                    status = f"{n} regions"
                elif n > 0:
                    status = f"{n} regions (below min_regions={args.min_regions})"
                elif tf in cv_to_jaspar:
                    status = "0 regions (motif not found in enhancers)"
                else:
                    status = "skipped (no JASPAR match)"
                print(f"  {tf}: {status}")

        except Exception as e:
            print(f"ERROR: Motif scanning failed: {e}")
            import traceback
            traceback.print_exc()
            print("No region sets will be written.")

    # Build motif scan matrix for output
    scan_records = []
    for tf in sorted(motif_hits.keys()):
        for idx in sorted(motif_hits[tf]):
            row = enh_df.iloc[idx]
            scan_records.append({
                'chr': row['chr'],
                'start': row['start'],
                'end': row['end'],
                'CCAN_id': row['CCAN_id'],
                'TF': tf,
            })

    scan_df = pd.DataFrame(scan_records)
    scan_path = outdir / "enhancer_motif_scan.tsv.gz"
    scan_df.to_csv(scan_path, sep='\t', index=False, compression='gzip')
    print(f"\nWrote motif scan results: {len(scan_df)} hits")

    # Write per-(cell_type, TF) BED files (4-col: chr, start, end, linked_genes)
    manifest_entries = []

    for ct_name, ct_data in ct_tfs.items():
        top_tfs = ct_data.get('top_tfs', [])
        for tf in top_tfs:
            if tf not in motif_hits or len(motif_hits[tf]) < args.min_regions:
                continue

            indices = sorted(motif_hits[tf])
            bed_records = []
            for idx in indices:
                row = enh_df.iloc[idx]
                genes = str(row.get('linked_genes', ''))
                if genes == 'nan':
                    genes = ''
                bed_records.append({
                    'Chromosome': row['chr'],
                    'Start': int(row['start']),
                    'End': int(row['end']),
                    'linked_genes': genes,
                })

            bed_df = pd.DataFrame(bed_records)

            # Sanitize names for filesystem (:: breaks Nextflow file() path handling)
            ct_safe = ct_name.replace(' ', '_').replace('/', '_').replace('+', 'pos')
            tf_safe = tf.replace(' ', '_').replace('/', '_').replace('::', '_x_')
            bed_fname = f"{ct_safe}__{tf_safe}.bed"
            bed_path = region_sets_dir / bed_fname
            bed_df.to_csv(bed_path, sep='\t', index=False, header=False)

            manifest_entries.append({
                'cell_type': ct_name,
                'tf': tf,
                'bed_file': bed_fname,
                'n_regions': len(bed_df),
                'source': 'jaspar_motif_scan',
            })

    # Write manifest
    manifest = {
        'region_sets': manifest_entries,
        'total_cell_types': len(ct_tfs),
        'total_tfs': len(all_tfs),
        'min_regions': args.min_regions,
        'matched_tfs': list(cv_to_jaspar.keys()),
        'unmatched_tfs': [tf for tf in all_tfs if tf not in cv_to_jaspar],
    }
    manifest_path = outdir / "region_set_manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWrote {len(manifest_entries)} region sets to manifest")

    # Summary
    summary_lines = [
        "Motif Scan Enhancers Summary",
        "=" * 40,
        f"Enhancer peaks scanned: {len(enh_df)}",
        f"TFs requested: {len(all_tfs)}",
        f"TFs matched to JASPAR: {len(cv_to_jaspar)}",
        f"TFs unmatched: {[tf for tf in all_tfs if tf not in cv_to_jaspar]}",
        f"TFs with motif hits (>= {args.min_regions}): "
        f"{sum(1 for tf in motif_hits if len(motif_hits[tf]) >= args.min_regions)}",
        f"Region sets written: {len(manifest_entries)}",
        "",
        "Per-cell-type region sets:",
    ]
    for ct_name in sorted(ct_tfs.keys()):
        ct_sets = [e for e in manifest_entries if e['cell_type'] == ct_name]
        tfs_str = ', '.join(f"{e['tf']}({e['n_regions']})" for e in ct_sets) if ct_sets else 'none'
        summary_lines.append(f"  {ct_name}: {len(ct_sets)} TFs — {tfs_str}")

    summary_path = outdir / "motif_scan_summary.txt"
    summary_path.write_text('\n'.join(summary_lines))
    print('\n'.join(summary_lines))


if __name__ == '__main__':
    main()
