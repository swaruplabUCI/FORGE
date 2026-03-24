#!/usr/bin/env python3
"""
extract_signaling_targets.py — Recipe Step C4

For validated (receiver, TF) pairs, look up enhancer targets with 3-tier
fallback: eRegulon regions → DORC-motif intersection → CCAN-motif fallback.

Inputs:
    - validated_receiver_tf_pairs.tsv
    - eRegulon region sets (optional)
    - CCAN motif scan region sets (optional)

Outputs:
    - signaling_footprint_targets/ (per-pathway__receiver__TF BED files)
    - signaling_targets_manifest.json
    - signaling_metadata.tsv
"""

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Extract signaling target regions")
    p.add_argument("--validated-pairs", required=True, help="Validated receiver-TF pairs TSV")
    p.add_argument("--eregulon-regions", default=None, help="eRegulon region sets directory")
    p.add_argument("--eregulon-manifest", default=None, help="eRegulon manifest JSON")
    p.add_argument("--ccan-regions", default=None, help="CCAN motif scan region sets directory")
    p.add_argument("--ccan-manifest", default=None, help="CCAN region set manifest JSON")
    p.add_argument("--min-regions", type=int, default=10)
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def sanitize_name(name):
    """Sanitize a name for filesystem use."""
    return str(name).replace(' ', '_').replace('/', '_').replace('+', 'pos').replace(':', '_')


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    targets_dir = outdir / "signaling_footprint_targets"
    targets_dir.mkdir(parents=True, exist_ok=True)

    # Load validated pairs
    print(f"Loading validated pairs from {args.validated_pairs}")
    pairs_df = pd.read_csv(args.validated_pairs, sep='\t')
    print(f"  {len(pairs_df)} validated pairs")

    if pairs_df.empty:
        print("No validated pairs. Writing empty outputs.")
        with open(outdir / "signaling_targets_manifest.json", 'w') as f:
            json.dump({"targets": []}, f)
        pd.DataFrame(columns=['pathway', 'sender', 'receiver', 'tf', 'source', 'bed_file', 'n_regions']).to_csv(
            outdir / "signaling_metadata.tsv", sep='\t', index=False)
        return

    # Load available region sets
    ereg_lookup = {}
    if args.eregulon_manifest and Path(args.eregulon_manifest).exists():
        with open(args.eregulon_manifest) as f:
            ereg_manifest = json.load(f)
        ereg_dir = Path(args.eregulon_regions) if args.eregulon_regions else None
        for entry in ereg_manifest.get('region_sets', []):
            tf = entry['tf']
            bed_file = entry['bed_file']
            if ereg_dir and (ereg_dir / bed_file).exists():
                ereg_lookup[tf.upper()] = ereg_dir / bed_file
        print(f"  eRegulon region sets available: {len(ereg_lookup)}")

    ccan_lookup = {}  # (ct, tf) -> path
    if args.ccan_manifest and Path(args.ccan_manifest).exists():
        with open(args.ccan_manifest) as f:
            ccan_manifest = json.load(f)
        ccan_dir = Path(args.ccan_regions) if args.ccan_regions else None
        for entry in ccan_manifest.get('region_sets', []):
            ct = entry.get('cell_type', '')
            tf = entry['tf']
            bed_file = entry['bed_file']
            if ccan_dir and (ccan_dir / bed_file).exists():
                ccan_lookup[(ct.upper(), tf.upper())] = ccan_dir / bed_file
        print(f"  CCAN region sets available: {len(ccan_lookup)}")

    # Process each validated pair
    manifest_entries = []
    metadata_records = []

    for _, row in pairs_df.iterrows():
        pathway = str(row.get('pathway', 'unknown'))
        sender = str(row.get('sender', 'unknown'))
        receiver = str(row.get('receiver', 'unknown'))
        tf = str(row.get('downstream_TF', row.get('tf', '')))

        if not tf or tf == 'nan':
            continue

        # 3-tier fallback
        source = None
        src_path = None

        # Tier 1: eRegulon regions
        if tf.upper() in ereg_lookup:
            src_path = ereg_lookup[tf.upper()]
            source = 'eregulon'

        # Tier 2: CCAN-motif scan (cell-type specific)
        if src_path is None:
            key = (receiver.upper(), tf.upper())
            if key in ccan_lookup:
                src_path = ccan_lookup[key]
                source = 'ccan_motif'
            else:
                # Try any cell type
                for (ct, t), path in ccan_lookup.items():
                    if t == tf.upper():
                        src_path = path
                        source = 'ccan_motif_any_ct'
                        break

        if src_path is None:
            print(f"  No regions found for TF={tf} (receiver={receiver})")
            continue

        # Check minimum regions
        try:
            region_df = pd.read_csv(src_path, sep='\t', header=None)
            if len(region_df) < args.min_regions:
                print(f"  Skipping {tf}: only {len(region_df)} regions (< {args.min_regions})")
                continue
        except Exception as e:
            print(f"  WARNING: Failed to read {src_path}: {e}")
            continue

        # Copy to targets directory with descriptive name
        target_fname = f"{sanitize_name(pathway)}__{sanitize_name(receiver)}__{sanitize_name(tf)}.bed"
        target_path = targets_dir / target_fname
        shutil.copy2(src_path, target_path)

        manifest_entries.append({
            'pathway': pathway,
            'sender': sender,
            'receiver': receiver,
            'tf': tf,
            'source': source,
            'bed_file': target_fname,
            'n_regions': len(region_df),
        })

        metadata_records.append({
            'pathway': pathway,
            'sender': sender,
            'receiver': receiver,
            'tf': tf,
            'source': source,
            'bed_file': target_fname,
            'n_regions': len(region_df),
        })

    # Write manifest
    manifest = {'targets': manifest_entries}
    with open(outdir / "signaling_targets_manifest.json", 'w') as f:
        json.dump(manifest, f, indent=2)

    # Write metadata
    meta_df = pd.DataFrame(metadata_records) if metadata_records else pd.DataFrame(
        columns=['pathway', 'sender', 'receiver', 'tf', 'source', 'bed_file', 'n_regions'])
    meta_df.to_csv(outdir / "signaling_metadata.tsv", sep='\t', index=False)

    print(f"\nWrote {len(manifest_entries)} signaling target region sets")
    print(f"Sources: {pd.Series([e['source'] for e in manifest_entries]).value_counts().to_dict() if manifest_entries else 'none'}")


if __name__ == '__main__':
    main()
