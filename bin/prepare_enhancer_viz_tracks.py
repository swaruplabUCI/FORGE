#!/usr/bin/env python3
"""
prepare_enhancer_viz_tracks.py — Recipe D Steps D1-D2

Convert upstream pipeline outputs (Cicero connections, CCAN enhancers,
pseudobulk bigWigs, GTF) into pyGenomeTracks-compatible formats:
  - cicero_links.bedpe   (CCAN co-accessibility arcs)
  - enhancer_annotations.bed (highlighted enhancer regions)
  - tracks_<gene>.ini    (pyGenomeTracks INI config per locus)

Inputs:
  --cicero-connections   Cicero connections CSV (Peak1, Peak2, coaccess)
  --enhancer-peaks       CCAN enhancer peaks BED (from EXTRACT_CCAN_ENHANCERS)
  --bigwig-dir           Directory containing per-cell-type .bw files
  --gtf                  Gene annotation GTF (for gene models track)
  --target-genes         Comma-separated gene list (or JSON file)
  --motif-scan-manifest  region_set_manifest.json from MOTIF_SCAN_ENHANCERS
  --coaccess-threshold   Minimum co-accessibility score for links (default 0.25)
  --window-kb            Window around TSS in kb (default 100)
  --outdir               Output directory (default .)
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np


def peak_to_coords(peak_str):
    """Parse 'chr1-1000-2000' or 'chr1:1000-2000' to (chr, start, end)."""
    parts = peak_str.replace(':', '-').split('-')
    return parts[0], int(parts[1]), int(parts[2])


def export_cicero_links(cicero_csv, threshold, outdir):
    """Convert Cicero connections to BEDPE format for pyGenomeTracks arcs."""
    conns = pd.read_csv(cicero_csv, sep='\t', compression='infer')
    links_rows = []
    for _, row in conns.iterrows():
        if row['coaccess'] < threshold:
            continue
        chr1, s1, e1 = peak_to_coords(row['Peak1'])
        chr2, s2, e2 = peak_to_coords(row['Peak2'])
        if chr1 == chr2:
            links_rows.append([chr1, s1, e1, chr2, s2, e2, row['coaccess']])

    links_df = pd.DataFrame(links_rows,
        columns=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'score'])
    out_path = os.path.join(outdir, 'cicero_links.bedpe')
    links_df.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"[D1b] Wrote {len(links_df)} CCAN links to {out_path}")
    return out_path


def export_enhancer_bed(enhancer_peaks_path, outdir):
    """Convert enhancer peaks to BED6 for pyGenomeTracks annotation track."""
    # Read BED — may be gzipped
    enh = pd.read_csv(enhancer_peaks_path, sep='\t', header=None,
                       comment='#', compression='infer',
                       usecols=[0, 1, 2],
                       names=['chr', 'start', 'end'])
    enh = enh[['chr', 'start', 'end']].copy()
    enh['name'] = 'enhancer'
    enh['score'] = 1000
    enh['strand'] = '.'
    out_path = os.path.join(outdir, 'enhancer_annotations.bed')
    enh.to_csv(out_path, sep='\t', header=False, index=False)
    print(f"[D1d] Wrote {len(enh)} enhancer annotations to {out_path}")
    return out_path


def discover_bigwigs(bigwig_dir):
    """Find all .bw files in directory, map cell_type_name -> path."""
    bw_map = {}
    if not os.path.isdir(bigwig_dir):
        print(f"[WARN] bigwig_dir not found: {bigwig_dir}")
        return bw_map
    for fname in sorted(os.listdir(bigwig_dir)):
        if fname.endswith('.bw') or fname.endswith('.bigwig') or fname.endswith('.bigWig'):
            # Cell type name = filename without extension
            ct_name = os.path.splitext(fname)[0]
            bw_map[ct_name] = os.path.join(bigwig_dir, fname)
    print(f"[D1a] Found {len(bw_map)} bigWig files in {bigwig_dir}")
    return bw_map


def parse_tss_from_gtf(gtf_path):
    """Extract gene name -> (chr, tss) from GTF gene records."""
    tss_map = {}
    import gzip
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            # Extract gene_name from attributes
            attrs = fields[8]
            gene_name = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                    break
            if gene_name:
                tss = start if strand == '+' else end
                tss_map[gene_name] = (chrom, tss, strand)
    return tss_map


def write_pygenometracks_ini(
    config_path, gene_gtf, bigwig_files, links_file,
    enhancer_bed, links_label='CCANs', colors=None,
):
    """Generate a pyGenomeTracks INI config file."""
    lines = []

    # Gene model track
    lines.append('[genes]')
    lines.append(f'file = {gene_gtf}')
    lines.append('title = Genes')
    lines.append('height = 3')
    lines.append('fontsize = 10')
    lines.append('style = UCSC')
    lines.append('merge_transcripts = true')
    lines.append('')

    lines.append('[spacer]')
    lines.append('height = 0.5')
    lines.append('')

    # CCAN arc links
    if links_file and os.path.exists(links_file):
        lines.append('[links]')
        lines.append(f'file = {links_file}')
        lines.append(f'title = {links_label}')
        lines.append('height = 3')
        lines.append('color = darkblue')
        lines.append('line_width = 1')
        lines.append('links_type = arcs')
        lines.append('orientation = inverted')
        lines.append('')

    # Enhancer annotation
    if enhancer_bed and os.path.exists(enhancer_bed):
        lines.append('[enhancers]')
        lines.append(f'file = {enhancer_bed}')
        lines.append('title = Enhancers')
        lines.append('height = 0.5')
        lines.append('color = #ff6b35')
        lines.append('border_color = none')
        lines.append('display = collapsed')
        lines.append('')

    # BigWig coverage tracks per cell type
    if colors is None:
        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab10')
        colors = {}
        for i, k in enumerate(bigwig_files.keys()):
            c = cmap(i % 10)
            colors[k] = f'#{int(c[0]*255):02x}{int(c[1]*255):02x}{int(c[2]*255):02x}'

    for cell_type, bw_path in bigwig_files.items():
        section_name = cell_type.replace(' ', '_').replace('/', '_').replace('+', 'pos').replace('-', 'neg')
        lines.append(f'[{section_name}]')
        lines.append(f'file = {bw_path}')
        lines.append(f'title = {cell_type}')
        lines.append('height = 2')
        lines.append(f'color = {colors.get(cell_type, "#333333")}')
        lines.append('min_value = 0')
        lines.append('number_of_bins = 500')
        lines.append('')

    # X-axis
    lines.append('[x-axis]')
    lines.append('fontsize = 10')
    lines.append('')

    with open(config_path, 'w') as f:
        f.write('\n'.join(lines))
    return config_path


def resolve_target_genes(target_genes_arg, motif_scan_manifest):
    """Resolve target gene list from arg or motif scan manifest."""
    genes = []
    if target_genes_arg:
        # Could be a JSON file or comma-separated string
        if os.path.isfile(target_genes_arg):
            with open(target_genes_arg) as f:
                data = json.load(f)
                if isinstance(data, list):
                    genes = data
                elif isinstance(data, dict) and 'genes' in data:
                    genes = data['genes']
        else:
            genes = [g.strip() for g in target_genes_arg.split(',') if g.strip()]

    # Fallback: derive from motif scan manifest
    if not genes and motif_scan_manifest and os.path.isfile(motif_scan_manifest):
        with open(motif_scan_manifest) as f:
            manifest = json.load(f)
        # Collect unique TF names as gene targets
        seen = set()
        for entry in manifest.get('region_sets', []):
            tf = entry.get('tf', '')
            if tf and tf not in seen:
                genes.append(tf)
                seen.add(tf)
        print(f"[D1] Derived {len(genes)} target genes from motif scan manifest")

    return genes


def main():
    parser = argparse.ArgumentParser(description='Prepare enhancer visualization tracks (Recipe D Steps D1-D2)')
    parser.add_argument('--cicero-connections', required=True, help='Cicero connections CSV')
    parser.add_argument('--enhancer-peaks', required=True, help='CCAN enhancer peaks BED')
    parser.add_argument('--bigwig-dir', required=True, help='Directory with per-cell-type bigWig files')
    parser.add_argument('--gtf', required=True, help='Gene annotation GTF')
    parser.add_argument('--target-genes', default='', help='Comma-separated genes or JSON file')
    parser.add_argument('--motif-scan-manifest', default='', help='region_set_manifest.json')
    parser.add_argument('--coaccess-threshold', type=float, default=0.25)
    parser.add_argument('--window-kb', type=int, default=100)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Step D1b: Export Cicero links as BEDPE
    links_file = export_cicero_links(
        args.cicero_connections, args.coaccess_threshold, args.outdir)

    # Step D1d: Export enhancer annotations as BED
    enhancer_bed = export_enhancer_bed(args.enhancer_peaks, args.outdir)

    # Step D1a: Discover bigWig files
    bw_map = discover_bigwigs(args.bigwig_dir)

    # Parse TSS from GTF for per-gene region targeting
    tss_map = parse_tss_from_gtf(args.gtf)
    # Build case-insensitive lookup: JASPAR TF names may differ in case from GTF gene names
    tss_map_upper = {k.upper(): k for k in tss_map}
    print(f"[D1] Parsed {len(tss_map)} gene TSS positions from GTF")

    # Resolve target genes
    genes = resolve_target_genes(args.target_genes, args.motif_scan_manifest)
    if not genes:
        print("[WARN] No target genes specified and none derivable from manifest. "
              "Writing a single global tracks.ini.")
        genes = ['_global']

    # Step D2: Write per-gene INI configs
    gene_regions = {}
    for gene in genes:
        if gene == '_global':
            ini_path = os.path.join(args.outdir, 'tracks_global.ini')
        else:
            ini_path = os.path.join(args.outdir, f'tracks_{gene}.ini')

        # FIX-95: Use absolute paths in INI so it works from any Nextflow work dir
        bw_map_abs = {k: os.path.abspath(v) for k, v in bw_map.items()}
        write_pygenometracks_ini(
            config_path=ini_path,
            gene_gtf=os.path.abspath(args.gtf),
            bigwig_files=bw_map_abs,
            links_file=os.path.abspath(links_file) if links_file else links_file,
            enhancer_bed=os.path.abspath(enhancer_bed) if enhancer_bed else enhancer_bed,
            links_label='Cicero CCANs',
        )

        # Compute region string for this gene (case-insensitive GTF lookup)
        # Normalize gene name to GTF case (e.g. Atf3 -> ATF3)
        gene_gtf = tss_map_upper.get(gene.upper())
        if gene_gtf and gene != gene_gtf:
            print(f"[D2] Normalized '{gene}' -> '{gene_gtf}' (GTF case)")
            gene = gene_gtf
        if gene in tss_map:
            chrom, tss, strand = tss_map[gene]
            window = args.window_kb * 1000
            region_str = f"{chrom}:{max(0, tss - window)}-{tss + window}"
            gene_regions[gene] = {
                'ini': ini_path,
                'region': region_str,
                'chrom': chrom,
                'tss': tss,
                'strand': strand,
            }
        elif gene != '_global':
            print(f"[WARN] Gene '{gene}' not found in GTF — skipping region targeting")
            continue
        else:
            gene_regions[gene] = {'ini': ini_path, 'region': '', 'chrom': '', 'tss': 0, 'strand': '+'}

    # Write manifest JSON for downstream composite_enhancer_viz.py
    manifest = {
        'links_file': links_file,
        'enhancer_bed': enhancer_bed,
        'bigwig_dir': args.bigwig_dir,
        'bigwig_files': bw_map,
        'gene_regions': gene_regions,
        'gtf': args.gtf,
        'window_kb': args.window_kb,
        'coaccess_threshold': args.coaccess_threshold,
    }
    manifest_path = os.path.join(args.outdir, 'track_manifest.json')
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"[D2] Wrote track manifest for {len(gene_regions)} genes to {manifest_path}")

    print("[prepare_enhancer_viz_tracks] Done.")


if __name__ == '__main__':
    main()
