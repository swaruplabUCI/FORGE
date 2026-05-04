#!/usr/bin/env python3
"""
plot_marker_coverage_tracks.py — Shi Fig 1E equivalent.

Render pyGenomeTracks coverage panels at canonical marker gene promoters,
one row per broad cell-type bigWig.

Inputs:
  --bigwig-manifest    JSON mapping lineage_label -> bigwig_path. Either:
                       (a) flat dict: {"Astrocyte": "/path/to/Astro-Epen.bw", ...}
                       (b) snap.ex manifest with key "bigwigs" => same map.
  --markers            JSON {"Astrocyte": "Gfap", ...} or comma list "Astrocyte=Gfap,..."
  --gtf                GTF for TSS lookup
  --window             ± bp around TSS to display (default 5000)

Outputs:
  <gene>_coverage.{pdf,png}    — one figure per marker gene
  marker_tracks_manifest.json
"""

import argparse
import gzip
import json
import re
import subprocess
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bigwig-manifest", required=True)
    p.add_argument("--markers", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--window", type=int, default=5000)
    p.add_argument("--outdir", default=".")
    p.add_argument("--width", type=int, default=18)
    p.add_argument("--dpi", type=int, default=200)
    return p.parse_args()


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "r")


def parse_markers(spec):
    p = Path(spec)
    if p.exists():
        return json.loads(p.read_text())
    out = {}
    for kv in spec.split(","):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_bigwigs(spec):
    obj = json.loads(Path(spec).read_text())
    if isinstance(obj, dict) and "bigwigs" in obj:
        return obj["bigwigs"]
    return obj


def gene_window(gtf_path, gene_name, window):
    pat = re.compile(rf'gene_name "{re.escape(gene_name)}"')
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            if not pat.search(parts[8]):
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            tss = start if strand == "+" else end - 1
            return chrom, max(0, tss - window), tss + window
    return None


def write_ini(bigwigs, ini_path):
    lines = []
    for label, path in bigwigs.items():
        # pyGenomeTracks section names cannot contain '/'; sanitize.
        name = re.sub(r"[^A-Za-z0-9_]", "_", label)
        lines.extend([
            f"[bw_{name}]",
            f"file = {path}",
            f"title = {label}",
            "height = 1.5",
            "color = #1f77b4",
            "min_value = 0",
            "show_data_range = true",
            "number_of_bins = 700",
            "type = fill",
            "",
        ])
    lines.extend([
        "[x-axis]",
        "fontsize = 8",
        "",
    ])
    Path(ini_path).write_text("\n".join(lines))


def render(ini_path, region, out_path, width, dpi):
    cmd = [
        "pyGenomeTracks",
        "--tracks", ini_path,
        "--region", region,
        "--outFileName", out_path,
        "--width", str(width),
        "--dpi", str(dpi),
        "--fontSize", "10",
    ]
    try:
        r = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"[1E][WARN] pyGenomeTracks failed for {region}: {e.stderr[:400]}", flush=True)
        return False


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bigwigs = parse_bigwigs(args.bigwig_manifest)
    markers = parse_markers(args.markers)
    if not bigwigs:
        raise SystemExit("[1E] empty bigwig manifest")
    if not markers:
        raise SystemExit("[1E] empty marker list")

    ini_path = outdir / "marker_tracks.ini"
    write_ini(bigwigs, ini_path)

    rendered = []
    for lineage, gene in markers.items():
        coords = gene_window(args.gtf, gene, args.window)
        if not coords:
            print(f"[1E] gene {gene} not found in GTF; skipping", flush=True)
            continue
        chrom, start, end = coords
        region = f"{chrom}:{start}-{end}"
        for ext in ("pdf", "png"):
            out = outdir / f"{gene}_coverage.{ext}"
            ok = render(str(ini_path), region, str(out), args.width, args.dpi)
            if ok:
                rendered.append({"lineage": lineage, "gene": gene,
                                "region": region, "file": str(out)})

    if not rendered:
        raise SystemExit("[1E] no marker tracks rendered (check bigwig paths and gene names)")

    (outdir / "marker_tracks_manifest.json").write_text(json.dumps({
        "rendered": rendered,
        "bigwig_manifest": bigwigs,
        "markers": markers,
        "window": args.window,
    }, indent=2))
    print(f"[1E] rendered {len(rendered)} marker track files", flush=True)


if __name__ == "__main__":
    main()
