#!/usr/bin/env python3
"""
annotate_peak_types.py — classify peaks as Promoter / Exonic / Intronic / Distal.

Promoter = within +/- promoter-window of any TSS (Shi et al. convention, default 2 kb).
Exonic   = overlaps an exon (and not Promoter).
Intronic = inside a gene body but not Promoter or Exonic.
Distal   = everything else.

Inputs:
  --peaks    BED, BED.gz, TSV, or h5ad. For h5ad, peak ids are pulled from var_names.
  --gtf      GTF (or .gtf.gz) used to derive TSS, exons, gene bodies.

Outputs:
  peaks_annotated.tsv (peak_id, chr, start, end, peakType, nearestGene, distance_to_TSS)
  peaks_annotated.summary (one line per category with counts and percentages)
"""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--peaks", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--promoter-window", type=int, default=2000)
    p.add_argument("--out", default="peaks_annotated.tsv")
    p.add_argument("--summary", default="peaks_annotated.summary")
    return p.parse_args()


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def load_peaks(path):
    p = Path(path)
    if p.suffix == ".h5ad":
        import anndata
        adata = anndata.read_h5ad(p, backed="r")
        peaks = pd.DataFrame({"peak_id": list(adata.var_names)})
    else:
        with _open(p) as f:
            first = f.readline().rstrip("\n")
        has_header = not first.startswith(("chr", "Chr", "1", "2", "3"))
        peaks = pd.read_csv(p, sep="\t", header=0 if has_header else None,
                             dtype=str, comment="#")
        # Normalize col names. Accept BED-like (chr,start,end) or peak_id-bearing TSV.
        cols = list(peaks.columns)
        if cols[:3] == [0, 1, 2] or set(cols[:3]) == {0, 1, 2}:
            peaks = peaks.rename(columns={0: "chr", 1: "start", 2: "end"})
        elif "chr" in cols and "start" in cols and "end" in cols:
            pass
        elif "peak_id" in cols or "peak" in cols:
            pid_col = "peak_id" if "peak_id" in cols else "peak"
            peaks = peaks.rename(columns={pid_col: "peak_id"})
        else:
            peaks = peaks.iloc[:, :3]
            peaks.columns = ["chr", "start", "end"]

    if "peak_id" not in peaks.columns:
        peaks["peak_id"] = (
            peaks["chr"].astype(str) + ":" +
            peaks["start"].astype(str) + "-" +
            peaks["end"].astype(str)
        )
    if "chr" not in peaks.columns:
        m = peaks["peak_id"].astype(str).str.extract(r"^([^:]+):(\d+)-(\d+)$")
        m.columns = ["chr", "start", "end"]
        peaks = pd.concat([peaks[["peak_id"]], m], axis=1)

    peaks = peaks.dropna(subset=["chr", "start", "end"]).copy()
    peaks["start"] = peaks["start"].astype(int)
    peaks["end"] = peaks["end"].astype(int)
    return peaks[["peak_id", "chr", "start", "end"]].drop_duplicates("peak_id").reset_index(drop=True)


def parse_gtf(gtf_path):
    tss_rows, exon_rows, gene_rows = [], [], []
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("gene", "exon"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', parts[8]))
            gname = attrs.get("gene_name") or attrs.get("gene_symbol") or attrs.get("gene_id")
            if not gname:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            if feat == "gene":
                tss = start if strand == "+" else end - 1
                tss_rows.append((chrom, tss, tss + 1, gname, strand))
                gene_rows.append((chrom, start, end, gname, strand))
            else:
                exon_rows.append((chrom, start, end, gname, strand))

    def to_pr(rows):
        df = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "Name", "Strand"])
        return pr.PyRanges(df)

    return to_pr(tss_rows), to_pr(exon_rows), to_pr(gene_rows)


def main():
    args = parse_args()
    print(f"[annotate_peaks] reading peaks from {args.peaks}", flush=True)
    peaks = load_peaks(args.peaks)
    print(f"[annotate_peaks]   {len(peaks)} peaks", flush=True)

    print(f"[annotate_peaks] reading GTF {args.gtf}", flush=True)
    tss_pr, exon_pr, gene_pr = parse_gtf(args.gtf)
    print(f"[annotate_peaks]   tss={len(tss_pr.df)} exon={len(exon_pr.df)} gene={len(gene_pr.df)}",
          flush=True)

    # Harmonize chr-prefix between peaks and GTF
    peak_has_chr = peaks["chr"].astype(str).str.startswith("chr").any()
    gtf_has_chr = tss_pr.df["Chromosome"].astype(str).str.startswith("chr").any()
    if peak_has_chr and not gtf_has_chr:
        for prx in (tss_pr, exon_pr, gene_pr):
            prx.Chromosome = "chr" + prx.Chromosome.astype(str)
    elif (not peak_has_chr) and gtf_has_chr:
        peaks["chr"] = peaks["chr"].astype(str).str.replace("^chr", "", regex=True)

    pr_peaks = pr.PyRanges(peaks.rename(
        columns={"chr": "Chromosome", "start": "Start", "end": "End"}))

    # Promoter = TSS extended +/- window
    prom_df = tss_pr.df.copy()
    prom_df["Start"] = (prom_df["Start"] - args.promoter_window).clip(lower=0)
    prom_df["End"] = prom_df["End"] + args.promoter_window
    prom_pr = pr.PyRanges(prom_df)

    in_prom = set(pr_peaks.overlap(prom_pr).df["peak_id"].astype(str).unique())
    in_exon = set(pr_peaks.overlap(exon_pr).df["peak_id"].astype(str).unique())
    in_gene = set(pr_peaks.overlap(gene_pr).df["peak_id"].astype(str).unique())

    def classify(pid):
        if pid in in_prom:
            return "Promoter"
        if pid in in_exon:
            return "Exonic"
        if pid in in_gene:
            return "Intronic"
        return "Distal"

    peaks["peakType"] = peaks["peak_id"].map(classify)

    nearest = pr_peaks.nearest(tss_pr).df
    nearest = nearest.rename(columns={"Name": "nearestGene", "Distance": "distance_to_TSS"})
    nearest = nearest[["peak_id", "nearestGene", "distance_to_TSS"]].drop_duplicates("peak_id")
    out = peaks.merge(nearest, on="peak_id", how="left")
    out.to_csv(args.out, sep="\t", index=False)

    counts = out["peakType"].value_counts()
    n = len(out)
    lines = []
    for k in ("Promoter", "Exonic", "Intronic", "Distal"):
        cnt = int(counts.get(k, 0))
        pct = (cnt / n * 100.0) if n else 0.0
        lines.append(f"{k}\t{cnt}\t{pct:.2f}")
    Path(args.summary).write_text("peakType\tcount\tpercent\n" + "\n".join(lines) + "\n")
    print(f"[annotate_peaks] wrote {args.out}  +  {args.summary}", flush=True)
    print("[annotate_peaks] breakdown:")
    for ln in lines:
        print(f"  {ln}", flush=True)


if __name__ == "__main__":
    main()
