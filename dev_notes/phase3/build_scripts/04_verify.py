#!/usr/bin/env python3
"""
Step 1d — verify the built tutorial dataset is a valid FORGE input.

Checks the things that would otherwise fail hours into a run:
  1. scanpy can read the subset .h5 the way bin/rna_qc.py does (gex_only=True)
  2. the matrix has real cells AND an ambient background (CellBender needs both)
  3. the fragments file is bgzf, tabix-queryable, and well-formed
  4. RNA and ATAC barcodes actually overlap (they are the cross-modality join key)
  5. the GTF subset is non-empty and confined to the target chromosomes

Run in snapatac_extended.sif (needs pysam AND scanpy).
"""
import gzip
import os
import sys
import numpy as np
import pysam
import scanpy as sc

W = "/dfs7/swaruplab/lesolano/FORGE/oneOff/20260806_tutorial/out"
H5 = f"{W}/samples/TUTORIAL_PBMC_raw_feature_bc_matrix.h5"
FRAG = f"{W}/samples/TUTORIAL_PBMC_atac_fragments.tsv.gz"
GTF = f"{W}/refs/gencode_chr21_22.gtf"
CHROMS = {"chr21", "chr22"}

fails = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" — {detail}" if detail else ""))
    if not ok:
        fails.append(label)


print("=" * 64)
print("1. scanpy read_10x_h5 (as bin/rna_qc.py does)")
print("=" * 64)
ad = sc.read_10x_h5(H5)          # gex_only=True by default
ad.var_names_make_unique()
print(f"  shape: {ad.shape[0]:,} barcodes x {ad.shape[1]:,} genes")
check("GEX-only read works", ad.n_vars > 30000, f"{ad.n_vars:,} genes")
check("barcode count as built", ad.n_obs == 20000, f"{ad.n_obs:,}")

print("\n" + "=" * 64)
print("2. cells + ambient background present (CellBender requirement)")
print("=" * 64)
umi = np.asarray(ad.X.sum(axis=1)).ravel()
n_cells = int((umi >= 1000).sum())
n_bg = int(((umi > 0) & (umi < 100)).sum())
print(f"  UMI range: {umi.min():.0f} – {umi.max():.0f}")
check("has real cells (UMI>=1000)", n_cells >= 900, f"{n_cells:,}")
check("has ambient background (0<UMI<100)", n_bg >= 10000, f"{n_bg:,}")
check("background outnumbers cells", n_bg > n_cells, f"{n_bg:,} vs {n_cells:,}")

print("\n" + "=" * 64)
print("3. fragments: bgzf + tabix + well-formed")
print("=" * 64)
check("fragments exist", os.path.exists(FRAG), f"{os.path.getsize(FRAG)/1048576:.1f} MB")
check("tabix index exists", os.path.exists(FRAG + ".tbi"))
tbx = pysam.TabixFile(FRAG)
contigs = set(tbx.contigs)
check("only target chromosomes indexed", contigs == CHROMS, str(sorted(contigs)))

# Derive the probe window from the DATA rather than hardcoding one. A fixed
# window like chr21:1-2Mb lands in the acrocentric short arm, which is N-rich and
# unmappable, so it legitimately holds zero fragments — that produced a spurious
# FAIL, and worse, it made the two checks below pass vacuously on an empty list.
probe = sorted(contigs)[0]
all_starts = [int(r.split("\t")[1]) for r in tbx.fetch(probe)]
lo, hi = min(all_starts), max(all_starts)
mid = (lo + hi) // 2
rows = list(tbx.fetch(probe, mid, mid + 1_000_000))
print(f"  {probe} spans {lo:,}–{hi:,}; probing {mid:,}–{mid + 1_000_000:,}")
check("tabix range query returns rows", len(rows) > 0,
      f"{len(rows):,} rows in a populated 1 Mb window")

# Guard the next two so they cannot pass on an empty result.
if not rows:
    check("5 columns per fragment row", False, "no rows to check")
    check("coordinates sorted", False, "no rows to check")
else:
    bad = [r for r in rows if len(r.split("\t")) != 5]
    check("5 columns per fragment row", not bad,
          f"{len(bad)} malformed of {len(rows):,}")
    starts = [int(r.split("\t")[1]) for r in rows]
    check("coordinates sorted", starts == sorted(starts), f"{len(starts):,} rows")
    check("whole-contig ordering sorted", all_starts == sorted(all_starts),
          f"{len(all_starts):,} fragments on {probe}")

print("\n" + "=" * 64)
print("4. RNA / ATAC barcode overlap (cross-modality join key)")
print("=" * 64)
frag_bcs = set()
with gzip.open(FRAG, "rt") as fh:
    for i, line in enumerate(fh):
        frag_bcs.add(line.split("\t")[3])
        if i > 400000:
            break
rna_bcs = set(ad.obs_names)
inter = rna_bcs & frag_bcs
print(f"  RNA barcodes {len(rna_bcs):,} | ATAC (sampled) {len(frag_bcs):,}")
check("ATAC barcodes are a subset of RNA", frag_bcs <= rna_bcs,
      f"{len(frag_bcs - rna_bcs)} ATAC-only")
check("substantial overlap", len(inter) > 900, f"{len(inter):,} shared")

print("\n" + "=" * 64)
print("5. GTF subset")
print("=" * 64)
chroms, n = set(), 0
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        chroms.add(line.split("\t", 1)[0])
        n += 1
check("GTF non-empty", n > 10000, f"{n:,} records")
check("GTF confined to target chromosomes", chroms <= CHROMS, str(sorted(chroms)))

total = sum(os.path.getsize(os.path.join(dp, f))
            for dp, _, fs in os.walk(W) for f in fs)
print("\n" + "=" * 64)
print(f"TOTAL DATASET SIZE: {total/1048576:.1f} MB")
print("=" * 64)
print("ALL CHECKS PASSED" if not fails else f"FAILURES: {fails}")
sys.exit(1 if fails else 0)
