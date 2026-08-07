#!/usr/bin/env python3
"""
Step 1a — extract chr21+chr22 fragments and tally fragments per barcode.

Run in snapatac_extended.sif (has pysam 0.23.3; the scgpu container does not).

Outputs into <out>/:
  frags_chr21_22.tsv     plain-text fragments for the target chromosomes
  frag_counts.tsv        barcode <TAB> n_fragments   (used to pick cells in 01b)

Barcodes are NOT filtered here: which barcodes are "cells" is decided in 01b from
RNA UMI counts intersected with these fragment counts.
"""
import argparse
import collections
import os
import pysam



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", default="chr21,chr22",
                    help="Comma-separated chromosomes to keep (default: chr21,chr22)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    CHROMS = [c.strip() for c in a.chroms.split(",") if c.strip()]

    tag = "_".join(CHROMS)
    tsv = os.path.join(a.out, f"frags_{tag}.tsv")
    counts = collections.Counter()
    n = 0

    tbx = pysam.TabixFile(a.fragments)
    present = [c for c in CHROMS if c in tbx.contigs]
    missing = [c for c in CHROMS if c not in tbx.contigs]
    if missing:
        raise SystemExit(f"Requested chromosomes absent from index: {missing}")
    print(f"Extracting {present} from {a.fragments}", flush=True)

    with open(tsv, "w") as out:
        for chrom in present:
            before = n
            for line in tbx.fetch(chrom):
                out.write(line + "\n")
                # fragments TSV: chrom, start, end, barcode, count
                counts[line.split("\t")[3]] += 1
                n += 1
            print(f"  {chrom}: {n - before:,} fragments", flush=True)

    with open(os.path.join(a.out, "frag_counts.tsv"), "w") as fh:
        for bc, c in counts.most_common():
            fh.write(f"{bc}\t{c}\n")

    print(f"\ntotal fragments : {n:,}")
    print(f"distinct barcodes: {len(counts):,}")
    print(f"plain tsv        : {tsv} "
          f"({os.path.getsize(tsv) / 1048576:.1f} MB uncompressed)")
    if counts:
        top = counts.most_common(5)
        print(f"top barcodes     : {top}")


if __name__ == "__main__":
    main()
