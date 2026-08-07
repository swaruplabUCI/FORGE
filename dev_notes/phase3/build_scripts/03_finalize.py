#!/usr/bin/env python3
"""
Step 1c — filter fragments to the selected barcodes, bgzip + index, subset the
GTF and blacklist, and write the manifest.

Run in snapatac_extended.sif (pysam supplies bgzf writing and tabix indexing;
neither container ships the bgzip/tabix CLIs).
"""
import argparse
import gzip
import os
import pysam

CHROMS = {"chr21", "chr22"}   # replaced from --chroms in main()


def human(p):
    return f"{os.path.getsize(p) / 1048576:.1f} MB"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frags-tsv", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--gtf-in", required=True)
    ap.add_argument("--blacklist-in", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--sample-id", default="TUTORIAL_PBMC")
    ap.add_argument("--chroms", default="chr21,chr22",
                    help="Comma-separated chromosomes; must match 01/02 "
                         "(default: chr21,chr22)")
    a = ap.parse_args()
    global CHROMS
    CHROMS = {c.strip() for c in a.chroms.split(",") if c.strip()}

    samples = os.path.join(a.out_dir, "samples")
    refs = os.path.join(a.out_dir, "refs")
    os.makedirs(samples, exist_ok=True)
    os.makedirs(refs, exist_ok=True)

    keep = set()
    with open(a.barcodes) as fh:
        for line in fh:
            keep.add(line.split("\t")[0])
    print(f"selected barcodes: {len(keep):,}")

    # ---- fragments -> bgzf + tabix ---------------------------------------
    out_frag = os.path.join(samples, f"{a.sample_id}_atac_fragments.tsv.gz")
    kept = total = 0
    with pysam.BGZFile(out_frag, "wb") as out, open(a.frags_tsv) as src:
        for line in src:
            total += 1
            if line.split("\t")[3] in keep:
                out.write(line.encode())
                kept += 1
    print(f"fragments kept: {kept:,} / {total:,} ({100*kept/total:.1f}%)")
    print(f"  {out_frag} ({human(out_frag)})")

    pysam.tabix_index(out_frag, preset="bed", force=True, keep_original=True)
    tbi = out_frag + ".tbi"
    print(f"  {tbi} ({'ok' if os.path.exists(tbi) else 'MISSING'})")

    # ---- GTF subset -------------------------------------------------------
    out_gtf = os.path.join(refs, "gencode_chr21_22.gtf")
    n_in = n_out = 0
    op = gzip.open if a.gtf_in.endswith(".gz") else open
    with op(a.gtf_in, "rt") as src, open(out_gtf, "w") as out:
        for line in src:
            if line.startswith("#"):
                out.write(line)
                continue
            n_in += 1
            if line.split("\t", 1)[0] in CHROMS:
                out.write(line)
                n_out += 1
    print(f"GTF records: {n_out:,} / {n_in:,} kept -> {out_gtf} ({human(out_gtf)})")

    # ---- blacklist subset -------------------------------------------------
    out_bl = os.path.join(refs, "blacklist_chr21_22.bed")
    b_in = b_out = 0
    op = gzip.open if a.blacklist_in.endswith(".gz") else open
    with op(a.blacklist_in, "rt") as src, open(out_bl, "w") as out:
        for line in src:
            if line.startswith(("#", "track")):
                continue
            b_in += 1
            if line.split("\t", 1)[0] in CHROMS:
                out.write(line)
                b_out += 1
    print(f"blacklist: {b_out:,} / {b_in:,} kept -> {out_bl}")

    # ---- manifest ---------------------------------------------------------
    man = os.path.join(a.out_dir, "manifest.csv")
    with open(man, "w") as fh:
        fh.write("sample_id,batch,sample_type,original_lane_id,rna_file,"
                 "fragment_file,condition_group,data_dir\n")
        fh.write(f"{a.sample_id},tutorial,lane,L1,"
                 f"{a.sample_id}_raw_feature_bc_matrix.h5,"
                 f"{a.sample_id}_atac_fragments.tsv.gz,"
                 f"ConditionA,{samples}\n")
    print(f"manifest -> {man}")


if __name__ == "__main__":
    main()
