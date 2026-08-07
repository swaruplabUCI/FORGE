#!/usr/bin/env python3
"""
Step 1b — pick barcodes and write the subset 10x multiome .h5.

Run in scgpu_extended.sif (h5py + scipy).

Design decisions, and why:

* The source .h5 is a COMBINED multiome matrix: 36,601 'Gene Expression' +
  111,743 'Peaks' features over 733,612 barcodes. We keep that structure so the
  subset is a drop-in for FORGE (bin/rna_qc.py calls sc.read_10x_h5, whose
  gex_only=True default filters to Gene Expression itself; CellBender consumes
  the raw .h5 directly).

* Features: keep ALL Gene Expression (annotation markers are genome-wide), and
  keep only Peaks on the target chromosomes so the peak axis is consistent with
  the chr21/22 fragment subset.

* Barcodes: keep N_CELLS real cells PLUS N_BACKGROUND low-count barcodes.
  Keeping only cells would break CellBender, which needs the ambient background
  distribution to estimate contamination. This is the single easiest thing to get
  wrong here.

* Cells are chosen by GEX UMI rank, intersected with barcodes that have at least
  MIN_FRAGS fragments on the target chromosomes, so both modalities cover the
  same cells (sample_id/barcode is the join key throughout FORGE).
"""
import argparse
import os
import h5py
import numpy as np
import scipy.sparse as sp

CHROMS = {"chr21", "chr22"}   # replaced from --chroms in main()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5", required=True)
    ap.add_argument("--frag-counts", required=True)
    ap.add_argument("--out-h5", required=True)
    ap.add_argument("--out-barcodes", required=True)
    ap.add_argument("--n-cells", type=int, default=1000)
    ap.add_argument("--n-background", type=int, default=19000)
    ap.add_argument("--min-frags", type=int, default=200)
    ap.add_argument("--chroms", default="chr21,chr22",
                    help="Comma-separated chromosomes; must match 01/03 "
                         "(default: chr21,chr22)")
    a = ap.parse_args()
    global CHROMS
    CHROMS = {c.strip() for c in a.chroms.split(",") if c.strip()}

    # ---- fragment counts per barcode ------------------------------------
    frag = {}
    with open(a.frag_counts) as fh:
        for line in fh:
            bc, c = line.rstrip("\n").split("\t")
            frag[bc] = int(c)
    print(f"barcodes with chr21/22 fragments: {len(frag):,}")

    with h5py.File(a.in_h5, "r") as f:
        g = f["matrix"]
        barcodes = np.array([b.decode() for b in g["barcodes"][:]])
        ftype = np.array([b.decode() for b in g["features"]["feature_type"][:]])
        interval = np.array([b.decode() for b in g["features"]["interval"][:]])
        n_feat, n_bc = g["shape"][:]
        print(f"source matrix: {n_feat:,} features x {n_bc:,} barcodes")

        # ---- feature mask ------------------------------------------------
        is_gex = ftype == "Gene Expression"
        peak_chrom = np.array([iv.split(":")[0] if ":" in iv else "" for iv in interval])
        keep_feat = is_gex | (~is_gex & np.isin(peak_chrom, list(CHROMS)))
        print(f"  keeping {is_gex.sum():,} GEX + "
              f"{(keep_feat & ~is_gex).sum():,} peaks on {sorted(CHROMS)} "
              f"= {keep_feat.sum():,} features")

        # ---- per-barcode GEX UMI (CSC: one column per barcode) -----------
        indptr = g["indptr"][:]
        indices = g["indices"]
        data = g["data"]
        gex_idx = np.flatnonzero(is_gex)
        gex_lo, gex_hi = gex_idx.min(), gex_idx.max() + 1
        assert np.array_equal(gex_idx, np.arange(gex_lo, gex_hi)), \
            "GEX features are not contiguous; per-barcode UMI shortcut invalid"

        umi = np.zeros(n_bc, dtype=np.int64)
        CH = 20000
        for lo in range(0, n_bc, CH):
            hi = min(lo + CH, n_bc)
            s, e = indptr[lo], indptr[hi]
            idx = indices[s:e]
            val = data[s:e]
            col = np.repeat(np.arange(lo, hi), np.diff(indptr[lo:hi + 1]))
            m = (idx >= gex_lo) & (idx < gex_hi)
            np.add.at(umi, col[m], val[m])
            if lo % 200000 == 0:
                print(f"    UMI scan {lo:,}/{n_bc:,}", flush=True)

        # ---- choose cells -------------------------------------------------
        fragvec = np.array([frag.get(bc, 0) for bc in barcodes])
        eligible = (fragvec >= a.min_frags) & (umi > 0)
        print(f"  eligible (>={a.min_frags} frags & UMI>0): {eligible.sum():,}")
        order = np.argsort(-umi)
        cells = [i for i in order if eligible[i]][:a.n_cells]
        if len(cells) < a.n_cells:
            print(f"  WARNING: only {len(cells)} cells met the criteria")

        # background: low-UMI barcodes, not already chosen. Take a spread from
        # the ambient regime so CellBender sees a realistic background profile.
        chosen = set(cells)
        bg_pool = [i for i in order if i not in chosen and 0 < umi[i] < 100]
        step = max(1, len(bg_pool) // a.n_background)
        background = bg_pool[::step][:a.n_background]
        print(f"  cells={len(cells):,} background={len(background):,} "
              f"(bg UMI range {umi[background].min()}-{umi[background].max()})")
        print(f"  cell UMI range {umi[cells].min()}-{umi[cells].max()}")

        keep_bc = np.array(sorted(cells + background))
        # ---- slice the matrix (CSC over barcodes) -------------------------
        cols = []
        for j in keep_bc:
            s, e = indptr[j], indptr[j + 1]
            cols.append((indices[s:e], data[s:e]))
        sub_indptr = np.zeros(len(keep_bc) + 1, dtype=np.int64)
        # remap feature ids
        fmap = -np.ones(n_feat, dtype=np.int64)
        fmap[np.flatnonzero(keep_feat)] = np.arange(keep_feat.sum())
        all_idx, all_val = [], []
        for k, (idx, val) in enumerate(cols):
            m = keep_feat[idx]
            ni, nv = fmap[idx[m]], val[m]
            o = np.argsort(ni)
            all_idx.append(ni[o]); all_val.append(nv[o])
            sub_indptr[k + 1] = sub_indptr[k] + len(ni)
        sub_indices = np.concatenate(all_idx) if all_idx else np.array([], dtype=np.int64)
        sub_data = np.concatenate(all_val) if all_val else np.array([], dtype=np.int32)
        print(f"  nonzeros: {len(sub_data):,} "
              f"(source {len(data):,}, {100*len(sub_data)/len(data):.2f}%)")

        # ---- write ---------------------------------------------------------
        os.makedirs(os.path.dirname(a.out_h5), exist_ok=True)
        with h5py.File(a.out_h5, "w") as o:
            m = o.create_group("matrix")
            m.create_dataset("barcodes",
                             data=np.array([barcodes[j].encode() for j in keep_bc]))
            m.create_dataset("data", data=sub_data.astype(np.int32), compression="gzip")
            m.create_dataset("indices", data=sub_indices.astype(np.int64),
                             compression="gzip")
            m.create_dataset("indptr", data=sub_indptr, compression="gzip")
            m.create_dataset("shape",
                             data=np.array([keep_feat.sum(), len(keep_bc)], dtype=np.int32))
            fo = m.create_group("features")
            src = g["features"]
            fkeep = np.flatnonzero(keep_feat)
            for k in src.keys():
                if k == "_all_tag_keys":
                    fo.create_dataset(k, data=src[k][:])
                else:
                    fo.create_dataset(k, data=src[k][:][fkeep], compression="gzip")

    with open(a.out_barcodes, "w") as fh:
        for j in keep_bc:
            kind = "cell" if j in chosen else "background"
            fh.write(f"{barcodes[j]}\t{kind}\t{umi[j]}\t{frag.get(barcodes[j], 0)}\n")

    mb = os.path.getsize(a.out_h5) / 1048576
    print(f"\nwrote {a.out_h5} ({mb:.1f} MB)")
    print(f"wrote {a.out_barcodes}")


if __name__ == "__main__":
    main()
