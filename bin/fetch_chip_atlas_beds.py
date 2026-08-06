#!/usr/bin/env python3
"""
fetch_chip_atlas_beds.py  —  build the per-TF ChIP reference for TF_CHIP_CORROBORATION

Produces a directory of <TF>.bed.gz files (merged ChIP-seq peaks per TF) that
TF_CHIP_CORROBORATION consumes as its empirical reference. This is the FORGE
analog of MAESTRO's CistromeDB dependency, but sourced from ChIP-Atlas, which
offers per-TF, per-assembly merged peak BEDs for mm10/mm39/hg38.

RUN THIS ONCE, OUTSIDE the pipeline (it needs network; pipeline processes do
not). Point params.tf_chip.chip_ref_dir at the output directory.

Two modes:
  --from-antibody-list : download ChIP-Atlas "Oanti" merged-peak BEDs per TF
                         (needs a TF list; queries the ChIP-Atlas peak browser).
  --from-local-bigbed  : convert an already-downloaded ChIP-Atlas allPeaks file,
                         splitting by antigen (TF) into per-TF BEDs (no network).

ChIP-Atlas assembly codes: mm10, mm39, hg38, hg19.
NOTE: FORGE's AD dataset is GRCm38==mm10; ChIP-Atlas mm10 aligns directly. For
mm39 datasets, either use ChIP-Atlas mm39 or liftOver (UCSC_Liftover/, per
project_scenicplus_liftover). Chromosome-name normalization ('chr' prefix) is
handled by the consumer.

This is a scaffold: the download URLs are documented inline so the exact call
can be adjusted to ChIP-Atlas's current API without touching the pipeline.
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict

# ChIP-Atlas peak-browser per-antigen merged BED (all cell types), by assembly:
#   https://chip-atlas.dbcls.jp/data/<assembly>/assembled/Oana.<assembly>.<TF>.AllAg.AllCell.bed
# and the genome-wide allPeaks table (large):
#   https://chip-atlas.dbcls.jp/data/<assembly>/allPeaks_light/allPeaks_light.<assembly>.<threshold>.bed.gz
CHIP_ATLAS_BASE = "https://chip-atlas.dbcls.jp/data"


def from_local_bigbed(allpeaks_path, tfs, outdir, assembly):
    """Split a ChIP-Atlas allPeaks_light BED into per-TF merged BEDs.

    The allPeaks_light format carries the antigen (TF) in column 4 as a
    semicolon-delimited metadata string, e.g. 'ID=...;Name=SPI1;...'. We keep
    only the wanted TFs (case-insensitive) and write chrom/start/end.
    """
    want = {t.lower() for t in tfs} if tfs else None
    buckets = defaultdict(list)
    opener = gzip.open if allpeaks_path.endswith(".gz") else open
    with opener(allpeaks_path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            meta = f[3]
            antigen = None
            for kv in meta.split(";"):
                if kv.lower().startswith("name="):
                    antigen = kv.split("=", 1)[1]
                    break
            if not antigen:
                continue
            if want is not None and antigen.lower() not in want:
                continue
            buckets[antigen].append((f[0], f[1], f[2]))

    os.makedirs(outdir, exist_ok=True)
    written = 0
    for tf, rows in buckets.items():
        out = os.path.join(outdir, f"{tf}.bed.gz")
        with gzip.open(out, "wt") as w:
            for c, s, e in rows:
                w.write(f"{c}\t{s}\t{e}\n")
        written += 1
        print(f"  {tf}: {len(rows)} peaks -> {out}", flush=True)
    print(f"Wrote {written} per-TF BEDs to {outdir} (assembly={assembly})", flush=True)
    if want:
        missing = want - {t.lower() for t in buckets}
        if missing:
            print(f"  NOT FOUND in allPeaks for {len(missing)} TFs: {sorted(missing)}", flush=True)


def from_antibody_list(tfs, outdir, assembly):
    """Download per-TF merged BEDs from ChIP-Atlas (needs `requests` + network).

    Scaffold only: prints the URLs it would fetch so the exact endpoint can be
    confirmed against ChIP-Atlas's current layout before wiring a real download.
    """
    try:
        import requests  # noqa: F401
        have_requests = True
    except ImportError:
        have_requests = False

    os.makedirs(outdir, exist_ok=True)
    for tf in tfs:
        url = f"{CHIP_ATLAS_BASE}/{assembly}/assembled/Oana.{assembly}.{tf}.AllAg.AllCell.bed"
        dest = os.path.join(outdir, f"{tf}.bed")
        if not have_requests:
            print(f"  [dry-run — no `requests`] would GET {url} -> {dest}", flush=True)
            continue
        import requests
        r = requests.get(url, timeout=120)
        if r.status_code == 200 and r.content:
            with open(dest, "wb") as w:
                w.write(r.content)
            print(f"  {tf}: {len(r.content)} bytes -> {dest}", flush=True)
        else:
            print(f"  {tf}: FAILED (HTTP {r.status_code}) {url}", flush=True)


def main():
    ap = argparse.ArgumentParser(description="Build per-TF ChIP reference for TF_CHIP_CORROBORATION")
    ap.add_argument("--assembly", required=True, choices=["mm10", "mm39", "hg38", "hg19"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tf-list", help="file with one TF symbol per line (optional filter)")
    ap.add_argument("--from-local-bigbed", help="path to a ChIP-Atlas allPeaks_light BED[.gz]")
    ap.add_argument("--from-antibody-list", action="store_true",
                    help="download per-TF merged BEDs from ChIP-Atlas")
    args = ap.parse_args()

    tfs = None
    if args.tf_list:
        with open(args.tf_list) as f:
            tfs = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]

    if args.from_local_bigbed:
        from_local_bigbed(args.from_local_bigbed, tfs, args.outdir, args.assembly)
    elif args.from_antibody_list:
        if not tfs:
            sys.exit("--from-antibody-list requires --tf-list")
        from_antibody_list(tfs, args.outdir, args.assembly)
    else:
        sys.exit("choose --from-local-bigbed <path> or --from-antibody-list")


if __name__ == "__main__":
    main()
