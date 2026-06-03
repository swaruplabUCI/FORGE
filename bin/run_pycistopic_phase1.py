#!/usr/bin/env python
"""
run_pycistopic_phase1.py — Phase 1 of the pyCisTopic 3-phase pipeline.

Runs:  pseudobulk export, MACS2 peak calling, consensus peak generation,
       per-sample pycistopic QC (parquet stats), and TSS BED generation.

Does NOT call create_cistopic_object_from_fragments() — that's Phase 2.

Outputs (written to --outdir):
  consensus_peak_calling/consensus_regions.bed
  consensus_peak_calling/pseudobulk_bw_files/*.bw
  qc/tss.bed
  qc/{sample_id}.fragments_stats_per_cb.parquet  (one per sample)
  cell_metadata_for_pycistopic.safe.tsv
  fragments_map.tsv
  group_list.tsv   (cell_type_safe, condition, n_cells — fan-out manifest)
  blacklist.bed    (decompressed if input was .gz)
"""

import argparse
import gc
import gzip
import json
import logging
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycistopic_phase1")


def parse_args():
    p = argparse.ArgumentParser(
        description="pyCisTopic Phase 1: pseudobulk, peak calling, QC."
    )
    p.add_argument("--fragments-map", required=True,
                   help="TSV with columns sample_id and fragments_path.")
    p.add_argument("--cell-metadata", required=True,
                   help="Cell metadata TSV with columns: barcode (index), "
                        "sample_id, cell_type_safe, condition.")
    p.add_argument("--species", required=True,
                   help="'mmusculus' or 'hsapiens'.")
    p.add_argument("--sample-id-col", default="sample_id")
    p.add_argument("--cell-type-col", default="cell_type_safe")
    p.add_argument("--condition-col", default="condition")
    p.add_argument("--variable", default="cell_type_safe",
                   help="Column to group pseudobulk by (CT only; not CT×cond).")
    p.add_argument("--n-cpu", type=int, default=8)
    p.add_argument("--blacklist-bed", default=None)
    p.add_argument("--gtf", default=None)
    p.add_argument("--min-cells", type=int, default=200,
                   help="Min cells per CT×condition group for group_list.tsv.")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def generate_tss_bed_from_gtf(gtf_path, output_bed, logger):
    logger.info("Generating TSS BED from GTF: %s", gtf_path)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    records = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, start, end, strand, attrs = (
                fields[0], int(fields[3]) - 1, int(fields[4]), fields[6], fields[8]
            )
            gene_name, gene_type = None, "unknown"
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1]
            if gene_name is None:
                continue
            tss = start if strand == "+" else end - 1
            records.append((chrom, tss, tss + 1, gene_name, 0, strand, gene_type))

    tss_df = pd.DataFrame(records, columns=[
        "# Chromosome", "Start", "End", "Gene", "Score", "Strand", "Transcript_type"
    ])
    if not tss_df.empty and not tss_df["# Chromosome"].iloc[0].startswith("chr"):
        tss_df["# Chromosome"] = "chr" + tss_df["# Chromosome"]
    tss_df.drop_duplicates(inplace=True)
    tss_df.to_csv(output_bed, sep="\t", index=False)
    logger.info("Wrote %d TSS entries to %s", len(tss_df), output_bed)


def main():
    args = parse_args()
    logger = setup_logger()
    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # ── Load cell metadata ──────────────────────────────────────────────────
    logger.info("Loading cell metadata from %s", args.cell_metadata)
    cell_data = pd.read_table(args.cell_metadata, index_col=0, sep="\t")

    if "barcode" in cell_data.columns:
        cell_data.index = cell_data["barcode"].astype(str)

    raw_idx = cell_data.index.astype(str)
    if raw_idx.str.contains(":").any():
        stripped = raw_idx.str.split(":", n=1).str[-1]
        n_changed = (stripped != raw_idx).sum()
        logger.info("Stripped sample-prefix from %d / %d barcodes.", n_changed, len(raw_idx))
        cell_data.index = stripped

    # Ensure required columns exist
    for col in [args.sample_id_col, args.cell_type_col]:
        if col not in cell_data.columns:
            raise ValueError(f"Required column '{col}' not found in cell metadata. "
                             f"Found: {list(cell_data.columns)}")

    if args.condition_col not in cell_data.columns:
        logger.warning("Condition column '%s' not found; per-group fan-out will be CT-only.",
                       args.condition_col)
        cell_data["_condition"] = "all"
        condition_col = "_condition"
    else:
        condition_col = args.condition_col

    # ── Load fragments map ──────────────────────────────────────────────────
    logger.info("Loading fragments map from %s", args.fragments_map)
    frag_map_df = pd.read_table(args.fragments_map, sep="\t")
    path_to_fragments = dict(zip(
        frag_map_df["sample_id"].astype(str),
        frag_map_df["fragments_path"].astype(str),
    ))

    sample_ids_in_cells = set(cell_data[args.sample_id_col].astype(str))
    path_to_fragments = {sid: p for sid, p in path_to_fragments.items()
                         if sid in sample_ids_in_cells and os.path.exists(p)}
    cell_data = cell_data[cell_data[args.sample_id_col].astype(str).isin(path_to_fragments)].copy()
    logger.info("Using %d samples / %d cells after fragment-file filtering.",
                len(path_to_fragments), len(cell_data))

    # ── Chromsizes ──────────────────────────────────────────────────────────
    if args.species in ("mmusculus", "mouse"):
        chrom_url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
        chrom_names = ["Chromosome", "End"]
    else:
        chrom_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
        chrom_names = ["Chromosome", "End"]

    chrom_file = os.path.join(outdir, "chromsizes.tsv")
    if not os.path.exists(chrom_file):
        logger.info("Downloading chromsizes for %s...", args.species)
        chromsizes = pd.read_table(chrom_url, header=None, names=chrom_names)
        chromsizes.insert(1, "Start", 0)
        chromsizes.to_csv(chrom_file, sep="\t", index=False)
    else:
        chromsizes = pd.read_table(chrom_file)

    # ── Blacklist ───────────────────────────────────────────────────────────
    if args.blacklist_bed:
        if args.blacklist_bed.endswith(".gz"):
            bl_out = os.path.join(outdir, "blacklist.bed")
            subprocess.check_call(f"gunzip -c '{args.blacklist_bed}' > '{bl_out}'", shell=True)
            path_to_blacklist = bl_out
        else:
            path_to_blacklist = args.blacklist_bed
    else:
        path_to_blacklist = os.path.join(
            os.path.dirname(pycisTopic.__file__), "blacklist", "hg38-blacklist.v2.bed"
        )
    logger.info("Blacklist: %s", path_to_blacklist)

    # Write decompressed blacklist into outdir so downstream phases can find it
    bl_out = os.path.join(outdir, "blacklist.bed")
    if not os.path.exists(bl_out):
        import shutil
        shutil.copy(path_to_blacklist, bl_out)

    # ── TSS BED ─────────────────────────────────────────────────────────────
    qc_outdir = os.path.join(outdir, "qc")
    ensure_dir(qc_outdir)
    tss_bed = os.path.join(qc_outdir, "tss.bed")
    if not os.path.exists(tss_bed):
        if args.gtf:
            generate_tss_bed_from_gtf(args.gtf, tss_bed, logger)
        else:
            logger.info("Generating TSS via pycistopic tss get_tss for %s", args.species)
            subprocess.check_call([
                "pycistopic", "tss", "get_tss",
                "--output", tss_bed, "--name", "hsapiens_gene_ensembl",
                "--to-chrom-source", "ucsc", "--ucsc", "hg38",
            ])

    # ── Pseudobulk export ───────────────────────────────────────────────────
    pb_outdir = os.path.join(outdir, "consensus_peak_calling")
    ensure_dir(pb_outdir)
    ensure_dir(os.path.join(pb_outdir, "pseudobulk_bed_files"))
    ensure_dir(os.path.join(pb_outdir, "pseudobulk_bw_files"))
    dropped_log = os.path.join(pb_outdir, "dropped_groups.tsv")

    logger.info("Exporting pseudobulk profiles (grouped by '%s')...", args.variable)

    def _sanitize(s):
        s = re.sub(r"\s+", "_", str(s))
        s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
        return re.sub(r"_+", "_", s).strip("_")

    while True:
        try:
            bw_paths, bed_paths = export_pseudobulk(
                input_data=cell_data,
                variable=args.variable,
                sample_id_col=args.sample_id_col,
                chromsizes=chromsizes,
                bed_path=os.path.join(pb_outdir, "pseudobulk_bed_files"),
                bigwig_path=os.path.join(pb_outdir, "pseudobulk_bw_files"),
                path_to_fragments=path_to_fragments,
                n_cpu=args.n_cpu,
                normalize_bigwig=True,
                temp_dir=tempfile.gettempdir(),
                split_pattern="-",
            )
            break
        except ValueError as e:
            msg = str(e)
            if "Fragment file" in msg and "does not exist" in msg:
                m = re.search(r"/([^/]+)/([^/]+)\.fragments\.tsv\.gz", msg)
                if not m:
                    raise
                bad_sample, bad_ct = m.group(1), m.group(2)
                logger.warning("Dropping cells of '%s' in sample '%s'; retrying.", bad_ct, bad_sample)
                with open(dropped_log, "a") as df:
                    df.write(f"{bad_sample}\t{bad_ct}\tmissing_fragments_file\n")
                mask = cell_data[args.variable].astype(str).map(_sanitize).eq(bad_ct)
                if mask.sum() == 0:
                    mask = cell_data[args.sample_id_col].astype(str).eq(bad_sample)
                    path_to_fragments.pop(bad_sample, None)
                cell_data = cell_data[~mask].copy()
                if len(cell_data) == 0:
                    raise RuntimeError("No cells remaining after filtering.")
                continue
            if "must have the same keys" in msg:
                remaining = set(cell_data[args.sample_id_col].astype(str).unique())
                for s in set(path_to_fragments) - remaining:
                    path_to_fragments.pop(s)
                continue
            raise

    with open(os.path.join(pb_outdir, "bw_paths.tsv"), "w") as f:
        for v in bw_paths:
            f.write(f"{v}\t{bw_paths[v]}\n")
    with open(os.path.join(pb_outdir, "bed_paths.tsv"), "w") as f:
        for v in bed_paths:
            f.write(f"{v}\t{bed_paths[v]}\n")

    # ── MACS2 peak calling ──────────────────────────────────────────────────
    logger.info("Calling peaks with MACS2...")
    genome_size = "mm" if args.species in ("mmusculus", "mouse") else "hs"
    narrow_peak_dict = peak_calling(
        macs_path="macs2",
        bed_paths=bed_paths,
        outdir=os.path.join(pb_outdir, "MACS"),
        genome_size=genome_size,
        n_cpu=args.n_cpu,
        input_format="BEDPE",
        shift=73,
        ext_size=146,
        keep_dup="all",
        q_value=0.05,
        skip_empty_peaks=True,
        _temp_dir=tempfile.gettempdir(),
    )

    # ── Consensus peaks ─────────────────────────────────────────────────────
    logger.info("Inferring consensus peaks...")
    peak_half_width = 250
    from pycisTopic.iterative_peak_calling import calculate_peaks_and_extend
    filtered_dict = {
        ct: peaks for ct, peaks in narrow_peak_dict.items()
        if len(calculate_peaks_and_extend(peaks, peak_half_width, chromsizes, path_to_blacklist)) > 0
    }
    logger.info("Kept %d/%d CTs with peaks after blacklist filter.", len(filtered_dict), len(narrow_peak_dict))
    consensus_peaks = get_consensus_peaks(
        filtered_dict,
        peak_half_width=peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=path_to_blacklist,
    )
    consensus_bed_path = os.path.join(pb_outdir, "consensus_regions.bed")
    consensus_peaks.to_bed(path=consensus_bed_path, keep=True, compression="infer", chain=False)
    logger.info("Consensus peaks written to %s", consensus_bed_path)

    # ── Per-sample pycistopic QC ────────────────────────────────────────────
    logger.info("Running pycistopic qc for %d samples...", len(path_to_fragments))
    sample_ids = list(path_to_fragments.keys())
    for sid in sample_ids:
        parquet_out = os.path.join(qc_outdir, f"{sid}.fragments_stats_per_cb.parquet")
        if os.path.exists(parquet_out):
            logger.info("  %s: QC parquet already exists, skipping.", sid)
            continue
        frag_path = path_to_fragments[sid]
        logger.info("  Running pycistopic qc for sample '%s': %s", sid, frag_path)
        subprocess.check_call([
            "pycistopic", "qc",
            "--fragments", frag_path,
            "--regions", consensus_bed_path,
            "--tss", tss_bed,
            "--output", os.path.join(qc_outdir, sid),
        ])

    # ── group_list.tsv ──────────────────────────────────────────────────────
    logger.info("Building group_list.tsv (min_cells=%d)...", args.min_cells)
    group_counts = (
        cell_data.groupby([args.cell_type_col, condition_col])
        .size()
        .reset_index(name="n_cells")
    )
    group_counts.columns = ["cell_type_safe", "condition", "n_cells"]
    # Require ALL conditions for a cell type to independently meet the threshold.
    # A cell type where one condition has 600 cells but the other has 200 is not
    # useful for differential comparison — drop the whole cell type in that case.
    ct_min = group_counts.groupby("cell_type_safe")["n_cells"].min()
    passing_cts = ct_min[ct_min >= args.min_cells].index
    passing = group_counts[group_counts["cell_type_safe"].isin(passing_cts)].copy()
    n_total_cts = group_counts["cell_type_safe"].nunique()
    logger.info(
        "%d / %d cell types pass min_cells=%d in ALL conditions (%d / %d CT×condition groups).",
        len(passing_cts), n_total_cts, args.min_cells,
        len(passing), len(group_counts),
    )
    group_list_path = os.path.join(outdir, "group_list.tsv")
    passing.to_csv(group_list_path, sep="\t", index=False)
    logger.info("group_list.tsv written (%d groups): %s", len(passing), group_list_path)

    # ── Write fragments_map.tsv to outdir (for downstream phases) ───────────
    frag_map_out = os.path.join(outdir, "fragments_map.tsv")
    pd.DataFrame(
        list(path_to_fragments.items()), columns=["sample_id", "fragments_path"]
    ).to_csv(frag_map_out, sep="\t", index=False)

    logger.info("Phase 1 complete.")


if __name__ == "__main__":
    main()
