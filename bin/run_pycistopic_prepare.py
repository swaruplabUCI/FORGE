#!/usr/bin/env python

import argparse
import os
import sys
import logging
import numpy as np
import pandas as pd

import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features,
)
from pycisTopic.gene_activity import get_gene_activity

import polars as pl
import pyranges as pr
import scanpy as sc
import pickle
import subprocess
import tempfile
import json

def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("pycisTopic_prepare")


def parse_args():
    p = argparse.ArgumentParser(
        description="Run pycisTopic preparation: pseudobulk, consensus peaks, cisTopic, DARs, gene activity."
    )

    # Legacy single-sample option (still supported, but we will prefer --fragments-map in the pipeline)
    p.add_argument(
        "--fragments",
        help="Path to fragments.tsv.gz (single-sample or merged ATAC). "
             "For multi-sample usage, prefer --fragments-map.",
    )

    # New multi-sample mapping: sample_id<TAB>fragments_path
    p.add_argument(
        "--fragments-map",
        help="TSV file with at least two columns: 'sample_id' and 'fragments_path'. "
             "Each row maps a sample_id to its ATAC fragments file.",
    )

    p.add_argument(
        "--cell-metadata",
        required=True,
        help="Cell metadata TSV (barcode index or 'barcode' column).",
    )
    p.add_argument(
        "--species",
        required=True,
        help="Species code for TSS (e.g., 'hsapiens').",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory (will be created).",
    )
    p.add_argument(
        "--sample-id-col",
        default="sample_id",
        help="Column in cell metadata with sample id.",
    )
    p.add_argument(
        "--cell-type-col",
        default="cell_type",
        help="Column in cell metadata with cell type annotation.",
    )
    p.add_argument(
        "--variable",
        default="cell_type",
        help="Column to group pseudobulk by (e.g., cell_type).",
    )
    p.add_argument(
        "--n-cpu",
        type=int,
        default=8,
        help="Number of CPUs to use for pycisTopic operations.",
    )
    p.add_argument(
        "--mallet-path",
        default="Mallet-202108/bin/mallet",
        help="Path to MALLET binary (if not installed, this script will attempt wget+tar).",
    )
    p.add_argument(
        "--topics",
        default="10,20,30,40",
        help="Comma-separated list of topic numbers to try (e.g., 10,20,30,40).",
    )
    p.add_argument(
        "--selected-topics",
        type=int,
        default=40,
        help="Number of topics to select as final model (must be in --topics).",
    )
    p.add_argument(
        "--do-gene-activity",
        action="store_true",
        help="If set, compute gene activity and save pycistopic_gene_activity.h5ad",
    )
    p.add_argument(
        "--mudata-stats",
        default=None,
        help="Path to mudata_stats.json produced by MuData creation. "
             "Used to restrict pycisTopic to samples present in the MuData.",
    )
    p.add_argument(
    "--blacklist-bed",
    default=None,
    help="Path to ENCODE blacklist BED (can be .bed or .bed.gz). If provided, overrides pycisTopic internal default.",
    )
    p.add_argument(
    "--gtf",
    default=None,
    help="Path to a Gencode GTF file. If provided, TSS BED is generated locally from the GTF "
         "instead of querying BioMart (avoids ArrowTypeError from polars/pyarrow incompatibility).",
    )
    p.add_argument(
        "--phase1-only",
        action="store_true",
        help="Stop after per-sample QC and write group_list.tsv. Skip CistopicObject "
             "creation, LDA, and region_sets (Phase 1 of the 3-phase pipeline).",
    )
    p.add_argument(
        "--condition-col",
        default="condition",
        help="Column in cell metadata with condition label, used for CT×condition "
             "fan-out in group_list.tsv. Falls back to 'all' if column is absent.",
    )
    p.add_argument(
        "--min-cells",
        type=int,
        default=200,
        help="Minimum cells per CT×condition group for inclusion in group_list.tsv.",
    )
    args = p.parse_args()

    # Sanity check: require at least one of fragments or fragments-map
    if args.fragments is None and args.fragments_map is None:
        p.error("You must provide either --fragments (single file) or --fragments-map (multi-sample).")

    return args

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def _write_group_list(cell_data, cell_type_col, condition_col, min_cells, outdir, logger):
    """Write group_list.tsv for Phase 2 fan-out (CT × condition, min_cells filtered)."""
    if condition_col not in cell_data.columns:
        logger.warning("[phase1] Condition column '%s' absent — grouping by cell type only.", condition_col)
        cell_data = cell_data.copy()
        cell_data["_cond"] = "all"
        condition_col = "_cond"

    grp = (
        cell_data.groupby([cell_type_col, condition_col])
        .size().reset_index(name="n_cells")
    )
    grp.columns = ["cell_type_safe", "condition", "n_cells"]
    ct_min = grp.groupby("cell_type_safe")["n_cells"].min()
    passing_cts = ct_min[ct_min >= min_cells].index
    passing = grp[grp["cell_type_safe"].isin(passing_cts)].copy()

    logger.info(
        "[phase1] %d / %d cell types pass min_cells=%d in ALL conditions (%d / %d groups)",
        len(passing_cts), grp["cell_type_safe"].nunique(), min_cells, len(passing), len(grp),
    )
    out = os.path.join(outdir, "group_list.tsv")
    passing.to_csv(out, sep="\t", index=False)
    logger.info("[phase1] group_list.tsv → %s", out)


def ensure_mallet(mallet_path, logger):
    # If mallet binary not present, download and extract in current dir
    if os.path.exists(mallet_path):
        return mallet_path
    logger.info(f"Requested MALLET path not found: {mallet_path}. Trying to download Mallet-202108...")
    try:
        subprocess.check_call(
            [
                "wget",
                "https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz",
            ]
        )
        subprocess.check_call(["tar", "-xf", "Mallet-202108-bin.tar.gz"])
    except Exception as e:
        logger.error(f"Failed to download or extract MALLET: {e}")
        raise
    new_path = os.path.join("Mallet-202108", "bin", "mallet")
    if not os.path.exists(new_path):
        raise FileNotFoundError(f"MALLET binary not found at {new_path} after extraction.")
    logger.info(f"MALLET available at: {new_path}")
    return new_path


_SPECIES_TO_UCSC = {
    "hsapiens": "hg38",
    "mmusculus": "mm10",
}
_SPECIES_TO_MACS = {
    "hsapiens": "hs",
    "mmusculus": "mm",
}
_SPECIES_TO_BIOMART = {
    "hsapiens": "hsapiens_gene_ensembl",
    "mmusculus": "mmusculus_gene_ensembl",
}


def download_chromsizes(species):
    ucsc_genome = _SPECIES_TO_UCSC.get(species, "hg38")
    chromsizes = pd.read_table(
        f"http://hgdownload.cse.ucsc.edu/goldenPath/{ucsc_genome}/bigZips/{ucsc_genome}.chrom.sizes",
        header=None,
        names=["Chromosome", "End"],
    )
    chromsizes.insert(1, "Start", 0)
    return chromsizes, ucsc_genome


def generate_tss_bed_from_gtf(gtf_path, output_bed, logger):
    """Parse a Gencode GTF to produce a TSS BED file, bypassing BioMart entirely.

    Generates the canonical 7-column format that pycisTopic expects
    (matching write_tss_annotation_to_bed / get_tss_annotation_from_ensembl):
        # Chromosome  Start  End  Gene  Score  Strand  Transcript_type
    """
    logger.info("Generating TSS BED from GTF: %s", gtf_path)
    import gzip

    opener = gzip.open if gtf_path.endswith(".gz") else open
    records = []
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based; BED is 0-based
            end = int(fields[4])
            strand = fields[6]
            # Extract gene_name and gene_type from attributes
            attrs = fields[8]
            gene_name = None
            gene_type = "unknown"
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1]
            if gene_name is None:
                continue
            # TSS is the start for + strand, end for - strand
            if strand == "+":
                tss = start
            else:
                tss = end - 1
            records.append((chrom, tss, tss + 1, gene_name, 0, strand, gene_type))

    tss_df = pd.DataFrame(records, columns=[
        "# Chromosome", "Start", "End", "Gene", "Score", "Strand", "Transcript_type"
    ])
    # Add 'chr' prefix if missing (Gencode uses 'chr' but just in case)
    if not tss_df["# Chromosome"].iloc[0].startswith("chr"):
        tss_df["# Chromosome"] = "chr" + tss_df["# Chromosome"]
    tss_df.drop_duplicates(inplace=True)
    tss_df.to_csv(output_bed, sep="\t", index=False)
    logger.info("Wrote %d TSS entries to %s", len(tss_df), output_bed)


def main():
    args = parse_args()
    logger = setup_logger()

    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # 1. Load cell metadata
    logger.info(f"Loading cell metadata from {args.cell_metadata}")
    cell_data = pd.read_table(args.cell_metadata, index_col=0, sep="\t")
    # BD Rhapsody barcodes are integers; cast index to str before pyCisTopic sees it.
    cell_data.index = cell_data.index.astype(str)

    if "barcode" in cell_data.columns:
        cell_data.index = cell_data["barcode"].astype(str)

    # Strip sample-prefix from barcodes so they match raw fragment-file barcodes.
    # The RNA pipeline prefixes obs_names with "<sample_id>:" during concatenation
    # (e.g. "10k_PBMC:AAACAGCCACCGGCTA-1"), but fragment files contain only the
    # raw barcode ("AAACAGCCACCGGCTA-1").
    raw_idx = cell_data.index.astype(str)
    if raw_idx.str.contains(":").any():
        stripped = raw_idx.str.split(":", n=1).str[-1]
        n_changed = (stripped != raw_idx).sum()
        logger.info(
            "Stripped sample-prefix from %d / %d barcodes "
            "(e.g. '%s' -> '%s').",
            n_changed,
            len(raw_idx),
            raw_idx[0],
            stripped[0],
        )
        cell_data.index = stripped

    # Allow sample_id_col to be either a column OR the index name
    cols = list(cell_data.columns)
    index_names = [cell_data.index.name] if cell_data.index.name is not None else []

    if args.sample_id_col not in cols + index_names:
        raise ValueError(
            f"sample-id column '{args.sample_id_col}' not found in cell metadata "
            f"(columns={cols}, index={index_names})"
        )
    if args.sample_id_col not in cell_data.columns and args.sample_id_col in index_names:
        cell_data[args.sample_id_col] = cell_data.index.astype(str)

    if args.variable not in cell_data.columns:
        raise ValueError(f"Grouping variable '{args.variable}' not found in cell metadata.")

    # ------------------------------------------------------------------
    # OPTIONAL: restrict to samples present in MuData (mudata_stats.json)
    # ------------------------------------------------------------------
    if args.mudata_stats is not None:
        logger.info(f"Loading MuData stats from {args.mudata_stats}")
        with open(args.mudata_stats, "r") as f:
            mudata_stats = json.load(f)

        mudata_samples = {s["sample_id"] for s in mudata_stats.get("samples", [])}
        logger.info(
            "MuData contains %d samples. Example: %s",
            len(mudata_samples),
            list(sorted(mudata_samples))[:5],
        )

        before_cells = cell_data.shape[0]
        cell_data = cell_data[
            cell_data[args.sample_id_col].astype(str).isin(mudata_samples)
        ].copy()
        after_cells = cell_data.shape[0]
        logger.info(
            "Filtered cell metadata using MuData sample list: %d -> %d cells.",
            before_cells,
            after_cells,
        )

    # -------------------------------------------------------------------------
    # Build mapping from sample_id -> fragments path
    # -------------------------------------------------------------------------
    if args.fragments_map is not None:
        logger.info(f"Loading fragments map from {args.fragments_map}")
        frag_map_df = pd.read_table(args.fragments_map, sep="\t")

        required_cols = {"sample_id", "fragments_path"}
        missing_cols = required_cols - set(frag_map_df.columns)
        if missing_cols:
            raise ValueError(
                f"Fragments map is missing required columns: {missing_cols}. "
                f"Found columns: {list(frag_map_df.columns)}"
            )

        path_to_fragments = dict(
            zip(
                frag_map_df["sample_id"].astype(str),
                frag_map_df["fragments_path"].astype(str),
            )
        )

        logger.info(
            f"Fragments map loaded for {len(path_to_fragments)} samples. "
            f"Example: {list(path_to_fragments.items())[:3]}"
        )

        # Subset to sample_ids that actually appear in cell_data, to avoid unused entries
        sample_ids_in_cells = set(cell_data[args.sample_id_col].astype(str))

        # Sanity check: every sample_id in cell_data (after MuData filtering) must appear in fragments_map
        missing_for_samples = sample_ids_in_cells - set(path_to_fragments.keys())
        if missing_for_samples:
            raise ValueError(
                "The following sample_ids appear in cell metadata (after MuData filtering) "
                "but not in fragments map: "
                f"{sorted(list(missing_for_samples))[:10]} (showing up to 10). "
                "Ensure the fragments_map.tsv covers all MuData samples."
            )

        # Restrict fragments mapping to the intersection:
        #   - samples in cell_data (already restricted to MuData), and
        #   - samples with fragments listed in fragments_map
        path_to_fragments = {
            sid: path_to_fragments[sid]
            for sid in sample_ids_in_cells
        }

        # Further filter: keep only sample_ids whose fragment files actually exist
        filtered_path_to_fragments = {}
        dropped_missing = []
        for sid, path in path_to_fragments.items():
            if not os.path.exists(path):
                dropped_missing.append((sid, path))
            else:
                filtered_path_to_fragments[sid] = path

        # Further filter: keep only sample_ids whose fragment files actually exist
        filtered_path_to_fragments = {}
        dropped_missing = []
        for sid, path in path_to_fragments.items():
            if not os.path.exists(path):
                dropped_missing.append((sid, path))
            else:
                filtered_path_to_fragments[sid] = path

        if dropped_missing:
            logger.warning(
                f"Dropping {len(dropped_missing)} sample_ids whose fragment files are missing. "
                f"Examples: {dropped_missing[:3]}"
            )

        if not filtered_path_to_fragments:
            raise RuntimeError(
                "After filtering, no valid fragment files remain for pycisTopic. "
                "Check your fragment directories and sample metadata."
            )

        path_to_fragments = filtered_path_to_fragments

        logger.info(
            f"Using {len(path_to_fragments)} samples with valid fragment files for pseudobulk."
        )

        # Also restrict cell_data to only those sample_ids with valid fragments
        valid_ids = set(path_to_fragments.keys())
        before_cells = cell_data.shape[0]
        cell_data = cell_data[cell_data[args.sample_id_col].astype(str).isin(valid_ids)].copy()
        after_cells = cell_data.shape[0]
        logger.info(
            f"Filtered cell metadata from {before_cells} to {after_cells} cells "
            f"matching {len(valid_ids)} samples with valid fragments."
        )

    else:
        # Legacy single-sample mode: use one fragments file for all cells
        logger.warning(
            "Running in single-fragments mode (no --fragments-map provided). "
            "All cells are assumed to come from the same fragments file."
        )
        # Use sample_id_col values to build a single-entry dict,
        # mapping the first sample_id to args.fragments
        first_sample_id = cell_data[args.sample_id_col].astype(str).iloc[0]
        path_to_fragments = {first_sample_id: args.fragments}

    # 2. Chrom sizes — species-aware
    chromsizes, ucsc_genome = download_chromsizes(args.species)
    logger.info(f"Downloaded chromsizes for {ucsc_genome} (species={args.species})")
    chromsizes.to_csv(os.path.join(outdir, f"{ucsc_genome}.chrom.sizes"), sep="\t", index=False)

    # 3. Export pseudobulk profiles
    pb_outdir = os.path.join(outdir, "consensus_peak_calling")
    ensure_dir(pb_outdir)
    ensure_dir(os.path.join(pb_outdir, "pseudobulk_bed_files"))
    ensure_dir(os.path.join(pb_outdir, "pseudobulk_bw_files"))

    dropped_log_path = os.path.join(pb_outdir, "dropped_groups.tsv")

    logger.info("Exporting pseudobulk profiles...")

    # We may need to iteratively drop problematic (sample_id, cell_type) groups
    # We log every dropped group to dropped_groups.tsv so users can inspect later.
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
                # If we reach this point, export_pseudobulk succeeded
                break    

            except ValueError as e:
                msg = str(e)
                # Typical message:
                #   "Fragment file /tmp/.../S52/L2_3_IT.fragments.tsv.gz does not exist."
                if "Fragment file" in msg and "does not exist" in msg:
                    logger.warning(
                        "export_pseudobulk encountered missing per-cell-type fragment file: %s",
                        msg,
                    )
                    import re    

                    m = re.search(r"/([^/]+)/([^/]+)\.fragments\.tsv\.gz", msg)
                    if not m:
                        logger.error(
                            "Could not parse (sample_id, cell_type_pattern) from error message; re-raising."
                        )
                        raise    

                    bad_sample, bad_celltype_pattern = m.group(1), m.group(2)

                    logger.warning(
                        "Will drop cells of cell type '%s' (sample=%s) and retry.",
                        bad_celltype_pattern,
                        bad_sample,
                    )

                    # Log dropped group to a TSV file for later inspection
                    with open(dropped_log_path, "a") as df:
                        df.write(
                            f"{bad_sample}\t{bad_celltype_pattern}\tmissing_fragments_file\n"
                        )

                    before = cell_data.shape[0]

                    # Drop only the problematic cell type, not the entire sample
                    # Convert bad_celltype_pattern back to original cell type name
                    # (underscores may replace spaces/special chars)
                    ct_col = args.variable
                    ct_values = cell_data[ct_col].astype(str)
                    # Try exact match on sanitized cell type name
                    import re as _re
                    def _sanitize(s):
                        s = _re.sub(r"\s+", "_", str(s))
                        s = _re.sub(r"[^A-Za-z0-9_.-]", "_", s)
                        s = _re.sub(r"_+", "_", s).strip("_")
                        return s
                    mask = ct_values.map(_sanitize).eq(bad_celltype_pattern)
                    if mask.sum() == 0:
                        # Fallback: drop entire sample if cell type can't be matched
                        logger.warning(
                            "Could not match cell type '%s'; dropping entire sample '%s'.",
                            bad_celltype_pattern, bad_sample,
                        )
                        mask = cell_data[args.sample_id_col].astype(str).eq(bad_sample)
                        if bad_sample in path_to_fragments:
                            path_to_fragments.pop(bad_sample, None)
                    cell_data = cell_data[~mask].copy()
                    after = cell_data.shape[0]
                    dropped = before - after

                    if dropped == 0:
                        logger.error(
                            "Failed to drop any cells for cell type '%s'; "
                            "cannot proceed safely. Re-raising the original error.",
                            bad_celltype_pattern,
                        )
                        raise    

                    logger.info(
                        "Filtered out %d cells (cell type '%s'). Remaining cells: %d",
                        dropped,
                        bad_celltype_pattern,
                        after,
                    )

                    if after == 0:
                        logger.error("No cells remaining after filtering; cannot proceed.")
                        raise

                    # FIX-98: Sync path_to_fragments with remaining samples in cell_data.
                    # After dropping cells, some samples may have zero cells left.
                    # split_fragment_files_by_cell_type requires exact key match.
                    remaining_samples = set(cell_data[args.sample_id_col].unique())
                    orphan_samples = set(path_to_fragments.keys()) - remaining_samples
                    for s in orphan_samples:
                        logger.info("FIX-98: Removing orphan sample '%s' from fragment map (no cells remain).", s)
                        path_to_fragments.pop(s)

                    # Loop back and retry export_pseudobulk with the reduced cell_data
                    continue

                # FIX-98: Catch key mismatch after cell filtering (e.g., "Unknown" drop
                # leaves samples in fragment map with no cells)
                if "must have the same keys" in msg:
                    logger.warning("FIX-98: Key mismatch between fragment map and cell data. Syncing...")
                    remaining_samples = set(cell_data[args.sample_id_col].unique())
                    orphan_samples = set(path_to_fragments.keys()) - remaining_samples
                    if orphan_samples:
                        for s in orphan_samples:
                            logger.info("FIX-98: Removing orphan sample '%s' from fragment map.", s)
                            path_to_fragments.pop(s)
                        continue
                    else:
                        # Keys match but still erroring — could be the reverse (cell_data has samples not in fragments)
                        logger.warning("FIX-98: No orphan samples found. Checking reverse direction...")
                        frag_samples = set(path_to_fragments.keys())
                        extra_in_cells = remaining_samples - frag_samples
                        if extra_in_cells:
                            logger.info("FIX-98: Dropping %d samples from cell_data that lack fragment files: %s",
                                       len(extra_in_cells), extra_in_cells)
                            cell_data = cell_data[cell_data[args.sample_id_col].isin(frag_samples)].copy()
                            continue
                        raise

                # Not a known recoverable error; re-raise
                raise

    # Write out mapping of group -> pseudobulk bigWig and BED
    with open(os.path.join(pb_outdir, "bw_paths.tsv"), "wt") as f:
        for v in bw_paths:
            _ = f.write(f"{v}\t{bw_paths[v]}\n")

    with open(os.path.join(pb_outdir, "bed_paths.tsv"), "wt") as f:
        for v in bed_paths:
            _ = f.write(f"{v}\t{bed_paths[v]}\n")

    # 4. Peak calling with MACS2 on pseudobulk
    logger.info("Calling peaks on pseudobulk with MACS2...")
    from pycisTopic.pseudobulk_peak_calling import peak_calling

    macs_path = "macs2"
    narrow_peak_dict = peak_calling(
        macs_path=macs_path,
        bed_paths=bed_paths,
        outdir=os.path.join(pb_outdir, "MACS"),
        genome_size=_SPECIES_TO_MACS.get(args.species, "hs"),
        n_cpu=args.n_cpu,
        input_format="BEDPE",
        shift=73,
        ext_size=146,
        keep_dup="all",
        q_value=0.05,
        _temp_dir=tempfile.gettempdir(),
    )

    # 5. Consensus peaks using TGCA iterative filtering
    logger.info("Inferring consensus peaks...")
    peak_half_width = 250
    path_to_blacklist = args.blacklist_bed or os.path.join(os.path.dirname(pycisTopic.__file__), "blacklist", "hg38-blacklist.v2.bed")
    # Pre-filter: drop cell types whose peaks vanish after extension/blacklist filtering.
    # pycisTopic's iterative_peak_filtering crashes on empty PyRanges (no .Count column).
    from pycisTopic.iterative_peak_calling import calculate_peaks_and_extend
    import pyranges as _pr
    filtered_dict = {}
    for ct, peaks in narrow_peak_dict.items():
        extended = calculate_peaks_and_extend(peaks, peak_half_width, chromsizes, path_to_blacklist)
        if len(extended) > 0:
            filtered_dict[ct] = peaks
        else:
            logger.warning(f"Dropping '{ct}' — 0 peaks after extension/blacklist filtering")
    logger.info(f"Kept {len(filtered_dict)}/{len(narrow_peak_dict)} cell types with peaks")
    consensus_peaks = get_consensus_peaks(
        filtered_dict,
        peak_half_width=peak_half_width,
        chromsizes=chromsizes,
        path_to_blacklist=path_to_blacklist,
    )
    consensus_bed_path = os.path.join(pb_outdir, "consensus_regions.bed")
    consensus_peaks.to_bed(
        path=consensus_bed_path,
        keep=True,
        compression="infer",
        chain=False,
    )

    # 6. QC & create cisTopic object
    qc_outdir = os.path.join(outdir, "qc")
    ensure_dir(qc_outdir)

    # TSS annotation
    tss_bed = os.path.join(qc_outdir, "tss.bed")
    if not os.path.exists(tss_bed):
        if args.gtf:
            generate_tss_bed_from_gtf(args.gtf, tss_bed, logger)
        else:
            logger.info(f"Generating TSS annotation via pycistopic tss get_tss for {args.species}")
            cmd = [
                "pycistopic",
                "tss",
                "get_tss",
                "--output",
                tss_bed,
                "--name",
                _SPECIES_TO_BIOMART.get(args.species, "hsapiens_gene_ensembl"),
                "--to-chrom-source",
                "ucsc",
                "--ucsc",
                ucsc_genome,
            ]
            subprocess.check_call(cmd)

    # Barcode-level QC via CLI — iterate over ALL samples
    from pycisTopic.qc import get_barcodes_passing_qc_for_sample

    sample_ids = list(path_to_fragments.keys())
    logger.info(f"Running barcode-level QC for {len(sample_ids)} samples...")

    cistopic_objects = []
    total_cells = 0
    for sid in sample_ids:
        frag_path = path_to_fragments[sid]
        logger.info(f"Running pycistopic qc on fragments (barcode-level metrics) using sample '{sid}': {frag_path}")

        # QC output per sample to avoid overwriting
        sample_qc_prefix = os.path.join(qc_outdir, sid)
        cmd_qc = [
            "pycistopic",
            "qc",
            "--fragments",
            frag_path,
            "--regions",
            consensus_bed_path,
            "--tss",
            tss_bed,
            "--output",
            sample_qc_prefix,
        ]
        subprocess.check_call(cmd_qc)

        bc_passing, thresholds = get_barcodes_passing_qc_for_sample(
            sample_id=sid,
            pycistopic_qc_output_dir=qc_outdir,
            unique_fragments_threshold=None,
            tss_enrichment_threshold=None,
            frip_threshold=0,
            use_automatic_thresholds=True,
        )

        logger.info(f"  {sid}: {len(bc_passing)} barcodes passing QC")

        if len(bc_passing) == 0:
            logger.warning(f"  {sid}: 0 barcodes passing QC — skipping")
            continue

        # Phase 1-only: QC done; skip CistopicObject creation
        if args.phase1_only:
            continue

        sample_metrics = (
            pl.read_parquet(os.path.join(qc_outdir, f"{sid}.fragments_stats_per_cb.parquet"))
            .to_pandas()
            .set_index("CB")
            .loc[bc_passing]
        )

        obj = create_cistopic_object_from_fragments(
            path_to_fragments=frag_path,
            path_to_regions=consensus_bed_path,
            path_to_blacklist=path_to_blacklist,
            metrics=sample_metrics,
            valid_bc=bc_passing,
            n_cpu=1,
            project=sid,
            split_pattern="-",
        )
        cistopic_objects.append(obj)
        total_cells += len(bc_passing)
        logger.info(f"  {sid}: created cisTopic object with {len(obj.cell_names)} cells")

    # Phase 1-only exit: write group_list + supporting files then stop.
    if args.phase1_only:
        _write_group_list(cell_data, args.cell_type_col, args.condition_col,
                          args.min_cells, outdir, logger)
        pd.DataFrame(
            list(path_to_fragments.items()), columns=["sample_id", "fragments_path"]
        ).to_csv(os.path.join(outdir, "fragments_map.tsv"), sep="\t", index=False)
        bl_dest = os.path.join(outdir, "blacklist.bed")
        if path_to_blacklist and os.path.abspath(path_to_blacklist) != os.path.abspath(bl_dest):
            import shutil
            shutil.copy(path_to_blacklist, bl_dest)
        meta_dest = os.path.join(outdir, "cell_metadata_for_pycistopic.safe.tsv")
        cell_data.to_csv(meta_dest, sep="\t", index=True)
        logger.info("[phase1] cell_metadata_for_pycistopic.safe.tsv → %s", meta_dest)
        logger.info("[phase1] Done. Outputs in %s", outdir)
        sys.exit(0)

    if len(cistopic_objects) == 0:
        raise ValueError("No samples produced cells passing QC — cannot proceed.")

    # Merge per-sample cisTopic objects
    if len(cistopic_objects) == 1:
        cistopic_obj = cistopic_objects[0]
        logger.info(f"Single sample — {len(cistopic_obj.cell_names)} cells total")
    else:
        cistopic_obj = cistopic_objects[0].merge(
            cistopic_objects[1:],
            is_acc=1,
            project="merged",
            copy=True,
            split_pattern="-",
        )
        logger.info(f"Merged {len(cistopic_objects)} samples — {len(cistopic_obj.cell_names)} cells total")

    # Merge cell type annotations into cistopic_obj.cell_data
    # create_cistopic_object_from_fragments builds cell_data from fragment metrics
    # but does not include user-provided annotations (e.g. cell_type_safe).
    # pycisTopic may prefix barcodes with "project___" (e.g. "10k_PBMC___AAACAGCCACCGGCTA-1")
    # while cell_data uses stripped raw barcodes ("AAACAGCCACCGGCTA-1").
    logger.info("Merging cell type annotations into cisTopic object...")
    cto_barcodes = cistopic_obj.cell_data.index.astype(str)
    meta_barcodes = cell_data.index.astype(str)
    shared = cto_barcodes.intersection(meta_barcodes)
    logger.info(f"  cistopic barcodes: {len(cto_barcodes)}, metadata barcodes: {len(meta_barcodes)}, shared: {len(shared)}")
    logger.info(f"  cistopic examples: {list(cto_barcodes[:3])}")
    logger.info(f"  metadata examples: {list(meta_barcodes[:3])}")

    # If no overlap, try stripping project suffix from cistopic barcodes.
    # pycisTopic's create_cistopic_object (tag_cells=True) appends
    # split_pattern + project to each barcode, e.g.:
    #   AAACAGCCACCGGCTA-1  ->  AAACAGCCACCGGCTA-1-10k_PBMC
    # when split_pattern="-" and project="10k_PBMC".
    if len(shared) == 0 and len(cto_barcodes) > 0:
        sample_cto = str(cto_barcodes[0])
        sample_meta = str(meta_barcodes[0])
        logger.info(f"  Trying to reconcile barcode formats...")

        # Try stripping known project suffix pattern: barcode + sep + project
        # The project name used was sample_id; detect the suffix
        project_suffix = None
        for sep in ["-", "___", "__", "_"]:
            suffix = sep + sample_id
            if sample_cto.endswith(suffix):
                project_suffix = suffix
                break

        if project_suffix:
            logger.info(f"  Detected project suffix '{project_suffix}' in cistopic barcodes")
            cto_stripped = cto_barcodes.str[:-len(project_suffix)]
            shared = cto_stripped.intersection(meta_barcodes)
            logger.info(f"  After stripping: shared={len(shared)} (e.g. '{sample_cto}' -> '{cto_stripped[0]}')")
            if len(shared) > 0:
                # Build lookup from stripped barcode -> original cistopic barcode
                bc_map = pd.Series(cto_barcodes.values, index=cto_stripped.values)
                # Reindex cell_data to match cistopic's original index
                cell_data_reindexed = cell_data.copy()
                cell_data_reindexed.index = cell_data_reindexed.index.map(
                    lambda x: bc_map.get(x, x)
                )
                cell_data = cell_data_reindexed
                meta_barcodes = cell_data.index.astype(str)
                shared = cto_barcodes.intersection(meta_barcodes)
                logger.info(f"  After reindexing metadata: shared={len(shared)}")
        else:
            logger.warning(f"  Could not detect project suffix pattern in cistopic barcode '{sample_cto}'")

    cols_to_add = [c for c in cell_data.columns if c not in cistopic_obj.cell_data.columns]
    if cols_to_add and len(shared) > 0:
        cistopic_obj.cell_data = cistopic_obj.cell_data.join(
            cell_data[cols_to_add], how="left"
        )
        logger.info(f"  Added columns to cistopic_obj.cell_data: {cols_to_add}")
        # Verify the target variable column has non-null values
        if args.variable in cols_to_add:
            n_valid = cistopic_obj.cell_data[args.variable].notna().sum()
            logger.info(f"  Variable '{args.variable}': {n_valid}/{len(cistopic_obj.cell_data)} non-null")
    else:
        logger.warning(f"  No new columns to add (cols_to_add={cols_to_add}, shared={len(shared)})")
        if len(shared) == 0 and cols_to_add:
            raise ValueError(
                f"Cannot merge cell annotations: 0 barcode overlap between cisTopic object "
                f"and metadata. cisTopic barcodes look like '{cto_barcodes[0]}', "
                f"metadata barcodes look like '{meta_barcodes[0]}'. "
                f"Check barcode prefix handling."
            )

    # 7. Topic modeling via MALLET
    logger.info("Running MALLET LDA models...")
    mallet_path = ensure_mallet(args.mallet_path, logger)
    os.environ["MALLET_MEMORY"] = "200G"

    n_topics = [int(x) for x in args.topics.split(",") if x.strip()]
    mallet_tmp = os.path.join(outdir, "mallet_tmp")
    mallet_save = os.path.join(outdir, "mallet_models")
    os.makedirs(mallet_tmp, exist_ok=True)
    os.makedirs(mallet_save, exist_ok=True)
    models = run_cgs_models_mallet(
        cistopic_obj,
        n_topics=n_topics,
        n_cpu=args.n_cpu,
        n_iter=500,
        random_state=555,
        alpha=50,
        alpha_by_topic=True,
        eta=0.1,
        eta_by_topic=False,
        tmp_path=mallet_tmp,
        save_path=mallet_save,
        mallet_path=mallet_path,
    )

    logger.info("Evaluating models...")
    model = evaluate_models(
        models,
        select_model=args.selected_topics,
        return_model=True,
    )
    cistopic_obj.add_LDA_model(model)

    # 8. Clustering/UMAP optional (you can extend later if you want)
    # from pycisTopic.clust_vis import find_clusters, run_umap
    # find_clusters(cistopic_obj, target='cell', k=10, res=[0.6, 1.2, 3], prefix='pycisTopic_', scale=True)

    # 9. Topic binarization
    logger.info("Binarizing topics...")
    region_bin_topics_top_3k = binarize_topics(cistopic_obj, method="ntop", ntop=3000, plot=False)
    region_bin_topics_otsu = binarize_topics(cistopic_obj, method="otsu", plot=False)

    # 10. Impute accessibility & DARs
    logger.info("Imputing accessibility...")
    imputed_acc_obj = impute_accessibility(cistopic_obj, scale_factor=10**6)
    normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

    logger.info("Identifying highly variable regions...")
    variable_regions = find_highly_variable_features(
        normalized_imputed_acc_obj,
        min_disp=0.05,
        min_mean=0.0125,
        max_mean=3,
        max_disp=np.inf,
        n_bins=20,
        n_top_features=None,
        plot=False,
    )

    logger.info("Identifying DARs by cell type...")
    markers_dict = find_diff_features(
        cistopic_obj,
        imputed_acc_obj,
        variable=args.cell_type_col,
        var_features=variable_regions,
        contrasts=None,
        adjpval_thr=0.05,
        log2fc_thr=np.log2(1.5),
        n_cpu=args.n_cpu,
        _temp_dir=tempfile.gettempdir(),
        split_pattern="-",
    )

    # 11. Save region sets
    logger.info("Saving region sets as BED files...")
    region_sets_dir = os.path.join(outdir, "region_sets")
    ensure_dir(region_sets_dir)
    topics_otsu_dir = os.path.join(region_sets_dir, "Topics_otsu")
    topics_top_dir = os.path.join(region_sets_dir, "Topics_top_3k")
    dars_dir = os.path.join(region_sets_dir, "DARs_cell_type")
    ensure_dir(topics_otsu_dir)
    ensure_dir(topics_top_dir)
    ensure_dir(dars_dir)

    for topic in region_bin_topics_otsu:
        region_names_to_coordinates(region_bin_topics_otsu[topic].index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(topics_otsu_dir, f"{topic}.bed"),
            sep="\t",
            header=False,
            index=False,
        )

    for topic in region_bin_topics_top_3k:
        region_names_to_coordinates(region_bin_topics_top_3k[topic].index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(topics_top_dir, f"{topic}.bed"),
            sep="\t",
            header=False,
            index=False,
        )

    for cell_type in markers_dict:
        dar_df = markers_dict[cell_type]
        if len(dar_df) == 0:
            logger.info(f"Skipping empty DAR set for {cell_type}")
            continue
        region_names_to_coordinates(dar_df.index).sort_values(
            ["Chromosome", "Start", "End"]
        ).to_csv(
            os.path.join(dars_dir, f"{cell_type}.bed"),
            sep="\t",
            header=False,
            index=False,
        )

    # 12. Gene activity (optional)
    if args.do_gene_activity:
        logger.info("Computing gene activity...")
        # chromsizes again as PyRanges
        chromsizes_pr = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

        pr_annotation = pd.read_table(tss_bed).rename(
            {"Name": "Gene", "# Chromosome": "Chromosome"}, axis=1
        )
        pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
        pr_annotation = pr.PyRanges(pr_annotation)

        gene_act, weights = get_gene_activity(
            imputed_acc_obj,
            pr_annotation,
            chromsizes_pr,
            use_gene_boundaries=True,
            upstream=[1000, 100000],
            downstream=[1000, 100000],
            distance_weight=True,
            decay_rate=1,
            extend_gene_body_upstream=10000,
            extend_gene_body_downstream=500,
            gene_size_weight=False,
            gene_size_scale_factor="median",
            remove_promoters=False,
            average_scores=True,
            scale_factor=1,
            extend_tss=[10, 10],
            gini_weight=True,
            return_weights=True,
            project="Gene_activity",
        )

        # Save gene activity as AnnData
        adata = sc.AnnData(
            X=gene_act.mtx.T,
            obs=pd.DataFrame(index=gene_act.cell_names),
            var=pd.DataFrame(index=gene_act.feature_names),
        )
        adata.write_h5ad(os.path.join(outdir, "pycistopic_gene_activity.h5ad"))

    # 13. Save cisTopic object
    logger.info("Saving cisTopic object...")
    with open(os.path.join(outdir, "cistopic_obj.pkl"), "wb") as f:
        pickle.dump(cistopic_obj, f)

    logger.info("pycisTopic preparation complete.")


if __name__ == "__main__":
    main()
