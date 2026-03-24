#!/usr/bin/env python

import argparse
import os
import sys
import logging
import subprocess
import textwrap


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger("scenicplus_snakemake")


def parse_args():
    p = argparse.ArgumentParser(description="Run SCENIC+ Snakemake pipeline.")
    p.add_argument("--cistopic", required=True, help="Path to cistopic_obj.pkl from pycisTopic.")
    p.add_argument("--rna", required=True, help="Path to RNA AnnData (h5ad).")
    p.add_argument("--region-sets", required=True, help="Path to region_sets directory (Topics_otsu, Topics_top_3k, DARs_cell_type).")
    p.add_argument("--ctx-rankings", required=True, help="cisTarget rankings feather.")
    p.add_argument("--ctx-scores", required=True, help="cisTarget scores feather.")
    p.add_argument("--motif-annotations", required=True, help="motif annotation table (.tbl).")
    p.add_argument("--outdir", required=True, help="Output directory (this script will create outs/ under it).")
    p.add_argument("--n-cpu", type=int, default=20, help="Number of cores for snakemake.")
    p.add_argument("--species", default="hsapiens", help="Species name for genome annotation (e.g., hsapiens).")
    p.add_argument("--biomart-host", default="http://ensembl.org/", help="Biomart host for genome annotations.")
    p.add_argument("--gtf", default=None, help="Local GTF file for offline genome annotation (bypasses BioMart).")
    p.add_argument("--fai", default=None, help="Local .fai index for offline chromsizes (bypasses BioMart).")
    p.add_argument(
        "--bc-transform",
        default=None,
        help="Python lambda string for barcode transform (e.g. \"lambda x: '10k_PBMC:' + x\"). "
             "If not provided, uses the SCENIC+ tutorial default.",
    )
    return p.parse_args()


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def generate_genome_annotation_from_gtf(gtf_path, out_path, logger):
    """Parse a Gencode GTF to produce genome_annotation.tsv matching BioMart format.

    SCENIC+ gene_search_space.py requires 7 columns:
      Chromosome, Start, End, Strand, Gene, Transcription_Start_Site, Transcript_type
    Strand must be '+' or '-' (not 1/-1).
    Only protein_coding genes are included.
    """
    logger.info(f"Generating genome_annotation.tsv from {gtf_path}")
    genes = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand_char = fields[6]  # '+' or '-'
            attrs = fields[8]
            gene_name = None
            gene_type = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    gene_name = attr.split('"')[1]
                elif attr.startswith("gene_type"):
                    gene_type = attr.split('"')[1]
            if gene_name is None:
                continue
            # Only keep protein_coding genes
            if gene_type != "protein_coding":
                continue
            # TSS: start of transcription depending on strand
            tss = start if strand_char == "+" else end
            # FIX-82: Keep "chr" prefix — pycisTopic regions use UCSC-style "chr1"
            # naming. Stripping the prefix caused zero-overlap join in get_search_space.
            genes[gene_name] = (chrom, start, end, strand_char, gene_name, tss, "protein_coding")
    with open(out_path, "w") as fh:
        fh.write("Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site\tTranscript_type\n")
        for g in sorted(genes.values()):
            fh.write(f"{g[0]}\t{g[1]}\t{g[2]}\t{g[3]}\t{g[4]}\t{g[5]}\t{g[6]}\n")
    logger.info(f"Wrote {len(genes)} protein_coding genes to {out_path}")


def generate_chromsizes_from_fai(fai_path, out_path, logger):
    """Parse a .fai index to produce chromsizes.tsv matching BioMart format."""
    logger.info(f"Generating chromsizes.tsv from {fai_path}")
    rows = []
    with open(fai_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            length = int(fields[1])
            # FIX-82: Keep "chr" prefix to match pycisTopic region naming
            rows.append((chrom, 0, length))
    with open(out_path, "w") as fh:
        fh.write("Chromosome\tStart\tEnd\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
    logger.info(f"Wrote {len(rows)} chromosomes to {out_path}")


def write_config_yaml(cfg_path, input_data, output_data, params_general, params_data, params_motif, params_infer, bc_transform_func=None):
    if bc_transform_func is None:
        bc_transform_func = "lambda x: f'{x}-10x_multiome_brain'"
    cfg = textwrap.dedent(
        f"""
        input_data:
          cisTopic_obj_fname: "{input_data['cistopic_obj']}"
          GEX_anndata_fname: "{input_data['rna']}"
          region_set_folder: "{input_data['region_sets']}"
          ctx_db_fname: "{input_data['ctx_rankings']}"
          dem_db_fname: "{input_data['ctx_scores']}"
          path_to_motif_annotations: "{input_data['motif_annotations']}"

        output_data:
          combined_GEX_ACC_mudata: "{output_data['combined_GEX_ACC_mudata']}"
          dem_result_fname: "{output_data['dem_result_fname']}"
          ctx_result_fname: "{output_data['ctx_result_fname']}"
          output_fname_dem_html: "{output_data['output_fname_dem_html']}"
          output_fname_ctx_html: "{output_data['output_fname_ctx_html']}"
          cistromes_direct: "{output_data['cistromes_direct']}"
          cistromes_extended: "{output_data['cistromes_extended']}"
          tf_names: "{output_data['tf_names']}"
          genome_annotation: "{output_data['genome_annotation']}"
          chromsizes: "{output_data['chromsizes']}"
          search_space: "{output_data['search_space']}"
          tf_to_gene_adjacencies: "{output_data['tf_to_gene_adjacencies']}"
          region_to_gene_adjacencies: "{output_data['region_to_gene_adjacencies']}"
          eRegulons_direct: "{output_data['eRegulons_direct']}"
          eRegulons_extended: "{output_data['eRegulons_extended']}"
          AUCell_direct: "{output_data['AUCell_direct']}"
          AUCell_extended: "{output_data['AUCell_extended']}"
          scplus_mdata: "{output_data['scplus_mdata']}"

        params_general:
          temp_dir: "{params_general['temp_dir']}"
          n_cpu: {params_general['n_cpu']}
          seed: {params_general['seed']}

        params_data_preparation:
          bc_transform_func: "\\"{bc_transform_func}\\""
          is_multiome: True
          key_to_group_by: ""
          nr_cells_per_metacells: 10
          direct_annotation: "Direct_annot"
          extended_annotation: "Orthology_annot"
          species: "{params_data['species']}"
          biomart_host: "{params_data['biomart_host']}"
          search_space_upstream: "1000 150000"
          search_space_downstream: "1000 150000"
          search_space_extend_tss: "10 10"

        params_motif_enrichment:
          species: "homo_sapiens"
          annotation_version: "v10nr_clust"
          motif_similarity_fdr: 0.001
          orthologous_identity_threshold: 0.0
          annotations_to_use: "Direct_annot Orthology_annot"
          fraction_overlap_w_dem_database: 0.4
          dem_max_bg_regions: 500
          dem_balance_number_of_promoters: True
          dem_promoter_space: 1000
          dem_adj_pval_thr: 0.05
          dem_log2fc_thr: 1.0
          dem_mean_fg_thr: 0.0
          dem_motif_hit_thr: 3.0
          fraction_overlap_w_ctx_database: 0.4
          ctx_auc_threshold: 0.005
          ctx_nes_threshold: 3.0
          ctx_rank_threshold: 0.05

        params_inference:
          tf_to_gene_importance_method: "GBM"
          region_to_gene_importance_method: "GBM"
          region_to_gene_correlation_method: "SR"
          order_regions_to_genes_by: "importance"
          order_TFs_to_genes_by: "importance"
          gsea_n_perm: 1000
          quantile_thresholds_region_to_gene: "0.85 0.90 0.95"
          top_n_regionTogenes_per_gene: "5 10 15"
          top_n_regionTogenes_per_region: ""
          min_regions_per_gene: 0
          rho_threshold: 0.05
          min_target_genes: 10
        """
    ).strip() + "\n"

    with open(cfg_path, "w") as f:
        f.write(cfg)


def preprocess_rna_raw(rna_path, logger):
    """Ensure RNA AnnData has .raw set from .layers['counts'].

    SCENIC+ requires raw counts in adata.raw. Our integrated RNA h5ad
    has counts in .layers['counts'] but .raw is not set.
    Returns the path to the (possibly rewritten) h5ad.
    """
    import scanpy as sc

    logger.info("Checking RNA AnnData for .raw layer: %s", rna_path)
    adata = sc.read_h5ad(rna_path)

    if adata.raw is not None:
        logger.info("RNA AnnData already has .raw set — skipping preprocessing.")
        return rna_path

    if "counts" in adata.layers:
        logger.info("Setting .raw from .layers['counts']")
        import anndata
        adata.raw = anndata.AnnData(
            X=adata.layers["counts"],
            var=adata.var,
            obs=adata.obs,
        )
    else:
        logger.warning(
            "No .raw and no 'counts' layer found — SCENIC+ may fail. "
            "Available layers: %s", list(adata.layers.keys())
        )
        return rna_path

    preprocessed_path = rna_path.replace(".h5ad", "_with_raw.h5ad")
    if preprocessed_path == rna_path:
        preprocessed_path = rna_path + ".with_raw.h5ad"
    adata.write_h5ad(preprocessed_path)
    logger.info("Wrote preprocessed RNA with .raw to %s", preprocessed_path)
    return preprocessed_path


def main():
    args = parse_args()
    logger = setup_logger()

    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # Preprocess RNA: ensure .raw is set from counts layer
    rna_path = preprocess_rna_raw(os.path.abspath(args.rna), logger)

    # Paths relative to outdir
    snakemake_dir = os.path.join(outdir, "scplus_pipeline")
    ensure_dir(snakemake_dir)

    # Initialize snakemake
    logger.info("Initializing SCENIC+ Snakemake pipeline...")
    subprocess.check_call(
        ["scenicplus", "init_snakemake", "--out_dir", snakemake_dir]
    )

    config_dir = os.path.join(snakemake_dir, "Snakemake", "config")
    workflow_dir = os.path.join(snakemake_dir, "Snakemake", "workflow")
    cfg_path = os.path.join(config_dir, "config.yaml")

    # NOTE: The Snakefile uses "search_spance" (upstream typo) which matches
    # the installed CLI subcommand. Do NOT rename it — the CLI rejects "search_space".

    outs_dir = os.path.join(outdir, "outs")
    ensure_dir(outs_dir)

    # Build config dictionaries
    input_data = {
        "cistopic_obj": os.path.abspath(args.cistopic),
        "rna": rna_path,
        "region_sets": os.path.abspath(args.region_sets),
        "ctx_rankings": os.path.abspath(args.ctx_rankings),
        "ctx_scores": os.path.abspath(args.ctx_scores),
        "motif_annotations": os.path.abspath(args.motif_annotations),
    }

    output_data = {
        "combined_GEX_ACC_mudata": os.path.join(outs_dir, "ACC_GEX.h5mu"),
        "dem_result_fname": os.path.join(outs_dir, "dem_results.hdf5"),
        "ctx_result_fname": os.path.join(outs_dir, "ctx_results.hdf5"),
        "output_fname_dem_html": os.path.join(outs_dir, "dem_results.html"),
        "output_fname_ctx_html": os.path.join(outs_dir, "ctx_results.html"),
        "cistromes_direct": os.path.join(outs_dir, "cistromes_direct.h5ad"),
        "cistromes_extended": os.path.join(outs_dir, "cistromes_extended.h5ad"),
        "tf_names": os.path.join(outs_dir, "tf_names.txt"),
        "genome_annotation": os.path.join(outs_dir, "genome_annotation.tsv"),
        "chromsizes": os.path.join(outs_dir, "chromsizes.tsv"),
        "search_space": os.path.join(outs_dir, "search_space.tsv"),
        "tf_to_gene_adjacencies": os.path.join(outs_dir, "tf_to_gene_adj.tsv"),
        "region_to_gene_adjacencies": os.path.join(outs_dir, "region_to_gene_adj.tsv"),
        "eRegulons_direct": os.path.join(outs_dir, "eRegulon_direct.tsv"),
        "eRegulons_extended": os.path.join(outs_dir, "eRegulons_extended.tsv"),
        "AUCell_direct": os.path.join(outs_dir, "AUCell_direct.h5mu"),
        "AUCell_extended": os.path.join(outs_dir, "AUCell_extended.h5mu"),
        "scplus_mdata": os.path.join(outs_dir, "scplusmdata.h5mu"),
    }

    params_general = {
        "temp_dir": os.path.join(outdir, "tmp"),
        "n_cpu": args.n_cpu,
        "seed": 666,
    }

    params_data = {
        "species": args.species,
        "biomart_host": args.biomart_host,
    }

    params_motif = {}  # fields are embedded directly in YAML
    params_infer = {}  # same

    ensure_dir(params_general["temp_dir"])

    logger.info(f"Writing SCENIC+ config to {cfg_path}")
    write_config_yaml(
        cfg_path,
        input_data,
        output_data,
        params_general,
        params_data,
        params_motif,
        params_infer,
        bc_transform_func=args.bc_transform,
    )

    # Pre-generate genome annotation + chromsizes from local files (bypass BioMart)
    if args.gtf and args.fai:
        genome_annot_path = output_data["genome_annotation"]
        chromsizes_path = output_data["chromsizes"]
        # FIX-82: Always regenerate to ensure chr prefix matches region naming
        generate_genome_annotation_from_gtf(args.gtf, genome_annot_path, logger)
        generate_chromsizes_from_fai(args.fai, chromsizes_path, logger)

    # Remove empty BED files — pyranges.read_bed crashes on 0-byte files
    region_sets_abs = os.path.abspath(args.region_sets)
    empty_removed = 0
    for root, _dirs, files in os.walk(region_sets_abs):
        for fname in files:
            fpath = os.path.join(root, fname)
            if fname.endswith(".bed") and os.path.getsize(fpath) == 0:
                os.remove(fpath)
                empty_removed += 1
    if empty_removed:
        logger.info(f"Removed {empty_removed} empty BED file(s) from region_sets/")

    # FIX-59: The set.union() bug in cistarget_wrangling.py is fixed via
    # the local source overlay at tools/scenicplus-main/src/ which is
    # prepended to PYTHONPATH by the Nextflow module script block.

    # Run snakemake
    logger.info("Running snakemake...")
    cmd = [
        "snakemake",
        "--cores",
        str(args.n_cpu),
    ]
    subprocess.check_call(cmd, cwd=os.path.join(snakemake_dir, "Snakemake"))

    logger.info("SCENIC+ Snakemake pipeline finished successfully.")


if __name__ == "__main__":
    main()
