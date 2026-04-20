#!/usr/bin/env Rscript
# Cicero Phase 3 — step 2 of 3: per-chromosome generate_cicero_models + assemble_connections
#
# Runs once per chromosome in parallel (Nextflow fan-out). Uses the pinned
# global distance_parameter from step 1 to remain bit-identical to the
# whole-genome run_cicero() baseline (Pearson rho = 1.0 validated by
# tests/cicero_parallel_test/global_dp_20260418_023059/).
#
# Output: conns_{chrom}.tsv.gz with columns Peak1, Peak2, coaccess

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(monocle3)
  library(cicero)
})

option_list <- list(
  make_option(c("--chrom"), type = "character", help = "Chromosome name (e.g., chr1)"),
  make_option(c("--cicero_cds"), type = "character", help = "Shared CiceroCellDataSet .rds from step 1"),
  make_option(c("--gene_annotation"), type = "character", help = "Shared gene_annotation .rds from step 1"),
  make_option(c("--dp"), type = "double", help = "Global distance_parameter scalar from step 1"),
  make_option(c("--window_bp"), type = "integer", default = 500000, help = "Window (bp) for generate_cicero_models"),
  make_option(c("--s"), type = "double", default = 0.75, help = "graphical lasso s parameter"),
  make_option(c("--max_elements"), type = "integer", default = 200, help = "max_elements for generate_cicero_models"),
  make_option(c("--outdir"), type = "character", default = ".", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$chrom) || is.null(opt$cicero_cds) || is.null(opt$gene_annotation) || is.null(opt$dp)) {
  stop("Required: --chrom, --cicero_cds, --gene_annotation, --dp", call. = FALSE)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$outdir)

cat(sprintf("== Cicero Phase 3 step 2: chrom=%s, dp=%.10f ==\n", opt$chrom, opt$dp))
cat(sprintf("  window=%d, s=%.2f, max_elements=%d\n",
            opt$window_bp, opt$s, opt$max_elements))

cat("Loading shared cicero_cds + gene_annotation...\n")
t0 <- Sys.time()
cicero_cds <- readRDS(opt$cicero_cds)
gene_ann   <- readRDS(opt$gene_annotation)
cat(sprintf("  loaded in %.1fs (cicero_cds: %d peaks, %d metacells)\n",
            as.numeric(Sys.time() - t0, units = "secs"),
            nrow(cicero_cds), ncol(cicero_cds)))

peaks_ch <- rownames(gene_ann)[gene_ann$chr == opt$chrom]
cat(sprintf("  %d peaks on %s\n", length(peaks_ch), opt$chrom))

if (length(peaks_ch) == 0) {
  cat(sprintf("  no peaks on %s — writing empty connections file\n", opt$chrom))
  out_path <- sprintf("conns_%s.tsv.gz", opt$chrom)
  gz <- gzfile(out_path, "w")
  writeLines("Peak1\tPeak2\tcoaccess", gz)
  close(gz)
  cat(sprintf("  wrote empty %s\n", out_path))
  q("no", 0)
}

cicero_cds_ch     <- cicero_cds[peaks_ch, ]
genomic_coords_ch <- data.frame(
  chr = gene_ann[peaks_ch, "chr"],
  bp1 = gene_ann[peaks_ch, "bp1"],
  bp2 = gene_ann[peaks_ch, "bp2"],
  row.names = peaks_ch,
  stringsAsFactors = FALSE
)

cat("\nGenerating cicero models + assembling connections...\n")
t1 <- Sys.time()
conns_ch <- tryCatch({
  models_ch <- generate_cicero_models(
    cicero_cds_ch,
    distance_parameter = opt$dp,
    window = opt$window_bp,
    s = opt$s,
    max_elements = opt$max_elements,
    genomic_coords = genomic_coords_ch
  )
  assemble_connections(models_ch, silent = TRUE)
}, error = function(e) {
  message(sprintf("  chrom %s failed: %s", opt$chrom, e$message))
  NULL
})

wall <- as.numeric(Sys.time() - t1, units = "secs")

if (is.null(conns_ch) || nrow(conns_ch) == 0) {
  cat(sprintf("  no connections for %s (wall=%.1fs); writing empty file\n", opt$chrom, wall))
  out_path <- sprintf("conns_%s.tsv.gz", opt$chrom)
  gz <- gzfile(out_path, "w")
  writeLines("Peak1\tPeak2\tcoaccess", gz)
  close(gz)
} else {
  cat(sprintf("  %d connections generated (wall=%.1fs)\n", nrow(conns_ch), wall))
  out_path <- sprintf("conns_%s.tsv.gz", opt$chrom)
  gz <- gzfile(out_path, "w")
  write.table(conns_ch, gz, sep = "\t", row.names = FALSE, quote = FALSE)
  close(gz)
  cat(sprintf("  wrote %s\n", out_path))
}
cat("== step 2 complete ==\n")
