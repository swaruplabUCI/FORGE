#!/usr/bin/env Rscript
# Cicero Phase 3 — step 1 of 3: build shared CDS + estimate global distance_parameter
#
# Replaces the single-threaded run_cicero() bottleneck in the old CICERO_FULL.
# This script runs ONCE per pipeline run and persists its outputs to .rds for
# the per-chromosome fan-out downstream.
#
# Outputs:
#   cicero_cds_shared.rds   shared CiceroCellDataSet used by per-chrom tasks
#   gene_annotation.rds     peak -> (chr, bp1, bp2) table
#   distance_parameter.txt  scalar global dp (Phase 3 validated rho=1.0 vs C1)
#   input_cds_ordered.rds   monocle3 CDS with pseudotime (preserves downstream contract)
#   umap_clusters.pdf, umap_pseudotime.pdf, cds_summary.txt

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(monocle3)
  library(cicero)
})

option_list <- list(
  make_option(c("--triplets"), type = "character", help = "Path to cicero_triplets.tsv.gz"),
  make_option(c("--outdir"), type = "character", default = ".", help = "Output directory"),
  make_option(c("--num_dim"), type = "integer", default = 30, help = "Number of PCA dimensions"),
  make_option(c("--sample_num"), type = "integer", default = 100, help = "k for Cicero nearest-neighbor graph"),
  make_option(c("--window_bp"), type = "integer", default = 500000, help = "Window (bp) for distance_parameter estimation"),
  make_option(c("--distance_constraint"), type = "integer", default = 250000, help = "Max distance between coaccessible peaks"),
  make_option(c("--dp_convergence"), type = "double", default = 1e-10, help = "DP convergence threshold"),
  make_option(c("--use_partition"), action = "store_true", default = FALSE, help = "Use monocle3 partitions in learn_graph")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$triplets)) stop("Required: --triplets", call. = FALSE)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$outdir)

parse_peak <- function(p) {
  if (!(grepl(":", p, fixed = TRUE) && grepl("-", p, fixed = TRUE))) {
    return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  }
  parts <- strsplit(p, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  chr <- parts[1]
  coords <- strsplit(parts[2], "-", fixed = TRUE)[[1]]
  if (length(coords) < 2) return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  c(chr = chr,
    bp1 = suppressWarnings(as.numeric(coords[1])),
    bp2 = suppressWarnings(as.numeric(coords[2])))
}

cat("== Cicero Phase 3 step 1: global distance_parameter ==\n")
cat("triplets:           ", opt$triplets, "\n")
cat("num_dim:            ", opt$num_dim, "\n")
cat("sample_num:         ", opt$sample_num, "\n")
cat("window_bp:          ", opt$window_bp, "\n")
cat("distance_constraint:", opt$distance_constraint, "\n")

# -------- Load triplets; filter invalid peaks and chrUn_* scaffolds --------
cat("\nLoading triplets...\n")
trip <- read.delim(opt$triplets, header = TRUE, stringsAsFactors = FALSE)
required <- c("cell", "peak", "count")
missing_cols <- setdiff(required, colnames(trip))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
cat(sprintf("  loaded %d rows (%d cells, %d peaks)\n",
            nrow(trip), length(unique(trip$cell)), length(unique(trip$peak))))

# Edge-case filter: drop unplaced contigs chrUn_*
before_n <- nrow(trip)
trip <- trip[!grepl("^chrUn_", trip$peak), , drop = FALSE]
if (nrow(trip) < before_n) {
  cat(sprintf("  dropped %d rows on chrUn_* scaffolds\n", before_n - nrow(trip)))
}

all_peaks <- sort(unique(trip$peak))
peak_parsed <- lapply(all_peaks, parse_peak)
peak_df <- do.call(rbind, lapply(peak_parsed, function(x) data.frame(
  chr = x["chr"], bp1 = as.numeric(x["bp1"]), bp2 = as.numeric(x["bp2"]),
  stringsAsFactors = FALSE
)))
rownames(peak_df) <- all_peaks
valid_mask <- !is.na(peak_df$chr) & is.finite(peak_df$bp1) & is.finite(peak_df$bp2)
if (!all(valid_mask)) {
  cat(sprintf("  dropping %d peaks with invalid coords\n", sum(!valid_mask)))
  trip <- trip[trip$peak %in% all_peaks[valid_mask], , drop = FALSE]
}

# -------- Build sparse matrix + CDS --------
cells <- sort(unique(trip$cell))
peaks <- sort(unique(trip$peak))
mat <- sparseMatrix(
  i = match(trip$peak, peaks),
  j = match(trip$cell, cells),
  x = trip$count,
  dims = c(length(peaks), length(cells))
)
rownames(mat) <- peaks
colnames(mat) <- cells

pp2 <- lapply(peaks, parse_peak)
gene_ann <- data.frame(
  gene_short_name = peaks,
  chr = vapply(pp2, `[`, character(1), "chr"),
  bp1 = as.integer(as.numeric(vapply(pp2, `[`, character(1), "bp1"))),
  bp2 = as.integer(as.numeric(vapply(pp2, `[`, character(1), "bp2"))),
  row.names = peaks,
  stringsAsFactors = FALSE
)
if (!all(startsWith(gene_ann$chr, "chr"), na.rm = TRUE)) {
  gene_ann$chr <- paste0("chr", gene_ann$chr)
}

cell_meta <- data.frame(cell_id = cells, row.names = cells, stringsAsFactors = FALSE)

cat("\nCreating CellDataSet...\n")
cds <- new_cell_data_set(
  expression_data = mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_ann
)

cat(sprintf("Preprocessing (num_dim = %d)...\n", opt$num_dim))
cds <- preprocess_cds(cds, num_dim = opt$num_dim)
cat("Reducing dimensions (UMAP)...\n")
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cat("Clustering cells...\n")
cds <- cluster_cells(cds)
if (opt$use_partition && "partition" %in% colnames(colData(cds))) {
  cat("Learning graph using partitions...\n")
  cds <- learn_graph(cds, use_partition = TRUE)
} else {
  cat("Learning graph without partitions...\n")
  cds <- learn_graph(cds, use_partition = FALSE)
}
cat("Ordering cells (pseudotime)...\n")
root_cell <- colnames(cds)[1]
cds <- order_cells(cds, root_cells = root_cell)

saveRDS(cds, file = "input_cds_ordered.rds")
cat("  saved input_cds_ordered.rds\n")
writeLines(c(
  sprintf("Cells: %d", ncol(cds)),
  sprintf("Peaks: %d", nrow(cds)),
  sprintf("num_dim: %d", opt$num_dim)
), con = "cds_summary.txt")

# UMAP plots (preserve downstream artifacts from old CICERO_FULL)
pdf("umap_clusters.pdf", width = 6, height = 5)
print(plot_cells(cds, color_cells_by = "cluster",
                 label_groups_by_cluster = TRUE,
                 show_trajectory_graph = TRUE))
dev.off()
pdf("umap_pseudotime.pdf", width = 6, height = 5)
print(plot_cells(cds, color_cells_by = "pseudotime",
                 label_groups_by_cluster = FALSE,
                 show_trajectory_graph = TRUE))
dev.off()

# -------- Cicero CDS + global distance_parameter --------
cat("\nBuilding cicero CDS...\n")
umap_coords <- reducedDims(cds)$UMAP
if (is.null(umap_coords)) stop("Missing UMAP coords after reduce_dimension")
cicero_cds <- make_cicero_cds(cds, reduced_coordinates = umap_coords)

genomic_coords <- data.frame(
  chr = gene_ann$chr, bp1 = gene_ann$bp1, bp2 = gene_ann$bp2,
  row.names = rownames(gene_ann), stringsAsFactors = FALSE
)

cat("\nEstimating global distance_parameter...\n")
t0 <- Sys.time()
dp_vec <- estimate_distance_parameter(
  cicero_cds,
  window = opt$window_bp,
  maxit = 100,
  sample_num = opt$sample_num,
  distance_constraint = opt$distance_constraint,
  distance_parameter_convergence = opt$dp_convergence,
  genomic_coords = genomic_coords
)
dp <- mean(unlist(dp_vec))
wall <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("  dp = %.10f (from %d sampled windows, %.1fs)\n",
            dp, length(unlist(dp_vec)), wall))

# -------- Persist shared artifacts for per-chrom fan-out --------
saveRDS(cicero_cds, "cicero_cds_shared.rds")
saveRDS(gene_ann,   "gene_annotation.rds")
writeLines(sprintf("%.15g", dp), con = "distance_parameter.txt")
cat("\nSaved cicero_cds_shared.rds, gene_annotation.rds, distance_parameter.txt\n")
cat("== step 1 complete ==\n")
