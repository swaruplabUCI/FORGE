#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(monocle3)
  library(cicero)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("--triplets"), type = "character", help = "Path to cicero_triplets.tsv.gz"),
  make_option(c("--outdir"), type = "character", help = "Output directory"),
  make_option(c("--gtf"), type = "character", help = "Gene annotation GTF/GFF"),
  make_option(c("--num_dim"), type = "integer", default = 30, help = "Number of PCA dimensions [default %default]"),
  make_option(c("--sample_num"), type = "integer", default = 8, help = "k for Cicero nearest-neighbor graph [default %default]"),
  make_option(c("--connections_cutoff"), type = "double", default = 0.15, help = "Coaccessibility cutoff for retained connections [default %default]"),
  make_option(c("--ccan_min_coaccess"), type = "double", default = 0.15, help = "Minimum coaccessibility for CCANs [default %default]"),
  make_option(c("--use_partition"), action = "store_true", default = FALSE,
              help = "Use precomputed partitions in Monocle3 (if available) [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$triplets) || is.null(opt$outdir) || is.null(opt$gtf)) {
  stop("Required arguments: --triplets, --outdir, --gtf", call. = FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$outdir)

## ---------------------------------------------------------------------------
## Helper: parse peaks of the form "chrX:start-end"
## ---------------------------------------------------------------------------

# Parse peak IDs of the form "chrX:start-end"
parse_peak <- function(p) {
  # Require colon + dash
  if (!(grepl(":", p, fixed = TRUE) && grepl("-", p, fixed = TRUE))) {
    return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  }

  # Split "chr5:149184489-149188489"
  parts <- strsplit(p, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  }

  chr <- parts[1]
  coords <- strsplit(parts[2], "-", fixed = TRUE)[[1]]
  if (length(coords) < 2) {
    return(c(chr = NA_character_, bp1 = NA_real_, bp2 = NA_real_))
  }

  bp1 <- suppressWarnings(as.numeric(coords[1]))
  bp2 <- suppressWarnings(as.numeric(coords[2]))
  c(chr = chr, bp1 = bp1, bp2 = bp2)
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Cicero full run: Conns / CCAN / Pseudotime / Plots\n")
cat(rep("=", 60), "\n\n")

## ---------------------------------------------------------------------------
## 1. Load triplets and filter to peaks with valid genomic coordinates
## ---------------------------------------------------------------------------

cat("Loading triplets from: ", opt$triplets, "\n", sep = "")
trip <- read.delim(opt$triplets, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("cell", "peak", "count")
missing_cols <- setdiff(required_cols, colnames(trip))
if (length(missing_cols) > 0) {
  stop("Triplets file must contain columns: cell, peak, count. Missing: ",
       paste(missing_cols, collapse = ", "))
}

cat(sprintf("  ✓ Loaded %d triplets (%d unique cells, %d unique peaks)\n",
            nrow(trip),
            length(unique(trip$cell)),
            length(unique(trip$peak))))

# Parse all unique peaks once
all_peaks <- sort(unique(trip$peak))

peak_parsed_list <- lapply(all_peaks, parse_peak)
peak_df <- do.call(
  rbind,
  lapply(peak_parsed_list, function(x) {
    data.frame(
      chr = x["chr"],
      bp1 = as.numeric(x["bp1"]),
      bp2 = as.numeric(x["bp2"]),
      stringsAsFactors = FALSE
    )
  })
)
rownames(peak_df) <- all_peaks

valid_peaks <- !is.na(peak_df$chr) &
               is.finite(peak_df$bp1) &
               is.finite(peak_df$bp2)

if (!all(valid_peaks)) {
  n_drop <- sum(!valid_peaks)
  cat(sprintf("  ! Dropping %d peaks with invalid coordinates out of %d total peaks\n",
              n_drop, length(all_peaks)))
}

valid_peak_names <- all_peaks[valid_peaks]

# Filter triplets to only valid peaks
trip <- trip[trip$peak %in% valid_peak_names, , drop = FALSE]

cat(sprintf("  ✓ After filtering: %d triplets (%d unique cells, %d unique peaks)\n",
            nrow(trip),
            length(unique(trip$cell)),
            length(unique(trip$peak))))

## ---------------------------------------------------------------------------
## 2. Build sparse cell x peak matrix from filtered triplets
## ---------------------------------------------------------------------------

cells <- sort(unique(trip$cell))
peaks <- sort(unique(trip$peak))

cell_index <- match(trip$cell, cells)
peak_index <- match(trip$peak, peaks)

mat <- sparseMatrix(
  i = peak_index,
  j = cell_index,
  x = trip$count,
  dims = c(length(peaks), length(cells))
)
rownames(mat) <- peaks
colnames(mat) <- cells

cat("  ✓ Constructed sparse matrix: ", nrow(mat), " peaks x ",
    ncol(mat), " cells\n", sep = "")

## ---------------------------------------------------------------------------
## 3. Build monocle3 CellDataSet (with valid peak coordinates)
## ---------------------------------------------------------------------------

# Parse coordinates for the filtered peaks
peak_parsed2 <- lapply(peaks, parse_peak)
peak_df2 <- do.call(
  rbind,
  lapply(peak_parsed2, function(x) {
    data.frame(
      chr = x["chr"],
      bp1 = as.numeric(x["bp1"]),
      bp2 = as.numeric(x["bp2"]),
      stringsAsFactors = FALSE
    )
  })
)
rownames(peak_df2) <- peaks

gene_annotation <- data.frame(
  gene_short_name = peaks,
  chr = as.character(peak_df2$chr),
  bp1 = as.integer(peak_df2$bp1),
  bp2 = as.integer(peak_df2$bp2),
  row.names = peaks,
  stringsAsFactors = FALSE
)

# Ensure all chromosomes have "chr" prefix (they already should)
if (!all(startsWith(gene_annotation$chr, "chr"), na.rm = TRUE)) {
  gene_annotation$chr <- paste0("chr", gene_annotation$chr)
}

cell_metadata <- data.frame(
  cell_id = cells,
  row.names = cells,
  stringsAsFactors = FALSE
)

cat("Creating CellDataSet...\n")
cds <- new_cell_data_set(
  expression_data = mat,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

## ---------------------------------------------------------------------------
## 4. Preprocess, reduce dimensions, cluster, learn graph, pseudotime
## ---------------------------------------------------------------------------

cat("Preprocessing CDS (PCA, num_dim = ", opt$num_dim, ")...\n", sep = "")
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
cat("  ✓ Saved ordered CDS -> input_cds_ordered.rds\n")

summary_lines <- c(
  sprintf("Cells: %d", ncol(cds)),
  sprintf("Peaks: %d", nrow(cds)),
  sprintf("UMAP dims: %d", opt$num_dim)
)
writeLines(summary_lines, con = "cds_summary.txt")

## ---------------------------------------------------------------------------
## 5. Build Cicero CDS and compute co-accessibility
## ---------------------------------------------------------------------------

cat("\nBuilding Cicero CDS and computing co-accessibility...\n")

genomic_coords <- data.frame(
  chr = gene_annotation$chr,
  bp1 = gene_annotation$bp1,
  bp2 = gene_annotation$bp2,
  row.names = rownames(gene_annotation),
  stringsAsFactors = FALSE
)

umap_coords <- reducedDims(cds)$UMAP
if (is.null(umap_coords)) {
  stop("UMAP coordinates not found in cds@reducedDims. Check reduce_dimension call.")
}

cicero_cds <- make_cicero_cds(
  cds,
  reduced_coordinates = umap_coords
)

conns <- run_cicero(
  cicero_cds,
  genomic_coords,
  sample_num = opt$sample_num
)

cat(sprintf("  ✓ Cicero connections computed: %d rows\n", nrow(conns)))

conns_filtered <- conns %>% dplyr::filter(coaccess >= opt$connections_cutoff)
cat(sprintf("  ✓ Retained %d connections with coaccess >= %.3f\n",
            nrow(conns_filtered), opt$connections_cutoff))

gz_conns <- gzfile("cicero_connections.tsv.gz", "w")
write.table(conns_filtered, gz_conns, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz_conns)
cat("  ✓ Saved connections -> cicero_connections.tsv.gz\n")

# ---------------------------------------------------------------------------
# 6. CCANs
# ---------------------------------------------------------------------------

cat("\nFinding CCANs...\n")
ccans <- generate_ccans(
  conns_filtered,
  coaccess_cutoff = opt$ccan_min_coaccess
)
cat(sprintf("  ✓ Identified %d CCANs\n", length(unique(ccans$CCAN))))

gz_ccan <- gzfile("CCAN_assignments.tsv.gz", "w")
write.table(ccans, gz_ccan, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz_ccan)
cat("  ✓ Saved CCAN assignments -> CCAN_assignments.tsv.gz\n")

# CCAN size plots (top bar chart + distribution histogram)
if (!is.null(ccans) && NROW(ccans) > 0 && "CCAN" %in% colnames(ccans)) {
  ccan_df <- as.data.frame(ccans)
  ccan_sizes <- ccan_df %>% dplyr::count(CCAN) %>% dplyr::arrange(dplyr::desc(n))
  n_ccans <- nrow(ccan_sizes)
  top_n <- min(50, n_ccans)
  ccan_top <- ccan_sizes[1:top_n, ]
  ccan_top$CCAN <- factor(ccan_top$CCAN, levels = ccan_top$CCAN)

  pdf("ccan_top_graph.pdf", width = max(6, top_n * 0.15), height = 4)
  print(
    ggplot(ccan_top, aes(x = CCAN, y = n)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
      xlab(sprintf("CCAN ID (top %d of %d)", top_n, n_ccans)) +
      ylab("Number of peaks") +
      ggtitle(sprintf("CCAN Sizes (top %d of %d)", top_n, n_ccans))
  )
  dev.off()
  cat("  ✓ Saved CCAN top bar chart -> ccan_top_graph.pdf\n")

  pdf("ccan_size_distribution.pdf", width = 6, height = 4)
  print(
    ggplot(ccan_sizes, aes(x = n)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      xlab("Number of peaks per CCAN") + ylab("Count") +
      ggtitle(sprintf("CCAN Size Distribution (n = %d CCANs)", n_ccans))
  )
  dev.off()
  cat("  ✓ Saved CCAN size distribution -> ccan_size_distribution.pdf\n")
} else {
  cat("  ! CCAN object has no rows or no 'CCAN' column; skipping CCAN size plot\n")
}

## ---------------------------------------------------------------------------
## 7. UMAP & pseudotime plots
## ---------------------------------------------------------------------------

cat("\nGenerating UMAP and pseudotime plots...\n")

pdf("umap_clusters.pdf", width = 6, height = 5)
plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = TRUE,
  show_trajectory_graph = TRUE
)
dev.off()
cat("  ✓ Saved UMAP clusters plot -> umap_clusters.pdf\n")

pdf("umap_pseudotime.pdf", width = 6, height = 5)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  show_trajectory_graph = TRUE
)
dev.off()
cat("  ✓ Saved UMAP pseudotime plot -> umap_pseudotime.pdf\n")

pd <- as.data.frame(colData(cds))
if (!is.null(pd$pseudotime)) {
  pdf("accessibility_pseudotime.pdf", width = 6, height = 4)
  ggplot(pd, aes(x = pseudotime)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    theme_minimal() +
    xlab("Pseudotime") +
    ylab("Density") +
    ggtitle("Cell density along pseudotime")
  dev.off()
  cat("  ✓ Saved accessibility-pseudotime plot -> accessibility_pseudotime.pdf\n")
}

# Accessibility changes across pseudotime (Cicero plot_accessibility_in_pseudotime)
tryCatch({
  if (!is.null(pd$pseudotime) && sum(is.finite(pd$pseudotime)) > 50) {
    cat("Running graph_test for pseudotime-variable peaks...\n")
    pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
    sig_peaks <- pr_test_res[!is.na(pr_test_res$q_value) & pr_test_res$q_value < 0.05, ]
    sig_peaks <- sig_peaks[order(sig_peaks$q_value), ]

    if (nrow(sig_peaks) >= 3) {
      top_n <- min(10, nrow(sig_peaks))
      top_peak_ids <- rownames(sig_peaks)[1:top_n]
      cds_finite <- cds[, is.finite(pseudotime(cds))]
      cds_subset <- cds_finite[top_peak_ids, ]

      pdf("accessibility_in_pseudotime.pdf", width = 10, height = 2 * top_n)
      plot_accessibility_in_pseudotime(cds_subset, breaks = 10)
      dev.off()
      cat(sprintf("  ✓ Saved accessibility across pseudotime (top %d peaks) -> accessibility_in_pseudotime.pdf\n", top_n))
    } else {
      cat("  ! Fewer than 3 significant peaks; skipping accessibility-in-pseudotime plot\n")
    }
  }
}, error = function(e) {
  cat(sprintf("  ! Accessibility-in-pseudotime plot failed: %s\n", e$message))
})

## ---------------------------------------------------------------------------
## 8. Simple regional connections plot (using colon-dash peaks)
## ---------------------------------------------------------------------------

cat("\nGenerating simple regional connections plot...\n")

if (nrow(conns_filtered) > 0) {
  peak_counts <- c(
    table(conns_filtered$Peak1),
    table(conns_filtered$Peak2)
  )
  peak_counts <- sort(peak_counts, decreasing = TRUE)
  top_peak <- names(peak_counts)[1]

  pcs <- parse_peak(top_peak)
  chr   <- pcs["chr"]
  start <- as.numeric(pcs["bp1"])
  end   <- as.numeric(pcs["bp2"])

  if (!is.na(chr) && is.finite(start) && is.finite(end)) {
    window <- 250000
    minbp <- max(0, start - window)
    maxbp <- end + window

    cat("Loading gene annotation from: ", opt$gtf, "\n", sep = "")
    gene_anno <- rtracklayer::import(opt$gtf)
    gene_anno_df <- as.data.frame(gene_anno)
    gene_anno_df$chromosome <- as.character(seqnames(gene_anno))
    if (!all(startsWith(gene_anno_df$chromosome, "chr"))) {
      gene_anno_df$chromosome <- paste0("chr", gene_anno_df$chromosome)
    }
    gene_anno_df$gene <- if ("gene_id" %in% names(gene_anno_df)) {
      as.character(gene_anno_df$gene_id)
    } else NA_character_
    gene_anno_df$transcript <- if ("transcript_id" %in% names(gene_anno_df)) {
      as.character(gene_anno_df$transcript_id)
    } else NA_character_
    gene_anno_df$symbol <- if ("gene_name" %in% names(gene_anno_df)) {
      as.character(gene_anno_df$gene_name)
    } else NA_character_

    pdf("connections_region.pdf", width = 10, height = 5)
    try({
      plot_connections(
        conns_filtered,
        chr = chr,
        minbp = minbp,
        maxbp = maxbp,
        gene_model = gene_anno_df,
        coaccess_cutoff = opt$connections_cutoff,
        connection_width = 0.5,
        collapseTranscripts = "longest",
        viewpoint = top_peak
      )
    }, silent = TRUE)
    dev.off()
    cat("  ✓ Saved regional connections plot -> connections_region.pdf\n")
  } else {
    cat("  ! Could not parse top peak name for regional plot; skipping\n")
  }
} else {
  cat("  ! No filtered connections found; skipping regional plot\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Cicero full run complete. Outputs written to: ", opt$outdir, "\n", sep = "")
cat(rep("=", 60), "\n\n")
