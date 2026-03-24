#!/usr/bin/env Rscript
# MAST Differential Expression Analysis
# Adapted for Nextflow integration

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(tidyverse)
  library(optparse)
})

# Parse arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Input h5ad file (converted to RDS)"),
  make_option(c("--output_dir"), type="character", help="Output directory"),
  make_option(c("--cell_type"), type="character", help="Cell type to analyze"),
  make_option(c("--group1"), type="character", help="Group 1 (ident.1)"),
  make_option(c("--group2"), type="character", help="Group 2 (ident.2)"),
  make_option(c("--condition_key"), type="character", default="condition_group", 
             help="Column name for condition groups"),
  make_option(c("--cell_type_key"), type="character", default="scanvi_prediction",
             help="Column name for cell types"),
  make_option(c("--test_use"), type="character", default="MAST", help="Test method")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (any(sapply(opt[c("input", "output_dir", "cell_type", "group1", "group2")], is.null))) {
  print_help(opt_parser)
  stop("Missing required arguments")
}

# Create output directory
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== MAST Differential Expression ===\n")
cat("Input:", opt$input, "\n")
cat("Cell type:", opt$cell_type, "\n")
cat("Comparison:", opt$group1, "vs", opt$group2, "\n\n")

# Load Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(opt$input)

# Set cell type identity
cat("Setting cell type identities...\n")
Idents(seurat_obj) <- seurat_obj@meta.data[[opt$cell_type_key]]

# Subset for specific cell type
cat("Subsetting for", opt$cell_type, "...\n")
subset_obj <- subset(seurat_obj, idents = opt$cell_type)
cat("  Subset contains", ncol(subset_obj), "cells\n")

# Set condition identities
cat("Setting condition identities...\n")
Idents(subset_obj) <- subset_obj@meta.data[[opt$condition_key]]

# Validate groups exist
available_groups <- levels(Idents(subset_obj))
if (!(opt$group1 %in% available_groups)) {
  stop(paste("Group1", opt$group1, "not found. Available:", paste(available_groups, collapse=", ")))
}
if (!(opt$group2 %in% available_groups)) {
  stop(paste("Group2", opt$group2, "not found. Available:", paste(available_groups, collapse=", ")))
}

# Print cell counts per group
cat("\nCell counts per group:\n")
print(table(Idents(subset_obj)))

# Run MAST
cat("\nRunning FindMarkers with", opt$test_use, "...\n")
de_results <- tryCatch({
  FindMarkers(
    object = subset_obj,
    ident.1 = opt$group1,
    ident.2 = opt$group2,
    test.use = opt$test_use,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR in FindMarkers:", e$message, "\n")
  return(NULL)
})

if (is.null(de_results)) {
  stop("FindMarkers failed - check error messages above")
}

# Add gene names as column
de_results$gene <- rownames(de_results)

# Save results
output_file <- file.path(
  opt$output_dir,
  paste0(opt$cell_type, "_", opt$group1, "_vs_", opt$group2, "_DEResults.csv")
)

cat("\nSaving results to:", output_file, "\n")
write.csv(de_results, output_file, row.names = FALSE)

# Print summary
cat("\n=== Results Summary ===\n")
cat("Total genes tested:", nrow(de_results), "\n")
sig_genes <- de_results %>% filter(p_val_adj < 0.05)
cat("Significant genes (FDR < 0.05):", nrow(sig_genes), "\n")
up_genes <- sig_genes %>% filter(avg_log2FC > 0)
down_genes <- sig_genes %>% filter(avg_log2FC < 0)
cat("  Upregulated:", nrow(up_genes), "\n")
cat("  Downregulated:", nrow(down_genes), "\n")

cat("\n Analysis complete!\n")
