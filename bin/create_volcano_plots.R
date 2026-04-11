#!/usr/bin/env Rscript
# Create volcano plots from MAST DE results

suppressPackageStartupMessages({
  library(EnhancedVolcano)
  library(tidyverse)
  library(optparse)
})

# Parse arguments
option_list <- list(
  make_option(c("--input"), type="character", help="Input DE results CSV"),
  make_option(c("--output"), type="character", help="Output PDF path"),
  make_option(c("--species"), type="character", default="human", help="Species (human/mouse)"),
  make_option(c("--pval_cutoff"), type="numeric", default=0.05, help="P-value cutoff"),
  make_option(c("--fc_cutoff"), type="numeric", default=0.25, help="Log2FC cutoff")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("=== Volcano Plot Generation ===\n")
cat("Input:", opt$input, "\n")
cat("Output:", opt$output, "\n")
cat("Species:", opt$species, "\n\n")

# Load DE results
cat("Loading DE results...\n")
de_results <- read_csv(opt$input, show_col_types = FALSE)

# Handle zero p-values
non_zero_pvals <- de_results$p_val_adj[de_results$p_val_adj != 0]
if (length(non_zero_pvals) > 0) {
  min_pval <- min(non_zero_pvals)
  de_results$p_val_adj[de_results$p_val_adj == 0] <- min_pval
  cat("Replaced", sum(de_results$p_val_adj == min_pval), "zero p-values with minimum:", min_pval, "\n")
}

# Optionally filter to protein-coding genes via biomaRt (requires internet)
de_plot <- de_results
tryCatch({
  suppressPackageStartupMessages(library(biomaRt))
  if (opt$species == "human") {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  } else {
    ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
  }
  all_gene_biotypes <- getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
    mart = ensembl
  )
  protein_coding_genes <- all_gene_biotypes %>% filter(gene_biotype == "protein_coding")
  de_plot <- de_results %>% filter(gene %in% protein_coding_genes$external_gene_name)
  cat("Filtered to", nrow(de_plot), "protein-coding genes (via biomaRt)\n\n")
}, error = function(e) {
  cat("WARNING: biomaRt unavailable (offline node) — using all", nrow(de_results), "genes\n\n")
  de_plot <<- de_results
})

# Extract comparison info from filename
filename <- basename(opt$input)
comparison <- gsub("_DEResults.csv", "", filename)

# Create volcano plot
cat("Creating volcano plot...\n")
pdf(opt$output, width = 10, height = 8)

print(
  EnhancedVolcano(
    de_plot,
    lab = de_plot$gene,
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = comparison,
    pCutoff = opt$pval_cutoff,
    FCcutoff = opt$fc_cutoff,
    pointSize = 2.0,
    labSize = 4.0,
    colAlpha = 0.5,
    legendPosition = 'right',
    drawConnectors = TRUE,
    widthConnectors = 0.5
  )
)

dev.off()

cat(" Volcano plot saved to", opt$output, "\n")
