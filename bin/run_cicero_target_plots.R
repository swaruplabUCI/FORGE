#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(cicero)
  library(monocle3)
  library(rtracklayer)
  library(igraph)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

option_list <- list(
  make_option(c("--connections"), type = "character", help = "cicero_connections.tsv.gz"),
  make_option(c("--ccan"),        type = "character", help = "CCAN_assignments.tsv.gz"),
  make_option(c("--cds"),         type = "character", help = "input_cds_ordered.rds"),
  make_option(c("--gtf"),         type = "character", help = "GTF/GFF3 for gene annotations"),
  make_option(c("--genes"),       type = "character", help = "Comma-separated gene symbols, e.g. 'Hmgb1,Hspa8'"),
  make_option(c("--outdir"),      type = "character", help = "Output directory for gene-specific plots"),
  make_option(c("--windows"),     type = "character", default = "100000,250000",
              help = "Comma-separated window sizes in bp [default %default]"),
  make_option(c("--thresholds"),  type = "character", default = "0.10,0.05",
              help = "Comma-separated coaccessibility cutoffs [default %default]"),
  make_option(c("--labels"),      type = "character", default = "strict,permissive",
              help = "Comma-separated labels for each window/threshold pair [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Parse comma-separated arguments
windows    <- as.numeric(strsplit(opt$windows, ",")[[1]])
thresholds <- as.numeric(strsplit(opt$thresholds, ",")[[1]])
labels     <- strsplit(opt$labels, ",")[[1]]

if (length(windows) != length(thresholds) || length(windows) != length(labels)) {
  stop("Error: --windows, --thresholds, and --labels must have the same number of elements.", call. = FALSE)
}

plot_settings <- Map(function(w, t, l) list(window = w, threshold = t, label = l), windows, thresholds, labels)

if (is.null(opt$connections) || is.null(opt$gtf) || is.null(opt$genes)) {
  stop("Missing required arguments: --connections, --gtf, --genes", call. = FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("Loading connections...\n")
conns <- read.delim(opt$connections, stringsAsFactors = FALSE)

cat("Loading CDS...\n")
if (!is.null(opt$cds) && file.exists(opt$cds)) {
    input_cds <- readRDS(opt$cds)
} else {
    input_cds <- NULL
    cat("Warning: CDS file not provided or not found. Viewpoint plots might lack peak info if needed.\n")
}

cat("Loading Gene Annotations (GTF)...\n")
# plot_connections requires a specific gene_model format:
# chromosome, start, end, gene, transcript, strand, symbol
gene_anno <- rtracklayer::import(opt$gtf)
gene_anno_df <- as.data.frame(gene_anno)

# Standardize column names for Cicero plot_connections
# It expects: chromosome, start, end, gene, transcript, strand, symbol
# 'gene' is gene_id, 'symbol' is gene_name usually.

# Ensure chromosome starts with 'chr'
if (!all(startsWith(as.character(gene_anno_df$seqnames), "chr"))) {
    gene_anno_df$seqnames <- paste0("chr", as.character(gene_anno_df$seqnames))
}

gene_model <- gene_anno_df %>%
  dplyr::select(seqnames, start, end, strand, type, gene_id, gene_name, transcript_id) %>%
  dplyr::rename(
    chromosome = seqnames,
    gene = gene_id,
    symbol = gene_name,
    transcript = transcript_id
  ) %>%
  dplyr::filter(type %in% c("exon", "gene", "transcript"))

# Ensure symbol is not NA
gene_model$symbol <- ifelse(is.na(gene_model$symbol), gene_model$gene, gene_model$symbol)
gene_model$transcript <- ifelse(is.na(gene_model$transcript), gene_model$gene, gene_model$transcript)

target_genes <- unique(trimws(strsplit(opt$genes, ",")[[1]]))

cat(sprintf("Processing %d target genes...\n", length(target_genes)))

for (gene_sym in target_genes) {
    cat(sprintf("Plotting %s...\n", gene_sym))
    
    # Find gene location
    gene_dat <- gene_model %>% filter(symbol == gene_sym)
    if (nrow(gene_dat) == 0) {
        cat(sprintf("  Gene %s not found in annotation.\n", gene_sym))
        next
    }
    
    # Determine plotting range based on gene extent
    gene_chr <- as.character(gene_dat$chromosome[1])
    gene_start <- min(gene_dat$start)
    gene_end <- max(gene_dat$end)
    gene_tss <- ifelse(as.character(gene_dat$strand[1]) == "+", gene_start, gene_end)
    
    for (setting in plot_settings) {
        window_bp <- setting$window
        cutoff <- setting$threshold
        label <- setting$label
        
        plot_min <- max(1, gene_tss - window_bp)
        plot_max <- gene_tss + window_bp
        
        pdf_file <- file.path(opt$outdir, sprintf("%s_connections_%s.pdf", gene_sym, label))
        
        cat(sprintf("  %s: %s:%d-%d (cutoff %.2f)\n", label, gene_chr, plot_min, plot_max, cutoff))
        
        pdf(pdf_file, width = 10, height = 6)
        tryCatch({
            # plot_connections from Cicero
            plot_connections(
                conns,
                chr = gene_chr,
                minbp = plot_min,
                maxbp = plot_max,
                gene_model = gene_model,
                coaccess_cutoff = cutoff,
                connection_width = 0.5,
                collapseTranscripts = "longest",
                include_axis_track = TRUE,
                alpha_by_coaccess = TRUE
            )
        }, error = function(e) {
            cat(sprintf("  Error plotting %s: %s\n", gene_sym, e$message))
            # Create empty plot with error message
            plot(0, 0, type="n", axes=FALSE, xlab="", ylab="")
            text(0, 0, paste("Error plotting", gene_sym, ":\n", e$message), cex=0.8)
        })
        dev.off()
    }
}

cat("Done.\n")
