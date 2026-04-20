#!/usr/bin/env Rscript
# Cicero Phase 3 — step 3 of 3: join per-chromosome connections + generate CCANs + plots
#
# Reads all conns_*.tsv.gz files (from CICERO_FULL_CHROM fan-out), rbinds,
# filters by connections_cutoff, writes cicero_connections.tsv.gz matching
# the old CICERO_FULL output contract (columns: Peak1, Peak2, coaccess).
# Also generates CCANs and ancillary plots/TSVs downstream consumers expect.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(cicero)
  library(monocle3)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("--conns_glob"), type = "character", default = "conns_*.tsv.gz",
              help = "Glob for per-chrom connection TSVs"),
  make_option(c("--cds"), type = "character", default = "input_cds_ordered.rds",
              help = "Ordered monocle3 CDS from step 1 (for pseudotime plots)"),
  make_option(c("--gtf"), type = "character", default = NULL,
              help = "Gene annotation GTF (for regional plot)"),
  make_option(c("--connections_cutoff"), type = "double", default = 0.15),
  make_option(c("--ccan_min_coaccess"), type = "double", default = 0.15),
  make_option(c("--outdir"), type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$outdir)

cat("== Cicero Phase 3 step 3: join + CCAN + plots ==\n")

# -------- Load + rbind per-chrom connections --------
files <- Sys.glob(opt$conns_glob)
cat(sprintf("Loading %d per-chrom conn files...\n", length(files)))
if (length(files) == 0) {
  stop("No conns_*.tsv.gz files matched glob: ", opt$conns_glob, call. = FALSE)
}

conns_list <- lapply(files, function(f) {
  df <- tryCatch(
    read.delim(f, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) { message("  failed to read ", f, ": ", e$message); NULL }
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  required <- c("Peak1", "Peak2", "coaccess")
  if (!all(required %in% colnames(df))) {
    message("  ", f, " missing required columns; skipping")
    return(NULL)
  }
  df
})
conns_list <- conns_list[!vapply(conns_list, is.null, logical(1))]

if (length(conns_list) == 0) {
  stop("All per-chrom conn files were empty or malformed", call. = FALSE)
}

conns <- do.call(rbind, conns_list)
cat(sprintf("  combined %d connections from %d non-empty chrom files\n",
            nrow(conns), length(conns_list)))

# -------- Filter by connections_cutoff --------
conns_filt <- conns[!is.na(conns$coaccess) & conns$coaccess >= opt$connections_cutoff, ,
                    drop = FALSE]
cat(sprintf("  retained %d connections with coaccess >= %.3f\n",
            nrow(conns_filt), opt$connections_cutoff))

gz <- gzfile("cicero_connections.tsv.gz", "w")
write.table(conns_filt, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)
cat("  wrote cicero_connections.tsv.gz\n")

# -------- CCANs --------
cat("\nGenerating CCANs...\n")
ccans <- generate_ccans(conns_filt, coaccess_cutoff = opt$ccan_min_coaccess)
cat(sprintf("  identified %d CCANs\n", length(unique(ccans$CCAN))))

gz_ccan <- gzfile("CCAN_assignments.tsv.gz", "w")
write.table(ccans, gz_ccan, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz_ccan)
cat("  wrote CCAN_assignments.tsv.gz\n")

# CCAN plots
if (!is.null(ccans) && NROW(ccans) > 0 && "CCAN" %in% colnames(ccans)) {
  ccan_df <- as.data.frame(ccans)
  ccan_sizes <- ccan_df %>% dplyr::count(CCAN) %>% dplyr::arrange(dplyr::desc(n))
  n_ccans <- nrow(ccan_sizes)
  top_n <- min(50, n_ccans)
  ccan_top <- ccan_sizes[1:top_n, ]
  ccan_top$CCAN <- factor(ccan_top$CCAN, levels = ccan_top$CCAN)

  pdf("ccan_top_graph.pdf", width = max(6, top_n * 0.15), height = 4)
  print(ggplot(ccan_top, aes(x = CCAN, y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    xlab(sprintf("CCAN ID (top %d of %d)", top_n, n_ccans)) +
    ylab("Number of peaks") +
    ggtitle(sprintf("CCAN Sizes (top %d of %d)", top_n, n_ccans)))
  dev.off()

  pdf("ccan_size_distribution.pdf", width = 6, height = 4)
  print(ggplot(ccan_sizes, aes(x = n)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    xlab("Number of peaks per CCAN") + ylab("Count") +
    ggtitle(sprintf("CCAN Size Distribution (n = %d CCANs)", n_ccans)))
  dev.off()
}

# -------- Pseudotime-adjacent plots (load ordered CDS from step 1) --------
if (file.exists(opt$cds)) {
  cat(sprintf("\nLoading ordered CDS from %s...\n", opt$cds))
  cds <- tryCatch(readRDS(opt$cds), error = function(e) { message(e$message); NULL })

  if (!is.null(cds)) {
    pd <- as.data.frame(colData(cds))
    if (!is.null(pd$pseudotime)) {
      pdf("accessibility_pseudotime.pdf", width = 6, height = 4)
      print(ggplot(pd, aes(x = pseudotime)) +
        geom_density(fill = "steelblue", alpha = 0.5) +
        theme_minimal() +
        xlab("Pseudotime") + ylab("Density") +
        ggtitle("Cell density along pseudotime"))
      dev.off()

      tryCatch({
        if (sum(is.finite(pd$pseudotime)) > 50) {
          cat("Running graph_test for pseudotime-variable peaks...\n")
          pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
          sig_peaks <- pr_test_res[!is.na(pr_test_res$q_value) & pr_test_res$q_value < 0.05, ]
          sig_peaks <- sig_peaks[order(sig_peaks$q_value), ]
          if (nrow(sig_peaks) >= 3) {
            top_peaks <- rownames(sig_peaks)[1:min(10, nrow(sig_peaks))]
            cds_finite <- cds[, is.finite(pseudotime(cds))]
            cds_subset <- cds_finite[top_peaks, ]
            pdf("accessibility_in_pseudotime.pdf", width = 10, height = 2 * length(top_peaks))
            plot_accessibility_in_pseudotime(cds_subset, breaks = 10)
            dev.off()
          }
        }
      }, error = function(e) cat(sprintf("  accessibility-in-pseudotime plot failed: %s\n", e$message)))
    }

    # -------- Regional connections plot around the most-connected peak --------
    if (!is.null(opt$gtf) && file.exists(opt$gtf) && nrow(conns_filt) > 0) {
      peak_counts <- c(table(conns_filt$Peak1), table(conns_filt$Peak2))
      peak_counts <- sort(peak_counts, decreasing = TRUE)
      top_peak <- names(peak_counts)[1]
      parts <- strsplit(top_peak, ":", fixed = TRUE)[[1]]
      if (length(parts) == 2) {
        coords <- strsplit(parts[2], "-", fixed = TRUE)[[1]]
        if (length(coords) == 2) {
          chr <- parts[1]
          start <- suppressWarnings(as.numeric(coords[1]))
          end   <- suppressWarnings(as.numeric(coords[2]))
          if (!is.na(chr) && is.finite(start) && is.finite(end)) {
            window <- 250000
            minbp <- max(0, start - window)
            maxbp <- end + window
            cat(sprintf("Loading gene annotation from %s...\n", opt$gtf))
            gene_anno <- rtracklayer::import(opt$gtf)
            gene_anno_df <- as.data.frame(gene_anno)
            gene_anno_df$chromosome <- as.character(seqnames(gene_anno))
            if (!all(startsWith(gene_anno_df$chromosome, "chr"))) {
              gene_anno_df$chromosome <- paste0("chr", gene_anno_df$chromosome)
            }
            gene_anno_df$gene       <- if ("gene_id"       %in% names(gene_anno_df)) as.character(gene_anno_df$gene_id)       else NA_character_
            gene_anno_df$transcript <- if ("transcript_id" %in% names(gene_anno_df)) as.character(gene_anno_df$transcript_id) else NA_character_
            gene_anno_df$symbol     <- if ("gene_name"     %in% names(gene_anno_df)) as.character(gene_anno_df$gene_name)     else NA_character_

            pdf("connections_region.pdf", width = 10, height = 5)
            try(plot_connections(conns_filt, chr = chr, minbp = minbp, maxbp = maxbp,
                                 gene_model = gene_anno_df,
                                 coaccess_cutoff = opt$connections_cutoff,
                                 connection_width = 0.5,
                                 collapseTranscripts = "longest",
                                 viewpoint = top_peak), silent = TRUE)
            dev.off()
          }
        }
      }
    }
  }
} else {
  cat(sprintf("\nordered CDS not found at %s; skipping pseudotime plots\n", opt$cds))
}

cat("\n== step 3 complete ==\n")
