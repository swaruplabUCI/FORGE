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
# window sizes (bp) and thresholds, comma-separated
  make_option(c("--windows"),     type = "character", default = "100000,250000,500000",
              help = "Comma-separated window sizes in bp, e.g. '100000,250000,500000' [default %default]"),
  make_option(c("--thresholds"),  type = "character", default = "0.10,0.05,0.025",
              help = "Comma-separated coaccessibility cutoffs, e.g. '0.10,0.05,0.025' [default %default]"),
  make_option(c("--labels"),      type = "character", default = "strict,permissive,exploratory",
              help = "Comma-separated labels for each window/threshold pair [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

windows    <- as.numeric(strsplit(opt$windows, ",")[[1]])
thresholds <- as.numeric(strsplit(opt$thresholds, ",")[[1]])
labels     <- strsplit(opt$labels, ",")[[1]]

if (!(length(windows) == length(thresholds) && length(windows) == length(labels))) {
  stop("windows, thresholds, and labels must all have the same length.", call. = FALSE)
}

plot_settings <- Map(function(w, t, lab) {
  list(window = w, threshold = t, suffix = lab)
}, windows, thresholds, labels)

if (is.null(opt$connections) || is.null(opt$ccan) || is.null(opt$cds) ||
    is.null(opt$gtf) || is.null(opt$genes) || is.null(opt$outdir)) {
  stop("Required: --connections, --ccan, --cds, --gtf, --genes, --outdir", call. = FALSE)
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 60), "\n", sep = "")
cat("CICERO GENE-SPECIFIC CONNECTION PLOTS\n")
cat(rep("=", 60), "\n\n")

# ---------------------------------------------------------------------------
# STEP 1: Load Cicero objects
# ---------------------------------------------------------------------------

cat("Loading connections from: ", opt$connections, "\n", sep = "")
conns <- read.delim(opt$connections, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Loaded %d connections\n", nrow(conns)))

cat("Loading CCAN assignments from: ", opt$ccan, "\n", sep = "")
ccan_assigns <- read.delim(opt$ccan, stringsAsFactors = FALSE)
cat(sprintf("  ✓ Loaded %d peak assignments\n", nrow(ccan_assigns)))

cat("Loading ordered CDS from: ", opt$cds, "\n", sep = "")
input_cds <- readRDS(opt$cds)
cat(sprintf("  ✓ Loaded CDS: %d cells x %d peaks\n", ncol(input_cds), nrow(input_cds)))

# ---------------------------------------------------------------------------
# STEP 2: Load gene annotations & compute TSS per gene symbol
# ---------------------------------------------------------------------------

cat("\nLoading gene annotations from: ", opt$gtf, "\n", sep = "")
gene_anno <- rtracklayer::import(opt$gtf)
gene_anno_df <- as.data.frame(gene_anno)

# Standardize chromosome naming
gene_anno_df$chromosome <- as.character(seqnames(gene_anno))
if (!all(startsWith(gene_anno_df$chromosome, "chr"))) {
  gene_anno_df$chromosome <- paste0("chr", gene_anno_df$chromosome)
}

# Build a plot_connections()-compatible gene model.
# In your current cicero install, generate_gene_model is not available/exported,
# so we construct a minimal data.frame: chr, start, end, gene
cat("Building gene_model from imported annotations (no cicero::generate_gene_model available)...\n")

# Pick a gene label column robustly across GTF/GFF3 imports
gene_label <- rep(NA_character_, nrow(gene_anno_df))
for (col in c("gene_name", "Name", "gene", "gene_id")) {
  if (col %in% names(gene_anno_df)) {
    gene_label <- dplyr::coalesce(gene_label, as.character(gene_anno_df[[col]]))
  }
}

# Prefer 'gene' features; fallback to 'transcript'; else use all rows
if ("type" %in% names(gene_anno_df)) {
  gene_rows <- gene_anno_df %>% dplyr::filter(type == "gene")
  if (nrow(gene_rows) == 0) {
    gene_rows <- gene_anno_df %>% dplyr::filter(type == "transcript")
  }
  if (nrow(gene_rows) == 0) {
    gene_rows <- gene_anno_df
  }
} else {
  gene_rows <- gene_anno_df
}

gene_model <- gene_rows %>%
  dplyr::mutate(gene = gene_label[match(row.names(.), row.names(gene_anno_df))]) %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::transmute(
    chr   = chromosome,
    start = as.numeric(start),
    end   = as.numeric(end),
    gene  = as.character(gene)
  )

cat(sprintf("  ✓ gene_model rows: %d\n", nrow(gene_model)))

# Basic columns for TSS lookup
gene_anno_df$gene_id       <- if ("gene_id"       %in% names(gene_anno_df)) as.character(gene_anno_df$gene_id)       else NA_character_
gene_anno_df$gene_name     <- if ("gene_name"     %in% names(gene_anno_df)) as.character(gene_anno_df$gene_name)     else NA_character_
gene_anno_df$transcript_id <- if ("transcript_id" %in% names(gene_anno_df)) as.character(gene_anno_df$transcript_id) else NA_character_

cat("Formatting annotations to compute TSS per gene symbol...\n")

# Use 'gene' features preferentially; fallback to 'transcript' if necessary
if ("type" %in% names(gene_anno_df)) {
  genes_df <- gene_anno_df %>% dplyr::filter(type == "gene")
  if (nrow(genes_df) == 0) {
    genes_df <- gene_anno_df %>% dplyr::filter(type == "transcript")
  }
} else {
  genes_df <- gene_anno_df
}

genes_df <- genes_df %>%
  dplyr::filter(!is.na(gene_name)) %>%
  mutate(
    TSS = ifelse(strand == "+", start, end),
    gene_symbol = gene_name
  ) %>%
  dplyr::select(gene_symbol, chromosome, TSS, strand)

cat(sprintf("  ✓ Gene model ready: %d TSS entries\n", nrow(genes_df)))

get_gene_tss <- function(df, symbol) {
  sub <- df[df$gene_symbol == symbol, ]
  if (nrow(sub) == 0) {
    return(NULL)
  }
  sub[1, c("chromosome", "TSS")]
}

# ---------------------------------------------------------------------------
# STEP 3: Define target genes from CLI
# ---------------------------------------------------------------------------

gene_symbols <- strsplit(opt$genes, ",")[[1]]
gene_symbols <- unique(trimws(gene_symbols))
cat("\nTarget genes requested: ", paste(gene_symbols, collapse = ", "), "\n", sep = "")

target_genes_list <- list()
for (g in gene_symbols) {
  tss_info <- get_gene_tss(genes_df, g)
  if (is.null(tss_info)) {
    cat(sprintf("  ! No TSS found for gene '%s' in annotations; skipping\n", g))
    next
  }
  target_genes_list[[length(target_genes_list) + 1]] <- data.frame(
    gene = g,
    chr  = as.character(tss_info$chromosome),
    tss  = as.numeric(tss_info$TSS),
    stringsAsFactors = FALSE
  )
}

if (length(target_genes_list) == 0) {
  stop("No valid target genes found in annotations; nothing to plot.", call. = FALSE)
}

target_genes <- do.call(rbind, target_genes_list)

cat("Resolved target genes with TSS:\n")
print(target_genes)

# ---------------------------------------------------------------------------
# STEP 4: Identify promoter peaks for each target (viewpoints)
# ---------------------------------------------------------------------------

cat("\nIdentifying promoter peaks (viewpoints) for each target gene...\n")

if (!exists("target_genes")) {
  stop("Internal error: target_genes was not created before promoter peak search", call. = FALSE)
}

# Parse peak strings in either format:
#   - "chr10:100-200" (Cicero connections)
#   - "chr10_100_200" (some CDS/peak matrices)
parse_peak_coords <- function(peak) {
  peak <- as.character(peak)

  if (grepl(":", peak, fixed = TRUE) && grepl("-", peak, fixed = TRUE)) {
    # chr:start-end
    sp <- strsplit(peak, ":", fixed = TRUE)[[1]]
    chr <- sp[1]
    se <- strsplit(sp[2], "-", fixed = TRUE)[[1]]
    if (length(se) != 2) return(NULL)
    return(list(chr = chr, start = as.numeric(se[1]), end = as.numeric(se[2])))
  }

  if (grepl("_", peak, fixed = TRUE)) {
    # chr_start_end
    parts <- strsplit(peak, "_", fixed = TRUE)[[1]]
    if (length(parts) != 3) return(NULL)
    return(list(chr = parts[1], start = as.numeric(parts[2]), end = as.numeric(parts[3])))
  }

  NULL
}

# Use peaks from connections so viewpoint strings match conns$Peak1/Peak2
all_conn_peaks <- unique(c(conns$Peak1, conns$Peak2))
cat(sprintf("  Connections peak namespace: %d unique peaks\n", length(all_conn_peaks)))
cat(sprintf("  Example peaks: %s\n", paste(head(all_conn_peaks, 3), collapse = ", ")))

# Precompute coordinates for all peaks once (used for debugging + fast window queries)
cat("  Parsing peak coordinates for connection peak namespace...\n")
peak_list <- lapply(all_conn_peaks, parse_peak_coords)
ok <- !vapply(peak_list, is.null, logical(1))
if (!any(ok)) {
  stop("Could not parse any peaks from connections (unexpected peak format).", call. = FALSE)
}

peak_df <- data.frame(
  peak  = all_conn_peaks[ok],
  chr   = vapply(peak_list[ok], `[[`, character(1), "chr"),
  start = vapply(peak_list[ok], `[[`, numeric(1), "start"),
  end   = vapply(peak_list[ok], `[[`, numeric(1), "end"),
  stringsAsFactors = FALSE
)
peak_df$mid <- (peak_df$start + peak_df$end) / 2

# Named coordinate maps for fast lookup
peak_chr   <- setNames(peak_df$chr,   peak_df$peak)
peak_start <- setNames(peak_df$start, peak_df$peak)
peak_end   <- setNames(peak_df$end,   peak_df$peak)

cat(sprintf("  ✓ Parsed %d/%d peaks\n", nrow(peak_df), length(all_conn_peaks)))
cat(sprintf("  ✓ Unique chromosomes in connections: %d\n", length(unique(peak_df$chr))))
cat(sprintf("  ✓ Example chromosomes: %s\n", paste(head(sort(unique(peak_df$chr))), collapse = ", ")))

# Find a promoter peak in ±search_window bp around TSS.
find_promoter_peak <- function(chr, tss, peak_names, search_window = 2000) {
  chr_peaks <- peak_names[peak_chr[peak_names] == chr]
  if (length(chr_peaks) == 0) return(NA_character_)

  best_peak <- NA_character_
  best_dist <- Inf

  for (peak in chr_peaks) {
    ps <- peak_start[[peak]]
    pe <- peak_end[[peak]]
    if (is.null(ps) || is.null(pe) || is.na(ps) || is.na(pe)) next

    if ((ps <= tss + search_window) && (pe >= tss - search_window)) {
      mid  <- (ps + pe) / 2
      dist <- abs(mid - tss)
      if (dist < best_dist) {
        best_dist <- dist
        best_peak <- peak
      }
    }
  }

  best_peak
}

target_genes$promoter_peak <- NA_character_
for (i in seq_len(nrow(target_genes))) {
  gene_name <- target_genes$gene[i]
  chr       <- target_genes$chr[i]
  tss       <- target_genes$tss[i]

  pp <- find_promoter_peak(chr, tss, all_conn_peaks, search_window = 2000)
  target_genes$promoter_peak[i] <- ifelse(is.na(pp), NA_character_, pp)

  cat(sprintf(
    "  %s promoter peak: %s\n",
    gene_name,
    ifelse(is.na(pp) || !nzchar(pp), "NOT FOUND", pp)
  ))
}

# Helper: peaks overlapping a window on a chromosome
peaks_in_window <- function(chr, minbp, maxbp) {
  peak_df$peak[peak_df$chr == chr & peak_df$end >= minbp & peak_df$start <= maxbp]
}

# Helper: connections that have BOTH peaks on chr and overlapping the window, and coaccess >= cutoff
connections_in_window <- function(chr, minbp, maxbp, cutoff) {
  # restrict by coaccess first
  idx <- which(conns$coaccess >= cutoff)
  if (length(idx) == 0) return(conns[0, , drop = FALSE])

  p1 <- conns$Peak1[idx]
  p2 <- conns$Peak2[idx]

  c1 <- peak_chr[p1]
  c2 <- peak_chr[p2]
  s1 <- peak_start[p1]; e1 <- peak_end[p1]
  s2 <- peak_start[p2]; e2 <- peak_end[p2]

  ok1 <- !is.na(c1) & c1 == chr & !is.na(s1) & !is.na(e1) & (e1 >= minbp) & (s1 <= maxbp)
  ok2 <- !is.na(c2) & c2 == chr & !is.na(s2) & !is.na(e2) & (e2 >= minbp) & (s2 <= maxbp)

  sub_idx <- idx[ok1 & ok2]
  if (length(sub_idx) == 0) return(conns[0, , drop = FALSE])

  conns[sub_idx, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# STEP 5: Generate connection plots per gene and threshold setting
# ---------------------------------------------------------------------------

cat("\nGenerating connection plots for target genes...\n")

outdir <- opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(target_genes))) {
  gene_name      <- target_genes$gene[i]
  chr            <- target_genes$chr[i]
  tss            <- target_genes$tss[i]
  promoter_peak  <- target_genes$promoter_peak[i]

  cat(sprintf("\n--- Plotting %s ---\n", gene_name))

  for (setting in plot_settings) {
    window    <- setting$window
    threshold <- setting$threshold
    suffix    <- setting$suffix

    start_win <- max(1, tss - window)
    end_win   <- tss + window

    cat(sprintf("  Setting: %s (window=±%dkb, threshold=%.3f)\n",
                suffix, window / 1000, threshold))
    cat(sprintf("    Region: %s:%d-%d\n", chr, start_win, end_win))

    # Data availability diagnostics
    win_peaks <- peaks_in_window(chr, start_win, end_win)
    n_win_peaks <- length(win_peaks)

    gm_rows <- 0
    if (all(c("chr","start","end") %in% names(gene_model))) {
      gm_rows <- sum(gene_model$chr == chr & gene_model$end >= start_win & gene_model$start <= end_win, na.rm = TRUE)
    }

    sub_conns <- connections_in_window(chr, start_win, end_win, threshold)
    n_sub <- nrow(sub_conns)

    cat(sprintf("    Diagnostics: peaks_in_window=%d | gene_model_rows_in_window=%d | conns_in_window_at_cutoff=%d\n",
                n_win_peaks, gm_rows, n_sub))

    if (n_sub > 0) {
      cat(sprintf("    Coaccess stats (in-window): min=%.3f mean=%.3f max=%.3f\n",
                  min(sub_conns$coaccess), mean(sub_conns$coaccess), max(sub_conns$coaccess)))
      cat(sprintf("    Example connections:\n"))
      print(utils::head(sub_conns, 3))
    } else {
      cat("    NOTE: 0 connections in this window at this cutoff -> plot will be empty by biology/data, not rendering.\n")
    }

    # Always write a small diagnostic TSV so you can inspect after Nextflow publishDir
    diag_tsv <- sprintf("%s/%s_debug_%s.tsv", outdir, gene_name, suffix)
    diag_df <- data.frame(
      gene = gene_name,
      chr = chr,
      tss = tss,
      window_bp = window,
      start = start_win,
      end = end_win,
      cutoff = threshold,
      peaks_in_window = n_win_peaks,
      gene_model_rows_in_window = gm_rows,
      conns_in_window = n_sub,
      stringsAsFactors = FALSE
    )
    write.table(diag_df, diag_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    # Plot WITHOUT viewpoint
    pdf_name_no_vp <- sprintf("%s/%s_connections_%s_noviewpoint.pdf",
                              outdir, gene_name, suffix)

    cat("    Plotting WITHOUT viewpoint... ")
    pdf(pdf_name_no_vp, width = 12, height = 6)

    if (n_sub == 0) {
      # Write an explicit “no data” diagnostic page
      plot.new()
      text(0.5, 0.65, sprintf("%s (%s)", gene_name, suffix), cex = 1.3)
      text(0.5, 0.50, sprintf("%s:%d-%d", chr, start_win, end_win), cex = 1.1)
      text(0.5, 0.35, sprintf("NO CONNECTIONS at coaccess_cutoff >= %.3f", threshold), cex = 1.2)
      cat("✓ (diagnostic page)\n")
      dev.off()
    } else {
      tryCatch({
        plot_connections(
          conns,
          chr = chr,
          minbp = start_win,
          maxbp = end_win,
          gene_model = gene_model,
          coaccess_cutoff = threshold,
          connection_width = 0.5,
          collapseTranscripts = "longest"
        )
        cat("✓\n")
      }, error = function(e) {
        cat(sprintf("✗ Error: %s\n", e$message))
      })
      dev.off()
    }

    # Plot WITH viewpoint (if promoter peak exists)
    if (!is.na(promoter_peak) && nzchar(promoter_peak)) {
      pdf_name_vp <- sprintf("%s/%s_connections_%s_withviewpoint.pdf",
                             outdir, gene_name, suffix)

      cat(sprintf("    Plotting WITH viewpoint (%s)... ", promoter_peak))
      pdf(pdf_name_vp, width = 12, height = 6)

      if (n_sub == 0) {
        plot.new()
        text(0.5, 0.65, sprintf("%s (%s)", gene_name, suffix), cex = 1.3)
        text(0.5, 0.50, sprintf("%s:%d-%d", chr, start_win, end_win), cex = 1.1)
        text(0.5, 0.35, sprintf("NO CONNECTIONS at coaccess_cutoff >= %.3f", threshold), cex = 1.2)
        text(0.5, 0.20, sprintf("Viewpoint requested: %s", promoter_peak), cex = 0.9)
        cat("✓ (diagnostic page)\n")
        dev.off()
      } else {
        tryCatch({
          plot_connections(
            conns,
            chr = chr,
            minbp = start_win,
            maxbp = end_win,
            gene_model = gene_model,
            coaccess_cutoff = threshold,
            connection_width = 0.5,
            collapseTranscripts = "longest",
            viewpoint = promoter_peak
          )
          cat("✓\n")
        }, error = function(e) {
          cat(sprintf("✗ Error: %s\n", e$message))
        })
        dev.off()
      }
    } else {
      cat("    Skipping viewpoint plot (no promoter peak found)\n")
    }
  }
}

# ============================================================
# STEP X: Summary statistics per gene (optional)
# ============================================================

cat("\nSummary statistics per target gene (±250kb window)...\n")

default_window <- 250000
default_cutoff <- if (length(thresholds) >= 2 && !is.na(thresholds[2])) thresholds[2] else 0.05

summ_list <- list()

# Use peaks from connections for consistent parsing
all_peaks <- unique(c(conns$Peak1, conns$Peak2))

for (i in seq_len(nrow(target_genes))) {
  gene_name <- target_genes$gene[i]
  chr       <- target_genes$chr[i]
  tss       <- target_genes$tss[i]

  start_win <- max(1, tss - default_window)
  end_win   <- tss + default_window

  # Peaks on this chr and overlapping the window
  window_peaks <- character(0)
  for (pk in all_peaks) {
    coords <- parse_peak_coords(pk)
    if (is.null(coords)) next
    if (!identical(coords$chr, chr)) next

    overlaps <- (coords$start <= end_win) && (coords$end >= start_win)
    if (overlaps) {
      window_peaks <- c(window_peaks, pk)
    }
  }
  window_peaks <- unique(window_peaks)

  gene_conns <- conns[
    (conns$Peak1 %in% window_peaks | conns$Peak2 %in% window_peaks) &
      conns$coaccess >= default_cutoff,
  ]

  summ_list[[length(summ_list) + 1]] <- data.frame(
    gene            = gene_name,
    chr             = chr,
    tss             = tss,
    window_bp       = default_window,
    cutoff          = default_cutoff,
    n_peaks_window  = length(window_peaks),
    n_conns_window  = nrow(gene_conns),
    mean_coaccess   = if (nrow(gene_conns) > 0) mean(gene_conns$coaccess) else NA_real_,
    max_coaccess    = if (nrow(gene_conns) > 0) max(gene_conns$coaccess) else NA_real_,
    stringsAsFactors = FALSE
  )
}

if (length(summ_list) > 0) {
  summary_df <- do.call(rbind, summ_list)
  summary_path <- file.path(opt$outdir, "cicero_target_gene_summary.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  cat("  ✓ Wrote summary to ", summary_path, "\n", sep = "")
}

plot_files <- list.files(outdir, pattern = "*.pdf", full.names = FALSE)
cat("Generated files:\n")
for (f in plot_files) {
  cat(sprintf("  • %s\n", f))
}
cat("\n")
