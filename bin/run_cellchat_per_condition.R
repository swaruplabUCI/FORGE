#!/usr/bin/env Rscript
# ===========================================================================
# CellChat Per-Condition Analysis (FIX-29)
# Runs the full CellChat L-R inference pipeline on cells from ONE condition.
# Produces a serialized CellChat object (.rds) and publication-quality plots.
#
# Usage:
#   run_cellchat_per_condition.R \
#     --h5ad integrated.h5ad \
#     --condition "90plus" \
#     --cell_type_key scanvi_prediction \
#     --condition_key condition_group \
#     --species human \
#     --output_rds cellchat_90plus.rds \
#     --output_dir cellchat_90plus_plots
# ===========================================================================

suppressPackageStartupMessages({
    library(CellChat)
    library(zellkonverter)
    library(SingleCellExperiment)
    library(Matrix)
    library(patchwork)
    library(ggplot2)
    library(optparse)
})

# Exact integral test, chunked with early exit: a normalised matrix trips on
# the first chunk, and a 909M-nnz raw matrix never materialises a full logical
# vector. Mirrors the helper in bin/run_cellchat_infer.R -- the bin/ scripts
# are standalone by design and share no library.
is_integral_matrix <- function(m, chunk = 1e7) {
    x <- m@x
    n <- length(x)
    if (n == 0) return(TRUE)
    starts <- seq(1, n, by = chunk)
    for (s in starts) {
        e <- min(s + chunk - 1, n)
        if (any(x[s:e] != floor(x[s:e]))) return(FALSE)
    }
    TRUE
}

# CPM + log1p, written as one in-place pass over the dgCMatrix @x.
# log1p(t(t(m) / ls) * 1e4) materialises two full transposes; on a large
# object that is tens of GB of temporaries and the process is OOM-killed with
# no R traceback. A dgCMatrix is column-compressed, so every non-zero's column
# -- and therefore its library size -- is recoverable from @p. Structural
# zeros are untouched, which is correct: log1p(0) = 0.
normalise_cpm_log1p <- function(m) {
    if (!inherits(m, "dgCMatrix"))
        m <- as(as(m, "CsparseMatrix"), "dgCMatrix")
    library.size <- Matrix::colSums(m)
    if (any(library.size == 0))
        stop("Cannot normalise: ", sum(library.size == 0),
             " cells have zero total counts")
    m@x <- log1p(m@x / rep.int(library.size, diff(m@p)) * 1e4)
    m
}

## Match run_cellchat.R: allow large globals for future (CellChat parallelization)
options(future.globals.maxSize = +Inf)

# --- CLI args ---------------------------------------------------------------
option_list <- list(
    make_option("--h5ad",          type = "character"),
    make_option("--condition",     type = "character"),
    make_option("--cell_type_key", type = "character", default = "scanvi_prediction"),
    make_option("--condition_key", type = "character", default = "condition_group"),
    make_option("--species",       type = "character", default = "human"),
    make_option("--output_rds",    type = "character", default = "cellchat.rds"),
    make_option("--output_dir",    type = "character", default = "cellchat_plots"),
    make_option("--group_mapping", type = "character", default = NULL,
                help = "JSON file mapping sample_id → condition_group (used if condition_key column is missing from h5ad)")
)
opts <- parse_args(OptionParser(option_list = option_list))

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load & subset ----------------------------------------------------------
cat("Loading h5ad via zellkonverter...\n")

# Read h5ad as SingleCellExperiment (no Seurat conversion — matches run_cellchat.R)
sce <- readH5AD(opts$h5ad)
meta <- as.data.frame(colData(sce))

cat(sprintf("Total cells: %d\n", ncol(sce)))

# FIX-46: If condition_key column is missing, merge from group_mapping JSON
if (!opts$condition_key %in% colnames(meta)) {
    if (is.null(opts$group_mapping) || !file.exists(opts$group_mapping)) {
        stop(sprintf("Column '%s' not found in h5ad and no --group_mapping provided.",
                     opts$condition_key))
    }
    cat(sprintf("Column '%s' not found — merging from group_mapping JSON...\n", opts$condition_key))
    mapping <- jsonlite::fromJSON(opts$group_mapping)
    sample_col <- if ("sample" %in% colnames(meta)) "sample"
                  else if ("sample_id" %in% colnames(meta)) "sample_id"
                  else if ("batch" %in% colnames(meta)) "batch"
                  else stop("Cannot find sample/sample_id/batch column to join group mapping")
    meta[[opts$condition_key]] <- mapping[meta[[sample_col]]]
    assigned <- sum(!is.na(meta[[opts$condition_key]]))
    cat(sprintf("  Mapped %d/%d cells via '%s' column\n", assigned, ncol(sce), sample_col))
}

cat(sprintf("Subsetting to condition: %s (key: %s)\n", opts$condition, opts$condition_key))

# Subset to the requested condition
cells_keep <- which(meta[[opts$condition_key]] == opts$condition)
if (length(cells_keep) < 50) {
    stop(sprintf("Too few cells (%d) for condition '%s'. Minimum 50 required.",
                 length(cells_keep), opts$condition))
}
sce_sub <- sce[, cells_keep]
meta <- meta[cells_keep, , drop = FALSE]
meta[[opts$cell_type_key]] <- droplevels(factor(meta[[opts$cell_type_key]]))
cat(sprintf("Cells after subsetting: %d\n", ncol(sce_sub)))

# --- Extract expression data (matches run_cellchat.R pattern) ----------------
# zellkonverter surfaces h5ad layers as additional assays, so `logcounts` here
# is layers/logcounts. The seurat_rds on-ramp writes X = RAW COUNTS plus
# layers/logcounts = the user's own normalised data; preferring that layer
# keeps this consumer consistent with run_cellchat_infer.R and uses the
# normalisation the user actually chose rather than re-deriving one.
#
# GUARDED on the layer's presence, so the on-ramp and native paths each take
# the branch that is correct for them.
#
# 2026-10-30: the `X` branch is now GUARDED ON INTEGRALITY, which is what
# run_cellchat_infer.R has done since 4cf3975. It previously applied
# CPM + log1p unconditionally, and on the native convention X is ALREADY
# log-normalised -- so it normalised data that was already normalised. Worse,
# colSums() of log-normalised values is not a library size at all; it is a sum
# of logs, so every cell was divided by an arbitrary per-cell constant
# (measured spread on real data: 2.13 to 19.25) and log1p'd again.
#
# Resolution order, now identical in all three CellChat readers:
#   1. scvi_normalized, when present
#   2. logcounts, when the producer supplied it        (the seurat_rds on-ramp)
#   3. X as-is, when X is already non-integral          (the native path)
#   4. CPM + log1p, when X is integral and nothing else is available
if ("scvi_normalized" %in% assayNames(sce_sub)) {
    data.input <- assay(sce_sub, "scvi_normalized")
} else if ("logcounts" %in% assayNames(sce_sub)) {
    cat("Found logcounts assay -> using it as the normalised matrix\n")
    data.input <- assay(sce_sub, "logcounts")
} else if ("X" %in% assayNames(sce_sub)) {
    data.input <- assay(sce_sub, "X")
    if (is_integral_matrix(data.input)) {
        cat("X is integral (raw counts) -> normalising with CPM + log1p\n")
        data.input <- normalise_cpm_log1p(data.input)
    } else {
        cat("X is already non-integral -> using it as the normalised matrix",
            "(native-path convention)\n")
    }
} else {
    stop("Cannot find appropriate expression matrix in h5ad")
}

# --- Create CellChat object -------------------------------------------------
cat("Creating CellChat object...\n")

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta,
                           group.by = opts$cell_type_key)

# Set the L-R database (matches run_cellchat.R — full DB, not subset)
if (opts$species == "human") {
    CellChatDB <- CellChatDB.human
} else {
    CellChatDB <- CellChatDB.mouse
}
cellchat@DB <- CellChatDB

# --- Core L-R inference pipeline (matches run_cellchat.R) -------------------
cat("Running CellChat L-R inference...\n")

cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Pathway-level aggregation
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# --- Centrality analysis ----------------------------------------------------
cat("Computing network centrality...\n")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# --- Communication patterns -------------------------------------------------
cat("Identifying communication patterns...\n")
tryCatch({
    # Determine optimal k via NMF
    nPatterns_out <- selectK(cellchat, pattern = "outgoing")
    nPatterns_in  <- selectK(cellchat, pattern = "incoming")

    # Use k=4 as a reasonable default (capped at min cell types - 1)
    n_ct <- length(levels(cellchat@idents))
    k_out <- min(4, max(2, n_ct - 1))
    k_in  <- min(4, max(2, n_ct - 1))

    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = k_out)
    cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = k_in)
}, error = function(e) {
    cat(sprintf("Warning: Pattern analysis failed: %s\n", e$message))
})

# --- Visualization ----------------------------------------------------------
cat("Generating plots...\n")

# 1. Aggregated network circle plot
pdf(file.path(opts$output_dir, "network_circle.pdf"), width = 10, height = 10)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = TRUE, label.edge = FALSE,
                 title.name = paste0("Number of interactions - ", opts$condition))
dev.off()

# 2. Signaling role heatmap
pdf(file.path(opts$output_dir, "signaling_role_heatmap.pdf"), width = 12, height = 8)
tryCatch({
    netAnalysis_signalingRole_heatmap(cellchat, signaling = cellchat@netP$pathways,
                                      pattern = "outgoing", title = opts$condition)
}, error = function(e) cat(sprintf("Heatmap outgoing skipped: %s\n", e$message)))
tryCatch({
    netAnalysis_signalingRole_heatmap(cellchat, signaling = cellchat@netP$pathways,
                                      pattern = "incoming", title = opts$condition)
}, error = function(e) cat(sprintf("Heatmap incoming skipped: %s\n", e$message)))
dev.off()

# 3. Top pathways bubble plot
pdf(file.path(opts$output_dir, "bubble_top_pathways.pdf"), width = 14, height = 10)
tryCatch({
    netVisual_bubble(cellchat, remove.isolate = TRUE)
}, error = function(e) cat(sprintf("Bubble plot skipped: %s\n", e$message)))
dev.off()

# 4. Communication pattern river/dot plots
if (!is.null(cellchat@netP$pattern)) {
    pdf(file.path(opts$output_dir, "patterns_outgoing.pdf"), width = 12, height = 8)
    tryCatch({
        netAnalysis_river(cellchat, pattern = "outgoing")
        netAnalysis_dot(cellchat, pattern = "outgoing")
    }, error = function(e) cat(sprintf("Pattern outgoing plots skipped: %s\n", e$message)))
    dev.off()

    pdf(file.path(opts$output_dir, "patterns_incoming.pdf"), width = 12, height = 8)
    tryCatch({
        netAnalysis_river(cellchat, pattern = "incoming")
        netAnalysis_dot(cellchat, pattern = "incoming")
    }, error = function(e) cat(sprintf("Pattern incoming plots skipped: %s\n", e$message)))
    dev.off()
}

# --- Save -------------------------------------------------------------------
cat(sprintf("Saving CellChat object to %s\n", opts$output_rds))
saveRDS(cellchat, file = opts$output_rds)

cat("CellChat per-condition analysis complete.\n")
