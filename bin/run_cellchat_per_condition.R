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
    library(Seurat)
    library(zellkonverter)
    library(SingleCellExperiment)
    library(patchwork)
    library(ggplot2)
    library(optparse)
})

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

# Read h5ad as SingleCellExperiment, then convert to Seurat
sce <- readH5AD(opts$h5ad)
seurat_obj <- as.Seurat(sce, counts = "X", data = NULL)

cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))

# FIX-46: If condition_key column is missing, merge from group_mapping JSON
if (!opts$condition_key %in% colnames(seurat_obj@meta.data)) {
    if (is.null(opts$group_mapping) || !file.exists(opts$group_mapping)) {
        stop(sprintf("Column '%s' not found in h5ad and no --group_mapping provided.",
                     opts$condition_key))
    }
    cat(sprintf("Column '%s' not found — merging from group_mapping JSON...\n", opts$condition_key))
    mapping <- jsonlite::fromJSON(opts$group_mapping)
    sample_col <- if ("sample" %in% colnames(seurat_obj@meta.data)) "sample"
                  else if ("sample_id" %in% colnames(seurat_obj@meta.data)) "sample_id"
                  else if ("batch" %in% colnames(seurat_obj@meta.data)) "batch"
                  else stop("Cannot find sample/sample_id/batch column to join group mapping")
    seurat_obj@meta.data[[opts$condition_key]] <- mapping[seurat_obj@meta.data[[sample_col]]]
    assigned <- sum(!is.na(seurat_obj@meta.data[[opts$condition_key]]))
    cat(sprintf("  Mapped %d/%d cells via '%s' column\n", assigned, ncol(seurat_obj), sample_col))
}

cat(sprintf("Subsetting to condition: %s (key: %s)\n", opts$condition, opts$condition_key))

# Subset to the requested condition
cells_keep <- which(seurat_obj@meta.data[[opts$condition_key]] == opts$condition)
if (length(cells_keep) < 50) {
    stop(sprintf("Too few cells (%d) for condition '%s'. Minimum 50 required.",
                 length(cells_keep), opts$condition))
}
seurat_sub <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_keep])
cat(sprintf("Cells after subsetting: %d\n", ncol(seurat_sub)))

# --- Create CellChat object -------------------------------------------------
cat("Creating CellChat object...\n")

# Set cell type identities
Idents(seurat_sub) <- seurat_sub@meta.data[[opts$cell_type_key]]

# Get normalized data
data.input <- GetAssayData(seurat_sub, assay = "RNA", slot = "data")
meta <- seurat_sub@meta.data

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta,
                           group.by = opts$cell_type_key)

# Set the L-R database
if (opts$species == "human") {
    CellChatDB <- CellChatDB.human
} else {
    CellChatDB <- CellChatDB.mouse
}

# Use Secreted Signaling subset (most relevant for brain tissue)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# --- Core L-R inference pipeline --------------------------------------------
cat("Running CellChat L-R inference...\n")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
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
