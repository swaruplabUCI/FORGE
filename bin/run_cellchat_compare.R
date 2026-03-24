#!/usr/bin/env Rscript
# ===========================================================================
# CellChat Comparative Analysis (FIX-29)
# Merges per-condition CellChat objects and performs comparative analysis:
#   - Differential interaction counts and weights
#   - Information flow comparison
#   - Differential signaling pathway analysis
#   - Joint manifold learning (UMAP) of signaling pathways
#
# Usage:
#   run_cellchat_compare.R \
#     --rds_files "cellchat_90plus.rds,cellchat_Control.rds" \
#     --conditions "90plus,Control" \
#     --species human \
#     --output_rds cellchat_comparison.rds \
#     --output_dir comparison_plots
# ===========================================================================

suppressPackageStartupMessages({
    library(CellChat)
    library(patchwork)
    library(ggplot2)
    library(optparse)
})

# --- CLI args ---------------------------------------------------------------
option_list <- list(
    make_option("--rds_files",   type = "character"),
    make_option("--conditions",  type = "character"),
    make_option("--species",     type = "character", default = "human"),
    make_option("--output_rds",  type = "character", default = "cellchat_comparison.rds"),
    make_option("--output_dir",  type = "character", default = "comparison_plots")
)
opts <- parse_args(OptionParser(option_list = option_list))

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load per-condition CellChat objects ------------------------------------
rds_files  <- strsplit(opts$rds_files, ",")[[1]]
conditions <- strsplit(opts$conditions, ",")[[1]]

cat(sprintf("Loading %d CellChat objects...\n", length(rds_files)))
object.list <- lapply(rds_files, readRDS)
names(object.list) <- conditions

# Lift each object to a common cell type set
cat("Lifting cell types to shared set...\n")
cellchat_list <- lapply(object.list, function(cc) {
    tryCatch({
        liftCellChat(cc, group.by = "ident")
    }, error = function(e) {
        cat(sprintf("liftCellChat note: %s — using original object\n", e$message))
        cc
    })
})

# --- Merge ------------------------------------------------------------------
cat("Merging CellChat objects...\n")
cellchat_merged <- mergeCellChat(cellchat_list, add.names = conditions)

# --- Comparative Visualizations ---------------------------------------------
cat("Generating comparative plots...\n")

# 1. Differential number of interactions & interaction strength
pdf(file.path(opts$output_dir, "interaction_comparison.pdf"), width = 14, height = 7)
tryCatch({
    p1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
    p2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2),
                               measure = "weight")
    print(p1 + p2)
}, error = function(e) cat(sprintf("Interaction comparison skipped: %s\n", e$message)))
dev.off()

# 2. Differential circle plots (net change)
pdf(file.path(opts$output_dir, "differential_circle.pdf"), width = 14, height = 7)
tryCatch({
    par(mfrow = c(1, 2))
    weight.max <- getMaxWeight(cellchat_list, attribute = c("idents", "count"))
    for (i in seq_along(cellchat_list)) {
        netVisual_circle(cellchat_list[[i]]@net$count,
                         weight.scale = TRUE, label.edge = FALSE,
                         edge.weight.max = weight.max[2],
                         title.name = paste0("Number of interactions - ", conditions[i]))
    }
    par(mfrow = c(1, 1))
}, error = function(e) cat(sprintf("Circle plots skipped: %s\n", e$message)))
dev.off()

# 3. Information flow ranking comparison
pdf(file.path(opts$output_dir, "information_flow.pdf"), width = 10, height = 12)
tryCatch({
    p <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
    print(p)
}, error = function(e) cat(sprintf("Information flow skipped: %s\n", e$message)))
dev.off()

# 4. Outgoing/Incoming signaling pattern comparison (heatmaps)
pdf(file.path(opts$output_dir, "signaling_changes_heatmap.pdf"), width = 14, height = 10)
tryCatch({
    # Functional similarity
    cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
    cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
    cellchat_merged <- netClustering(cellchat_merged, type = "functional")
    netVisual_embeddingPairwise(cellchat_merged, type = "functional",
                                label.size = 3.5, title = "Functional similarity")
}, error = function(e) cat(sprintf("Functional embedding skipped: %s\n", e$message)))
tryCatch({
    # Structural similarity
    cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
    cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
    cellchat_merged <- netClustering(cellchat_merged, type = "structural")
    netVisual_embeddingPairwise(cellchat_merged, type = "structural",
                                label.size = 3.5, title = "Structural similarity")
}, error = function(e) cat(sprintf("Structural embedding skipped: %s\n", e$message)))
dev.off()

# 5. Specific signaling pathway changes (top differential pathways)
pdf(file.path(opts$output_dir, "pathway_specific_comparison.pdf"), width = 16, height = 10)
tryCatch({
    # Get shared pathways between conditions
    pathways1 <- cellchat_list[[1]]@netP$pathways
    pathways2 <- cellchat_list[[2]]@netP$pathways
    shared_pathways <- intersect(pathways1, pathways2)

    if (length(shared_pathways) > 0) {
        # Show top 6 most different pathways
        top_n <- min(6, length(shared_pathways))
        for (pathway in shared_pathways[1:top_n]) {
            tryCatch({
                netVisual_aggregate(cellchat_list[[1]], signaling = pathway,
                                    layout = "circle",
                                    title.name = paste0(pathway, " - ", conditions[1]))
                netVisual_aggregate(cellchat_list[[2]], signaling = pathway,
                                    layout = "circle",
                                    title.name = paste0(pathway, " - ", conditions[2]))
            }, error = function(e) NULL)
        }
    }
}, error = function(e) cat(sprintf("Pathway comparison skipped: %s\n", e$message)))
dev.off()

# 6. Outgoing/incoming signal changes per cell type
pdf(file.path(opts$output_dir, "signal_changes_scatter.pdf"), width = 12, height = 10)
tryCatch({
    num_cells <- length(levels(cellchat_list[[1]]@idents))
    for (i in seq_len(num_cells)) {
        tryCatch({
            netVisual_bubble(cellchat_merged, sources.use = i,
                           comparison = c(1, 2), max.dataset = 2,
                           remove.isolate = TRUE,
                           title.name = paste0("Outgoing from ", levels(cellchat_list[[1]]@idents)[i]))
        }, error = function(e) NULL)
    }
}, error = function(e) cat(sprintf("Signal scatter skipped: %s\n", e$message)))
dev.off()

# --- Save merged object -----------------------------------------------------
cat(sprintf("Saving comparison object to %s\n", opts$output_rds))
saveRDS(cellchat_merged, file = opts$output_rds)

cat("CellChat comparative analysis complete.\n")
