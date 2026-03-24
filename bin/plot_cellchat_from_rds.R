#!/usr/bin/env Rscript

# Plot CellChat results from an existing .rds object
# Generates overview, signaling-role, and top-pathway circle plots

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(patchwork)
})

# ------------------------------
# Parse command line arguments
# ------------------------------
option_list <- list(
    make_option(c("--input_rds"), type = "character",
                help = "Input CellChat RDS file (e.g. integrated_cellchat.rds)"),
    make_option(c("--output_prefix"), type = "character",
                help = "Output prefix (used for naming plot directory)"),
    make_option(c("--max_pathways"), type = "integer", default = 10,
                help = "Maximum number of top pathways to plot")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_rds) || is.null(opt$output_prefix)) {
    stop("Both --input_rds and --output_prefix must be provided.")
}

# ------------------------------
# Check RDS existence and load
# ------------------------------
if (!file.exists(opt$input_rds)) {
    message("Input RDS not found (", opt$input_rds, "). Skipping plotting.")
    quit(save = "no", status = 0)
}

message("Loading CellChat object from: ", opt$input_rds)
cellchat <- readRDS(opt$input_rds)

# ------------------------------
# Setup output directory
# ------------------------------
plot_dir <- paste0(opt$output_prefix, "_cellchat_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
message("Plots will be written to: ", plot_dir)

# ------------------------------
# 1) Overview circle plots
# ------------------------------
message("Generating overview circle plots...")

groupSize <- as.numeric(table(cellchat@idents))

pdf(file.path(plot_dir, paste0(opt$output_prefix, "_overview.pdf")),
    width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(
    cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
)

netVisual_circle(
    cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights/strength"
)

dev.off()

# ------------------------------
# 2) Signaling role analysis
# ------------------------------
message("Generating signaling role plots...")

pdf(file.path(plot_dir, paste0(opt$output_prefix, "_signaling_roles.pdf")),
    width = 10, height = 8)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
gg3 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
print(gg1)
print(gg2)
print(gg3)
dev.off()

# ------------------------------
# 3) Top pathways – circle plots only
# ------------------------------
message("Generating top-pathway circle plots...")

pathways.show <- cellchat@netP$pathways
if (length(pathways.show) > 0) {
    n_pathways <- min(opt$max_pathways, length(pathways.show))

    for (i in seq_len(n_pathways)) {
        pathway <- pathways.show[i]
        pdf(file.path(plot_dir,
                      paste0(opt$output_prefix, "_", pathway, "_circle.pdf")),
            width = 8, height = 8)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        dev.off()
    }
} else {
    message("No pathways in cellchat@netP$pathways; skipping pathway plots.")
}

message("CellChat plotting complete.")
message("Plots written to: ", plot_dir)
