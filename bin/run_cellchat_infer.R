#!/usr/bin/env Rscript

# CellChat Inference Script (no plotting)
# Compatible with h5ad input from scVI/SCANVI integration

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(patchwork)
    library(Seurat)
    library(schard)
    library(Matrix)
})

# ------------------------------
# Parse command line arguments
# ------------------------------
option_list <- list(
    make_option(c("--input"), type = "character", help = "Input h5ad file"),
    make_option(c("--output_prefix"), type = "character", help = "Output prefix"),
    make_option(c("--cell_type_key"), type = "character", default = "cell_type",
                help = "Cell type annotation key in h5ad"),
    make_option(c("--condition_key"), type = "character", default = "none",
                help = "Condition key for subsetting (optional)"),
    make_option(c("--species"), type = "character", default = "human",
                help = "Species: human or mouse"),
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## Allow large globals for future (CellChat parallelization)
options(future.globals.maxSize = +Inf)

# ------------------------------
# Load data from h5ad (using schard — no Python/basilisk needed)
# ------------------------------
message("Loading h5ad file: ", opt$input)
seurat_obj <- schard::h5ad2seurat(opt$input)

# Extract normalized data (schard puts X into 'data' layer)
assay_obj <- seurat_obj@assays[[DefaultAssay(seurat_obj)]]
available_layers <- Layers(assay_obj)
if ("data" %in% available_layers) {
    data.input <- GetAssayData(seurat_obj, layer = "data")
} else if ("counts" %in% available_layers) {
    counts <- GetAssayData(seurat_obj, layer = "counts")
    library.size <- Matrix::colSums(counts)
    data.input <- log1p(Matrix::t(Matrix::t(counts) / library.size) * 10000)
} else {
    stop("Cannot find appropriate expression matrix in h5ad")
}

# Extract metadata
meta <- seurat_obj@meta.data

# Check if cell type key exists
if (!opt$cell_type_key %in% colnames(meta)) {
    stop(
        "Cell type key '", opt$cell_type_key,
        "' not found in metadata. Available: ",
        paste(colnames(meta), collapse = ", ")
    )
}

# Subset by condition if specified (hook for future use)
if (opt$condition_key != "none" && opt$condition_key %in% colnames(meta)) {
    message("Subsetting by condition: ", opt$condition_key)
    # User can insert condition-based subsetting here if desired
}

# ------------------------------
# Create CellChat object
# ------------------------------
message("Creating CellChat object...")
cellchat <- createCellChat(
    object = data.input,
    meta   = meta,
    group.by = opt$cell_type_key
)

# Set database based on species
message("Loading CellChatDB for ", opt$species, "...")
if (opt$species == "human") {
    CellChatDB <- CellChatDB.human
} else if (opt$species == "mouse") {
    CellChatDB <- CellChatDB.mouse
} else {
    stop("Species must be 'human' or 'mouse'")
}
cellchat@DB <- CellChatDB

# ------------------------------
# Preprocessing and inference
# ------------------------------
message("Preprocessing data...")
cellchat <- subsetData(cellchat)

# FIX-78: Limit parallel workers to 4 — computeCommunProb() does sparse->dense
# coercions (~3 GB each); too many concurrent workers can exhaust memory
future::plan("multisession", workers = min(opt$threads, 4))

message("Identifying over-expressed genes and interactions...")
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

message("Computing communication probabilities...")
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

message("Computing pathway-level communication...")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

message("Performing network analysis...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# ------------------------------
# Export results (RDS + CSV)
# ------------------------------
message("Exporting results (RDS + CSV)...")

rds_file <- paste0(opt$output_prefix, "_cellchat.rds")
csv_file <- paste0(opt$output_prefix, "_cellchat_results.csv")

saveRDS(cellchat, file = rds_file)

df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = csv_file, row.names = FALSE)

message("CellChat inference complete!")
message("  RDS: ", rds_file)
message("  CSV: ", csv_file)
