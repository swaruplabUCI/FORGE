#!/usr/bin/env Rscript

# CellChat Analysis Script for Single Dataset
# Compatible with h5ad input from scVI/SCANVI integration

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(patchwork)
    library(zellkonverter)
    library(SingleCellExperiment)
    library(Matrix)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--input"), type="character", help="Input h5ad file"),
    make_option(c("--output_prefix"), type="character", help="Output prefix"),
    make_option(c("--cell_type_key"), type="character", default="cell_type", 
                help="Cell type annotation key in h5ad"),
    make_option(c("--condition_key"), type="character", default="none",
                help="Condition key for subsetting (optional)"),
    make_option(c("--species"), type="character", default="human",
                help="Species: human or mouse"),
    make_option(c("--threads"), type="integer", default=4, help="Number of threads")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Setup output directory
plot_dir <- paste0(opt$output_prefix, "_cellchat_plots")
dir.create(plot_dir, showWarnings = FALSE)

## Allow large globals for future (CellChat parallelization)
## FIX-58: Raised from 20 GB to unlimited — 287K cell dataset exceeds 20 GB
## (error: 20.43 GiB > 20.00 GiB, driven by 9.42 GiB ...future.FUN closure)
options(future.globals.maxSize = +Inf)

# Load data from h5ad
message("Loading h5ad file: ", opt$input)
sce <- readH5AD(file = opt$input)

# Extract normalized counts
# scVI/SCANVI typically stores normalized data in layers
if ("scvi_normalized" %in% assayNames(sce)) {
    data.input <- assay(sce, "scvi_normalized")
} else if ("X" %in% assayNames(sce)) {
    # If normalized data not available, normalize counts
    counts <- assay(sce, "X")
    library.size <- Matrix::colSums(counts)
    data.input <- log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000)
} else {
    stop("Cannot find appropriate expression matrix in h5ad")
}

# Extract metadata
meta <- as.data.frame(colData(sce))

# Check if cell type key exists
if (!opt$cell_type_key %in% colnames(meta)) {
    stop("Cell type key '", opt$cell_type_key, "' not found in metadata. Available: ", 
         paste(colnames(meta), collapse=", "))
}

# Subset by condition if specified
if (opt$condition_key != "none" && opt$condition_key %in% colnames(meta)) {
    message("Subsetting by condition: ", opt$condition_key)
    # User can modify this to select specific conditions
}

# Create CellChat object
message("Creating CellChat object...")
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = opt$cell_type_key)

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

# Preprocessing
message("Preprocessing data...")
cellchat <- subsetData(cellchat)

# Set up parallel processing
future::plan("multisession", workers = opt$threads)

# Identify over-expressed genes and interactions
message("Identifying over-expressed genes and interactions...")
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities
message("Computing communication probabilities...")
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer pathway-level communication
message("Computing pathway-level communication...")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Network analysis
message("Performing network analysis...")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Generate plots
message("Generating visualizations...")

# 1. Overview plots
pdf(file.path("cellchat_plots", paste0(opt$output_prefix, "_overview.pdf")), 
    width = 12, height = 6)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = TRUE, label.edge = FALSE, 
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = TRUE, label.edge = FALSE, 
                 title.name = "Interaction weights/strength")
dev.off()

# 2. Signaling role analysis
pdf(file.path("cellchat_plots", paste0(opt$output_prefix, "_signaling_roles.pdf")),
    width = 10, height = 8)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
gg3 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
print(gg1)
print(gg2)
print(gg3)
dev.off()

# 3. Top pathways
pathways.show <- cellchat@netP$pathways
if (length(pathways.show) > 0) {
    # Plot top 10 pathways
    n_pathways <- min(10, length(pathways.show))
    for (i in 1:n_pathways) {
        pathway <- pathways.show[i]

        ## Hierarchy plot – can trigger netVisual_hierarchy1 bug; protect with tryCatch
        tryCatch({
            pdf(file.path("cellchat_plots", 
                          paste0(opt$output_prefix, "_", pathway, "_hierarchy.pdf")),
                width = 10, height = 8)
            netVisual_aggregate(cellchat, signaling = pathway, layout = "hierarchy")
            dev.off()
        }, error = function(e) {
            message("Skipping hierarchy plot for pathway ", pathway,
                    " due to error in netVisual_aggregate/netVisual_hierarchy1: ", e$message)
            # Close any open device just in case
            if (names(dev.cur()) != "null device") dev.off()
        })

        ## Circle plot – usually safe
        pdf(file.path("cellchat_plots", 
                      paste0(opt$output_prefix, "_", pathway, "_circle.pdf")),
            width = 8, height = 8)
        netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
        dev.off()
    }
}

# Export results
message("Exporting results...")

# Save CellChat object
saveRDS(cellchat, file = paste0(opt$output_prefix, "_cellchat.rds"))

# Export communications as CSV
df.net <- subsetCommunication(cellchat)
write.csv(df.net, file = paste0(opt$output_prefix, "_cellchat_results.csv"), 
          row.names = FALSE)

# Generate HTML report
message("Generating HTML report...")
tryCatch({
    rmarkdown::render(
        input = system.file("templates", "cellchat_report.Rmd", package = "CellChat"),
        output_file = paste0(opt$output_prefix, "_cellchat_report.html"),
        params = list(cellchat_object = cellchat, prefix = opt$output_prefix)
    )
}, error = function(e) {
    message("HTML CellChat report generation failed, likely due to a plotting bug in CellChat (e.g. netVisual_hierarchy1): ", e$message)
    message("CellChat object (.rds) and CSV results were saved successfully; skipping HTML report.")
})

message("CellChat analysis complete!")
