#!/usr/bin/env Rscript

# CellChat Comparison Script
# Compares multiple CellChat objects (e.g., conditions, timepoints)

suppressPackageStartupMessages({
    library(optparse)
    library(CellChat)
    library(patchwork)
})

option_list <- list(
    make_option(c("--input"), type="character", 
                help="Comma-separated list of CellChat RDS files"),
    make_option(c("--output_prefix"), type="character", help="Output prefix"),
    make_option(c("--threads"), type="integer", default=4, help="Number of threads")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

dir.create("comparison_plots", showWarnings = FALSE)

# Load CellChat objects
input_files <- strsplit(opt$input, ",")[[1]]
object.list <- lapply(input_files, readRDS)
names(object.list) <- gsub("_cellchat.rds", "", basename(input_files))

message("Loaded ", length(object.list), " CellChat objects: ", 
        paste(names(object.list), collapse=", "))

# Merge objects
message("Merging CellChat objects...")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Comparative analysis
message("Performing comparative analysis...")

# 1. Compare interaction counts
pdf(file.path("comparison_plots", paste0(opt$output_prefix, "_interaction_counts.pdf")),
    width=8, height=6)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
print(gg1 + gg2)
dev.off()

# 2. Differential interaction heatmap
pdf(file.path("comparison_plots", paste0(opt$output_prefix, "_diff_interactions.pdf")),
    width=12, height=10)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
print(gg1 + gg2)
dev.off()

# 3. Compare major sources and targets
pdf(file.path("comparison_plots", paste0(opt$output_prefix, "_sources_targets.pdf")),
    width=12, height=6)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], 
                                                   title = names(object.list)[i],
                                                   weight.MinMax = weight.MinMax)
}
print(patchwork::wrap_plots(plots = gg))
dev.off()

# Save merged object
saveRDS(cellchat, file = paste0(opt$output_prefix, "_comparison.rds"))

message("Comparative analysis complete!")
