#!/usr/bin/env Rscript
# ===========================================================================
# hdWGCNA Enrichment & Network Visualization (FIX-29)
# Runs as part of the global (condition-agnostic) analysis tier:
#   1. Enrichr (GO Biological Process, Cellular Component, Molecular Function)
#   2. Module network plots (per-module hub gene networks)
#   3. Combined hub gene network
#   4. Co-expression UMAP (supervised + unsupervised)
#
# Expects an hdWGCNA-processed Seurat object with completed modules.
#
# Usage:
#   run_hdwgcna_enrichment.R \
#     --seurat_rds hdwgcna_INH.rds \
#     --cell_type INH \
#     --species human \
#     --output_prefix enrichment
# ===========================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(cowplot)
    library(patchwork)
    library(WGCNA)
    library(hdWGCNA)
    library(enrichR)
    library(igraph)
    library(optparse)
})

theme_set(theme_cowplot())
set.seed(12345)
# RunModuleUMAP exports per-cell distance matrices via the `future` parallel backend,
# which trips future.globals.maxSize (default 500 MiB) on cell types with many cells.
# Bump to 100 GB to match run_hdwgcna_celltype.R.
options(future.globals.maxSize = 100 * 1024^3)

# --- CLI args ---------------------------------------------------------------
option_list <- list(
    make_option("--seurat_rds",    type = "character"),
    make_option("--cell_type",     type = "character"),
    make_option("--species",       type = "character", default = "human"),
    make_option("--output_prefix", type = "character", default = "enrichment")
)
opts <- parse_args(OptionParser(option_list = option_list))

dir.create("enrichr_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("network_plots", recursive = TRUE, showWarnings = FALSE)

# --- Load -------------------------------------------------------------------
cat(sprintf("Loading Seurat object from %s\n", opts$seurat_rds))
seurat_obj <- readRDS(opts$seurat_rds)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
cat(sprintf("Modules found: %d\n", length(unique(modules$module))))
cat(sprintf("Total genes in modules: %d\n", nrow(modules)))

# --- Step 1: Enrichr -------------------------------------------------------
cat("Running Enrichr enrichment analysis...\n")

tryCatch({
    dbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023')

    seurat_obj <- RunEnrichr(seurat_obj, dbs = dbs, max_genes = 100)
    enrich_df <- GetEnrichrTable(seurat_obj)

    write.csv(enrich_df, "enrichr_results.csv", row.names = FALSE)
    cat(sprintf("Enrichr results: %d terms\n", nrow(enrich_df)))

    # Bar plots per module
    EnrichrBarPlot(
        seurat_obj,
        outdir = "enrichr_plots",
        n_terms = 10,
        plot_size = c(5, 7),
        logscale = TRUE
    )
    cat("Enrichr bar plots saved.\n")

    # Dot plot for biological process
    tryCatch({
        p_dot <- EnrichrDotPlot(
            seurat_obj,
            database = "GO_Biological_Process_2023",
            n_terms = 2
        )
        ggsave(file.path("enrichr_plots", "dotplot_BP.pdf"), p_dot, width = 14, height = 10)
    }, error = function(e) cat(sprintf("Enrichr dot plot failed: %s\n", e$message)))

}, error = function(e) {
    cat(sprintf("Enrichr failed: %s\n", e$message))
    cat("This may be due to network issues with the Enrichr API.\n")
})

# --- Step 2: Module Network Plots -------------------------------------------
cat("Generating module network plots...\n")

tryCatch({
    ModuleNetworkPlot(
        seurat_obj,
        outdir = "network_plots",
        n_inner = 10,
        n_outer = 15,
        vertex.label.cex = 1
    )
    cat("Module network plots saved.\n")
}, error = function(e) cat(sprintf("Module network plots failed: %s\n", e$message)))

# --- Step 3: Combined Hub Gene Network --------------------------------------
# HubGeneNetworkPlot() prints to the active graphics device and returns invisibly;
# wrapping with ggsave produces a blank ggplot. Use pdf() / dev.off() instead
# (the same pattern used by run_hdwgcna_celltype.R).
cat("Generating combined hub gene network...\n")

tryCatch({
    pdf(file.path("network_plots", "hub_gene_network.pdf"), width = 14, height = 12)
    HubGeneNetworkPlot(
        seurat_obj,
        n_hubs = 3,
        n_other = 5,
        edge_prop = 0.75,
        mods = 'all'
    )
    dev.off()
    cat("Hub gene network saved.\n")
}, error = function(e) {
    cat(sprintf("Hub gene network failed: %s\n", e$message))
    try(dev.off(), silent = TRUE)
})

# --- Step 4: Co-expression UMAP ---------------------------------------------
cat("Running co-expression UMAP...\n")

tryCatch({
    # Unsupervised UMAP
    seurat_obj <- RunModuleUMAP(
        seurat_obj,
        n_hubs = 10,
        n_neighbors = 15,
        min_dist = 0.1
    )

    p_umap <- ModuleUMAPPlot(
        seurat_obj,
        edge.alpha = 0.25,
        sample_edges = TRUE,
        edge_prop = 0.1,
        label_hubs = 2
    )
    ggsave("module_umap.pdf", p_umap, width = 10, height = 10)
    cat("Module UMAP saved.\n")
}, error = function(e) cat(sprintf("Module UMAP failed: %s\n", e$message)))

# Supervised UMAP
tryCatch({
    seurat_obj <- RunModuleUMAP(
        seurat_obj,
        n_hubs = 10,
        supervised = TRUE,
        target_weight = 0.5
    )

    p_umap_sup <- ModuleUMAPPlot(
        seurat_obj,
        edge.alpha = 0.25,
        sample_edges = TRUE,
        edge_prop = 0.1,
        label_hubs = 2
    )
    ggsave(file.path("network_plots", "module_umap_supervised.pdf"),
           p_umap_sup, width = 10, height = 10)
    cat("Supervised module UMAP saved.\n")
}, error = function(e) cat(sprintf("Supervised UMAP failed: %s\n", e$message)))

# --- Save enriched object ---------------------------------------------------
cat(sprintf("Saving enriched Seurat object...\n"))
saveRDS(seurat_obj, opts$seurat_rds)  # overwrite with enrichment results added

cat("hdWGCNA enrichment & network visualization complete.\n")
