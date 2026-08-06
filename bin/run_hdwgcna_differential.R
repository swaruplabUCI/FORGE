#!/usr/bin/env Rscript
# ===========================================================================
# hdWGCNA Differential Module Eigengene (DME) Analysis (FIX-29)
# Performs:
#   1. FindDMEs between control and treatment conditions
#   2. Module-trait correlation (if traits specified)
#   3. Lollipop and volcano plots for DMEs
#
# Expects an hdWGCNA-processed Seurat object (from HDWGCNA_PER_CELLTYPE)
# with completed network construction and module eigengenes.
#
# Usage:
#   run_hdwgcna_differential.R \
#     --seurat_rds hdwgcna_INH.rds \
#     --cell_type INH \
#     --cell_type_key scanvi_prediction \
#     --condition_key condition_group \
#     --control Control \
#     --treatment 90plus \
#     --traits "braaksc,pmi,msex,age_death,nCount_RNA,nFeature_RNA" \
#     --output_prefix hdwgcna_diff
# ===========================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(cowplot)
    library(patchwork)
    library(WGCNA)
    library(hdWGCNA)
    library(optparse)
})

theme_set(theme_cowplot())
set.seed(12345)

# --- CLI args ---------------------------------------------------------------
option_list <- list(
    make_option("--seurat_rds",    type = "character"),
    make_option("--cell_type",     type = "character"),
    make_option("--cell_type_key", type = "character", default = "scanvi_prediction"),
    make_option("--condition_key", type = "character", default = "condition_group"),
    make_option("--control",       type = "character", default = "Control"),
    make_option("--treatment",     type = "character", default = "90plus"),
    make_option("--traits",        type = "character", default = NULL),
    make_option("--group_mapping", type = "character", default = NULL,
                help = "JSON mapping sample -> condition_group (used if condition_key column is missing from the Seurat object)"),
    make_option("--output_prefix", type = "character", default = "hdwgcna_diff")
)
opts <- parse_args(OptionParser(option_list = option_list))

dir.create("dme_plots", recursive = TRUE, showWarnings = FALSE)

# --- Load -------------------------------------------------------------------
cat(sprintf("Loading Seurat object from %s\n", opts$seurat_rds))
seurat_obj <- readRDS(opts$seurat_rds)

# Verify hdWGCNA has been run
if (is.null(GetModules(seurat_obj))) {
    stop("No hdWGCNA modules found in Seurat object. Run HDWGCNA_PER_CELLTYPE first.")
}

# FIX-46 (ported from run_cellchat_per_condition.R): the hdWGCNA per-cell-type object
# is built from the annotated RNA object *before* condition_group is assigned, so the
# condition_key column is typically absent here. Recover it from the same
# sample -> condition_group JSON map the rest of the pipeline uses.
if (!opts$condition_key %in% colnames(seurat_obj@meta.data)) {
    if (is.null(opts$group_mapping) || !file.exists(opts$group_mapping)) {
        stop(sprintf("Condition column '%s' not found in Seurat meta.data and no valid --group_mapping provided.",
                     opts$condition_key))
    }
    cat(sprintf("Column '%s' not found — recovering from group_mapping JSON: %s\n",
                opts$condition_key, opts$group_mapping))
    mapping <- jsonlite::fromJSON(opts$group_mapping)   # named character vector: sample -> group
    sample_col <- if ("sample" %in% colnames(seurat_obj@meta.data)) "sample"
                  else if ("sample_id" %in% colnames(seurat_obj@meta.data)) "sample_id"
                  else if ("batch" %in% colnames(seurat_obj@meta.data)) "batch"
                  else stop("Cannot find a sample/sample_id/batch column to join the group mapping.")
    seurat_obj@meta.data[[opts$condition_key]] <-
        unname(mapping[as.character(seurat_obj@meta.data[[sample_col]])])
    assigned <- sum(!is.na(seurat_obj@meta.data[[opts$condition_key]]))
    cat(sprintf("  Mapped %d/%d cells to '%s' via '%s' column\n",
                assigned, nrow(seurat_obj@meta.data), opts$condition_key, sample_col))
    if (assigned == 0) {
        stop("group_mapping produced 0 assignments — sample names do not match the JSON keys.")
    }
}

# Sanitize cell type name for filesystem-safe filenames
safe_cell_type <- gsub("[^A-Za-z0-9._-]", "_", opts$cell_type)

cat(sprintf("Cell type: %s\n", opts$cell_type))
cat(sprintf("Condition key: %s\n", opts$condition_key))
cat(sprintf("Control: %s vs Treatment: %s\n", opts$control, opts$treatment))

# --- Step 1: FindDMEs between conditions ------------------------------------
cat("Running FindDMEs between conditions...\n")

# Get barcodes for each group (cells of this cell type in each condition)
group1_cells <- seurat_obj@meta.data %>%
    filter(!!sym(opts$cell_type_key) == opts$cell_type &
           !!sym(opts$condition_key) == opts$control) %>%
    rownames()

group2_cells <- seurat_obj@meta.data %>%
    filter(!!sym(opts$cell_type_key) == opts$cell_type &
           !!sym(opts$condition_key) == opts$treatment) %>%
    rownames()

cat(sprintf("  Control cells (%s): %d\n", opts$control, length(group1_cells)))
cat(sprintf("  Treatment cells (%s): %d\n", opts$treatment, length(group2_cells)))

if (length(group1_cells) < 10 || length(group2_cells) < 10) {
    cat("WARNING: Too few cells in one group (<10). Skipping DME analysis.\n")
    cat("DME analysis requires sufficient cells in both conditions.\n")
    writeLines("SKIPPED: insufficient cells", paste0(opts$output_prefix, "_SKIPPED.txt"))
    quit(save = "no", status = 0)
}

# Get wgcna_name (the active WGCNA experiment set during HDWGCNA_PER_CELLTYPE).
# NOTE: hdWGCNA exposes GetActiveWGCNAName() — there is no GetWGCNANames() (the prior
# call crashed every run once the condition_group fix let execution reach this point).
wgcna_name <- GetActiveWGCNAName(seurat_obj)
cat(sprintf("Active WGCNA name: %s\n", wgcna_name))

# Run DME test.
# barcodes1 = treatment (TG), barcodes2 = control (WT) so a POSITIVE avg_log2FC means the
# module eigengene is UP in the treatment condition — matching the MAST DEG "<treatment>_vs_
# <control>" convention and the hdWGCNA vignette's disease-vs-control orientation. (The prior
# order put control first, silently inverting every module's reported direction.)
DMEs <- FindDMEs(
    seurat_obj,
    barcodes1 = group2_cells,   # treatment (e.g. TG)
    barcodes2 = group1_cells,   # control   (e.g. WT)
    test.use = 'wilcox',
    wgcna_name = wgcna_name
)

# Save DME results
dme_file <- paste0("dme_results_", safe_cell_type, ".csv")
write.csv(DMEs, file = dme_file, row.names = TRUE)
cat(sprintf("DME results saved to %s\n", dme_file))

# Summarize significant DMEs
sig_dmes <- DMEs %>% filter(p_val_adj < 0.05)
cat(sprintf("Significant DMEs (padj < 0.05): %d / %d modules\n",
            nrow(sig_dmes), nrow(DMEs)))

# --- Step 2: DME Visualizations --------------------------------------------
cat("Generating DME plots...\n")

# Lollipop plot
tryCatch({
    p_lollipop <- PlotDMEsLollipop(
        seurat_obj, DMEs,
        wgcna_name = wgcna_name,
        pvalue = "p_val_adj"
    )
    ggsave(file.path("dme_plots", paste0("lollipop_", safe_cell_type, ".pdf")),
           p_lollipop, width = 8, height = 6)
}, error = function(e) cat(sprintf("Lollipop plot failed: %s\n", e$message)))

# Volcano plot
tryCatch({
    p_volcano <- PlotDMEsVolcano(
        seurat_obj, DMEs,
        wgcna_name = wgcna_name
    )
    ggsave(file.path("dme_plots", paste0("volcano_", safe_cell_type, ".pdf")),
           p_volcano, width = 8, height = 6)
}, error = function(e) cat(sprintf("Volcano plot failed: %s\n", e$message)))

# --- Step 3: Module-Trait Correlation (if traits provided) ------------------
if (!is.null(opts$traits) && opts$traits != "") {
    cat("Running module-trait correlation...\n")

    trait_list <- strsplit(opts$traits, ",")[[1]]
    trait_list <- trimws(trait_list)

    # Filter to traits that exist in metadata
    available_traits <- intersect(trait_list, colnames(seurat_obj@meta.data))
    missing_traits <- setdiff(trait_list, colnames(seurat_obj@meta.data))

    if (length(missing_traits) > 0) {
        cat(sprintf("  WARNING: Traits not found in metadata: %s\n",
                    paste(missing_traits, collapse = ", ")))
    }

    if (length(available_traits) > 0) {
        cat(sprintf("  Using traits: %s\n", paste(available_traits, collapse = ", ")))

        # Ensure proper data types for traits
        for (trait in available_traits) {
            vals <- seurat_obj@meta.data[[trait]]
            if (is.character(vals)) {
                unique_vals <- unique(vals[!is.na(vals)])
                if (length(unique_vals) == 2) {
                    cat(sprintf("  Converting binary trait '%s' to factor\n", trait))
                    seurat_obj@meta.data[[trait]] <- as.factor(seurat_obj@meta.data[[trait]])
                } else if (length(unique_vals) > 2) {
                    cat(sprintf("  WARNING: Skipping non-binary categorical trait '%s' (%d levels)\n",
                                trait, length(unique_vals)))
                    available_traits <- setdiff(available_traits, trait)
                }
            }
        }

        if (length(available_traits) > 0) {
            # Compute module-trait correlations
            tryCatch({
                seurat_obj <- ModuleTraitCorrelation(
                    seurat_obj,
                    traits = available_traits,
                    group.by = opts$cell_type_key
                )

                # Get and save correlation results
                mt_cor <- GetModuleTraitCorrelation(seurat_obj)

                # Plot heatmap
                p_trait <- PlotModuleTraitCorrelation(
                    seurat_obj,
                    label = 'fdr',
                    label_symbol = 'stars',
                    text_size = 2,
                    text_digits = 2,
                    text_color = 'white',
                    high_color = 'yellow',
                    mid_color = 'black',
                    low_color = 'purple',
                    plot_max = 0.2,
                    combine = TRUE
                )

                ggsave(paste0("module_trait_", safe_cell_type, ".pdf"),
                       p_trait, width = 12, height = 8)
                cat("Module-trait correlation heatmap saved.\n")
            }, error = function(e) {
                cat(sprintf("Module-trait correlation failed: %s\n", e$message))
            })
        }
    } else {
        cat("  No valid traits found in metadata — skipping module-trait correlation.\n")
    }
}

# --- Step 4: Enrichment on DME-significant modules --------------------------
cat("Running enrichment on significant DME modules...\n")

dir.create("enrichment_results", recursive = TRUE, showWarnings = FALSE)

tryCatch({
    # Run Enrichr on all modules
    dbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023')
    seurat_obj <- RunEnrichr(seurat_obj, dbs = dbs, max_genes = 100)
    enrich_df <- GetEnrichrTable(seurat_obj)

    write.csv(enrich_df, file.path("enrichment_results", paste0("enrichr_all_", safe_cell_type, ".csv")),
              row.names = FALSE)

    # Generate bar plots per module
    EnrichrBarPlot(
        seurat_obj,
        outdir = file.path("enrichment_results", "barplots"),
        n_terms = 10,
        plot_size = c(5, 7),
        logscale = TRUE
    )

    # Dot plot for biological process
    tryCatch({
        p_dot <- EnrichrDotPlot(
            seurat_obj,
            database = "GO_Biological_Process_2023",
            n_terms = 3
        )
        ggsave(file.path("enrichment_results", paste0("dotplot_BP_", safe_cell_type, ".pdf")),
               p_dot, width = 14, height = 10)
    }, error = function(e) cat(sprintf("Enrichr dot plot failed: %s\n", e$message)))

}, error = function(e) {
    cat(sprintf("Enrichment analysis failed: %s\n", e$message))
})

cat("hdWGCNA differential analysis complete.\n")
