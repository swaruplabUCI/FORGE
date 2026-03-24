#!/usr/bin/env Rscript
# hdWGCNA single cell-type analysis - RDS-based

suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(cowplot)
    library(patchwork)
    library(WGCNA)
    library(hdWGCNA)
    library(ggrepel)
    library(gtools)
    library(ggpubr)
    library(igraph)
    library(ggraph)
    library(tidygraph)
    library(png)
    library(pheatmap)
})

# ============================================================================
# SIMPLE ARGUMENT PARSER 
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    arg_list <- list()
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--seurat_rds" && i < length(args)) {
            arg_list$seurat_rds <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--cell_type" && i < length(args)) {
            arg_list$cell_type <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--cell_type_key" && i < length(args)) {
            arg_list$cell_type_key <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--output_dir" && i < length(args)) {
            arg_list$output_dir <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--metadata" && i < length(args)) {
            arg_list$metadata <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--threads" && i < length(args)) {
            arg_list$threads <- as.integer(args[i + 1])
            i <- i + 2
        } else {
            i <- i + 1
        }
    }
    return(arg_list)
}

opt <- parse_args(args)

# Set defaults
if (is.null(opt$cell_type_key)) opt$cell_type_key <- "scanvi_prediction"
if (is.null(opt$output_dir)) opt$output_dir <- "."
if (is.null(opt$threads)) opt$threads <- 32
if (is.null(opt$metadata)) opt$metadata <- NULL

# Validate required arguments
if (is.null(opt$seurat_rds) || !file.exists(opt$seurat_rds)) {
    stop("ERROR: Input Seurat RDS file not found or not specified")
}

if (is.null(opt$cell_type)) {
    stop("ERROR: Cell type not specified")
}

# ============================================================================
# SETUP
# ============================================================================

set.seed(42)
theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = opt$threads)
options(future.globals.maxSize = 100 * 1024^3)

output_dir <- file.path(opt$output_dir, "FiguresTables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("=== hdWGCNA Analysis ===")
message("Processing cell type: ", opt$cell_type)
message("Output directory: ", output_dir)
message("Using ", opt$threads, " threads")

# ============================================================================
# CELL TYPE ABBREVIATION MAPPING
# ============================================================================

get_cell_type_abbreviation <- function(cell_type) {
  cell_type_map <- c(
    "L2/3 IT" = "L23", "Oligodendrocyte" = "OLG", "Astrocyte" = "AST",
    "L5 IT" = "L5I", "L4 IT" = "L4I", "Microglia-PVM" = "MGL",
    "Pvalb" = "PVB", "OPC" = "OPC", "Vip" = "VIP",
    "L6 IT" = "L6I", "Sst" = "SST", "Lamp5" = "LM5",
    "Lamp5 Lhx6" = "LLH", "Sncg" = "SNC", "L6 CT" = "L6C",
    "L6b" = "L6B", "L6 IT Car3" = "LIC", "Chandelier" = "CHD",
    "Endothelial" = "END", "L5/6 NP" = "LNP", "VLMC" = "VLM",
    "Pax6" = "PAX", "L5 ET" = "L5E", "Sst Chodl" = "SCH",
    "Pericyte" = "PER"
  )
  
  if (cell_type %in% names(cell_type_map)) {
    return(cell_type_map[cell_type])
  } else {
    clean_name <- gsub("[^A-Za-z0-9]", "", cell_type)
    return(toupper(paste0(substr(clean_name, 1, 1), 
                         substr(clean_name, nchar(clean_name)-1, nchar(clean_name)))))
  }
}

cell_abbrev <- get_cell_type_abbreviation(opt$cell_type)
wgcna_name <- paste0("hdWGCNA_", gsub("\\s+|/", "_", opt$cell_type))

# ============================================================================
# LOAD AND SUBSET SEURAT OBJECT
# ============================================================================

message("Loading Seurat object from RDS...")
seurat_obj <- readRDS(opt$seurat_rds)
message("Loaded Seurat object with ", ncol(seurat_obj), " cells")

message("Subsetting to ", opt$cell_type, " cells...")
seurat_obj <- subset(seurat_obj, subset = !!sym(opt$cell_type_key) == opt$cell_type)
n_cells <- ncol(seurat_obj)
message("Subset contains ", n_cells, " cells")

# Check minimum cell threshold
if(n_cells < 100) {
    warning_msg <- paste("Insufficient cells (", n_cells, ") for analysis. Minimum is 100.")
    warning(warning_msg)
    
    skip_summary <- data.frame(
        Cell_Type = opt$cell_type,
        Total_Cells = n_cells,
        Status = "Skipped - Insufficient cells"
    )
    
    write.csv(skip_summary, 
              file = file.path(opt$output_dir, paste0(opt$cell_type, "_analysis_skipped.csv")),
              row.names = FALSE)
    quit(save = "no", status = 0)
}

# ============================================================================
# PREPROCESSING
# ============================================================================

print("Finding variable features...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 12000)

print("Running ScaleData...")
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))

# Run PCA and UMAP
print("Running PCA and UMAP...")
n_features <- length(VariableFeatures(seurat_obj))
max_pcs <- min(50, n_cells - 1, n_features - 1)
umap_dims <- min(30, max_pcs)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = max_pcs, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:umap_dims, verbose = FALSE)

# ============================================================================
# hdWGCNA SETUP
# ============================================================================

print("Setting up for WGCNA...")
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = wgcna_name
)

# ============================================================================
# CONSTRUCT METACELLS
# ============================================================================

print("Constructing metacells...")
# Use sample if available, otherwise use a dummy variable
if("sample" %in% colnames(seurat_obj@meta.data)) {
    group_by_vars <- c(opt$cell_type_key, "sample")
} else {
    seurat_obj$dummy_group <- "all"
    group_by_vars <- c(opt$cell_type_key, "dummy_group")
}

seurat_obj <- MetacellsByGroups(
    seurat_obj,
    group.by = group_by_vars,
    reduction = 'pca',
    k = 20,
    max_shared = 10,
    ident.group = opt$cell_type_key,
    min_cells = 50,
    wgcna_name = wgcna_name
)

# Normalize metacell expression
seurat_obj <- NormalizeMetacells(seurat_obj, wgcna_name = wgcna_name)

# Get metacell object and check count
metacell_obj <- GetMetacellObject(seurat_obj, wgcna_name = wgcna_name)
n_metacells <- ncol(metacell_obj)
print(paste("Number of metacells created:", n_metacells))

if(n_metacells < 20) {
    warning_msg <- paste("Too few metacells (", n_metacells, 
                        ") for robust WGCNA analysis. Minimum recommended is 20.")
    warning(warning_msg)
    
    skip_summary <- data.frame(
        Cell_Type = opt$cell_type,
        Total_Cells = n_cells,
        Metacells_Created = n_metacells,
        Status = "Skipped - Insufficient metacells",
        Recommendation = "Adjust k parameter or combine with similar cell types"
    )
    
    write.csv(skip_summary, 
              file = file.path(opt$output_dir, paste0(opt$cell_type, "_analysis_skipped.csv")),
              row.names = FALSE)
    quit(save = "no", status = 0)
}

# Run PCA on metacells
print("Running PCA on metacells...")
n_genes <- nrow(metacell_obj)
n_features_meta <- min(2000, n_genes - 1, n_metacells * 100)
n_pcs_meta <- min(50, n_metacells - 1, n_features_meta - 1)

metacell_obj <- FindVariableFeatures(metacell_obj, selection.method = "vst", nfeatures = n_features_meta)
metacell_obj <- ScaleData(metacell_obj)
metacell_obj <- RunPCA(metacell_obj, features = VariableFeatures(object = metacell_obj), npcs = n_pcs_meta)

seurat_obj <- SetMetacellObject(seurat_obj, metacell_obj, wgcna_name = wgcna_name)

# Run Harmony if sample info available and there are multiple samples
if("sample" %in% colnames(metacell_obj@meta.data)) {
    n_samples <- length(unique(metacell_obj@meta.data$sample))
    if(n_samples >= 2) {
        print(paste0("Running Harmony integration by sample (", n_samples, " samples)..."))
        seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars = 'sample', wgcna_name = wgcna_name)
    } else {
        print("Single sample detected, skipping Harmony integration.")
    }
}

# ============================================================================
# SETUP EXPRESSION MATRIX
# ============================================================================

print("Setting up expression matrix...")
seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = opt$cell_type,
    group.by = opt$cell_type_key,
    assay = 'RNA',
    slot = 'data',
    wgcna_name = wgcna_name
)

# ============================================================================
# SOFT POWER SELECTION
# ============================================================================

print("Testing soft power thresholds...")
seurat_obj <- TestSoftPowers(
    seurat_obj,
    powers = c(1:10, seq(12, 30, by = 2)),
    networkType = 'signed',
    wgcna_name = wgcna_name
)

# Plot and save soft power results
power_table <- GetPowerTable(seurat_obj, wgcna_name = wgcna_name)
print("Power table:")
print(power_table)
write.csv(power_table, 
          file = file.path(output_dir, paste0(opt$cell_type, "_power_table.csv")),
          row.names = FALSE)

plot_list <- PlotSoftPowers(seurat_obj, wgcna_name = wgcna_name)
soft_power_png <- file.path(output_dir, paste0(opt$cell_type, "_soft_powers.png"))
png(soft_power_png, width = 2000, height = 1500, res = 200)
print(wrap_plots(plot_list, ncol=2))
dev.off()

# Select soft power
good_powers <- power_table[power_table$SFT.R.sq > 0.8 & power_table$mean.k. < 100, ]

if(nrow(good_powers) == 0) {
    print("No power meets R² > 0.8 and mean.k < 100, relaxing to mean.k < 200")
    good_powers <- power_table[power_table$SFT.R.sq > 0.8 & power_table$mean.k. < 200, ]
}

if(nrow(good_powers) == 0) {
    print("Warning: Using power with highest R² regardless of connectivity")
    best_idx <- which.max(power_table$SFT.R.sq)
    selected_power <- power_table$Power[best_idx]
} else {
    selected_power <- min(good_powers$Power)
}

if(selected_power < 3) {
    print(paste("Power", selected_power, "is too low, setting to minimum of 3"))
    selected_power <- 3
}

print(paste("Selected soft power:", selected_power))

# ============================================================================
# CONSTRUCT NETWORK
# ============================================================================

print("Constructing co-expression network...")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
unique_tom_name <- paste0(opt$cell_type, "_TOM_", timestamp)

seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power = selected_power,
    setDatExpr = FALSE,
    tom_name = unique_tom_name,
    networkType = 'signed',
    wgcna_name = wgcna_name
)

# Plot dendrogram
pdf(file.path(output_dir, paste0(opt$cell_type, "_dendrogram.pdf")), width = 12, height = 6)
PlotDendrogram(seurat_obj, main = paste(opt$cell_type, 'hdWGCNA Dendrogram'), wgcna_name = wgcna_name)
dev.off()

# ============================================================================
# MODULE EIGENGENES & CONNECTIVITY
# ============================================================================

print("Computing module eigengenes...")
if("sample" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "sample", wgcna_name = wgcna_name)
} else {
    seurat_obj <- ModuleEigengenes(seurat_obj, wgcna_name = wgcna_name)
}

print("Computing module connectivity...")
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = opt$cell_type_key,
    group_name = opt$cell_type,
    wgcna_name = wgcna_name
)

# Rename modules
seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(cell_abbrev, "-M"),
    wgcna_name = wgcna_name
)

# Get modules
modules <- GetModules(seurat_obj, wgcna_name = wgcna_name)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']
print(paste("Number of modules:", length(mods)))

write.csv(modules, 
          file = file.path(output_dir, paste0(opt$cell_type, "_modules.csv")),
          row.names = FALSE)

# ============================================================================
# MODULE TRAIT CORRELATION (QC metrics only)
# ============================================================================

print("Computing module trait correlations...")
cur_traits <- c()

if("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    cur_traits <- c(cur_traits, "nCount_RNA")
}
if("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    cur_traits <- c(cur_traits, "nFeature_RNA")
}
if("pct_counts_mt" %in% colnames(seurat_obj@meta.data)) {
    cur_traits <- c(cur_traits, "pct_counts_mt")
}
if("pct_counts_ribo" %in% colnames(seurat_obj@meta.data)) {
    cur_traits <- c(cur_traits, "pct_counts_ribo")
}

if(length(cur_traits) > 0) {
    tryCatch({
        seurat_obj <- ModuleTraitCorrelation(
            seurat_obj,
            traits = cur_traits,
            group.by = opt$cell_type_key,
            wgcna_name = wgcna_name
        )
        
        p <- PlotModuleTraitCorrelation(
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
            combine = TRUE,
            wgcna_name = wgcna_name
        )
        
        pdf(file.path(output_dir, paste0(opt$cell_type, "_module_trait_correlation.pdf")), 
            width = max(10, length(cur_traits) * 1.5), 
            height = max(8, length(mods) * 0.5))
        print(p)
        dev.off()
    }, error = function(e) {
        warning(paste("Module trait correlation failed:", e$message))
    })
}

# ============================================================================
# HUB GENES
# ============================================================================

print("Identifying hub genes...")
hub_genes_list <- list()
for(mod in mods) {
    hub_genes <- GetHubGenes(seurat_obj, n_hubs = 25, mods = mod, wgcna_name = wgcna_name)
    hub_genes_list[[mod]] <- hub_genes
    
    write.csv(hub_genes, 
              file = file.path(output_dir, paste0(opt$cell_type, "_", mod, "_hub_genes.csv")),
              row.names = FALSE)
}

all_hub_genes <- do.call(rbind, lapply(names(hub_genes_list), function(x) {
    df <- hub_genes_list[[x]]
    df$module <- x
    return(df)
}))
write.csv(all_hub_genes, 
          file = file.path(output_dir, paste0(opt$cell_type, "_all_hub_genes.csv")),
          row.names = FALSE)

# ============================================================================
# NETWORK PLOTS
# ============================================================================

print("Creating module network plots...")
ModuleNetworkPlot(
    seurat_obj,
    outdir = file.path(output_dir, 'ModuleNetworks'),
    wgcna_name = wgcna_name
)

# Hub gene network for top modules (by connectivity)
print("Creating hub gene network plots...")
if(length(mods) > 0) {
    top_modules <- mods[1:min(4, length(mods))]
    
    pdf(file.path(output_dir, paste0(opt$cell_type, "_hub_gene_networks.pdf")), 
        width = 16, height = 12)
    HubGeneNetworkPlot(
        seurat_obj,
        n_hubs = 10,
        n_other = 5,
        edge_prop = 0.75,
        mods = top_modules,
        wgcna_name = wgcna_name
    )
    dev.off()
}

# ============================================================================
# MODULE SCORES & VISUALIZATION
# ============================================================================

print("Computing module expression scores...")
seurat_obj <- ModuleExprScore(
    seurat_obj,
    n_genes = 25,
    method = 'Seurat',
    wgcna_name = wgcna_name
)

# Plot module scores on UMAP
plot_list <- ModuleFeaturePlot(
    seurat_obj,
    features = 'hMEs',
    order = TRUE,
    ucell = FALSE,
    reduction = "umap",
    wgcna_name = wgcna_name
)

pdf(file.path(output_dir, paste0(opt$cell_type, "_module_scores_umap.pdf")), 
    width = 20, height = ceiling(length(mods)/4) * 5)
print(wrap_plots(plot_list, ncol = 4))
dev.off()

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("Creating summary report...")
summary_stats <- data.frame(
    Metric = c("Total cells", 
               "Number of modules", 
               "Number of metacells",
               "Soft power used",
               "Total hub genes identified"),
    Value = c(
        n_cells,
        length(mods),
        n_metacells,
        selected_power,
        nrow(all_hub_genes)
    )
)

write.csv(summary_stats, 
          file = file.path(output_dir, paste0(opt$cell_type, "_analysis_summary.csv")),
          row.names = FALSE)

# ============================================================================
# SAVE RESULTS
# ============================================================================

print("Saving results...")
saveRDS(seurat_obj, file = file.path(opt$output_dir, paste0(opt$cell_type, "_seurat_with_wgcna.rds")))

# Get module eigengenes
MEs <- GetMEs(seurat_obj, harmonized = TRUE, wgcna_name = wgcna_name)

results <- list(
    seurat_obj = seurat_obj,
    modules = modules,
    MEs = MEs,
    hub_genes = all_hub_genes,
    summary_stats = summary_stats,
    selected_power = selected_power,
    cell_type = opt$cell_type,
    output_dir = output_dir,
    n_cells = n_cells,
    n_modules = length(mods),
    n_metacells = n_metacells
)

saveRDS(results, file = file.path(opt$output_dir, paste0(opt$cell_type, "_complete_results.rds")))

print("============================================================")
print("ANALYSIS COMPLETE!")
print("============================================================")
print(paste("Cell type:", opt$cell_type))
print(paste("Total cells:", n_cells))
print(paste("Modules identified:", length(mods)))
print(paste("Hub genes:", nrow(all_hub_genes)))
print(paste("Output directory:", output_dir))
print("============================================================")
