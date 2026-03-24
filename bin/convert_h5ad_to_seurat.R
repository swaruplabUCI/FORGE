#!/usr/bin/env Rscript
# FIX-54b: h5ad to Seurat conversion using zellkonverter (Basilisk-managed Python)
#
# Previous version required CONDA_PREFIX to configure reticulate for schard,
# which fails under Singularity --contain where no conda env is active.
# zellkonverter uses Basilisk to manage its own Python environment internally,
# so no external CONDA_PREFIX or reticulate config is needed.

suppressPackageStartupMessages({
    library(zellkonverter)
    library(SingleCellExperiment)
    library(Seurat)
    library(SeuratObject)
    library(rhdf5)
})

message("=== All packages loaded successfully ===")

# ============================================================================
# Parse command-line arguments
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
    arg_list <- list()
    i <- 1
    while (i <= length(args)) {
        if (args[i] == "--h5ad" && i < length(args)) {
            arg_list$h5ad <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--output" && i < length(args)) {
            arg_list$output <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--cell_type_key" && i < length(args)) {
            arg_list$cell_type_key <- args[i + 1]
            i <- i + 2
        } else if (args[i] == "--assay_name" && i < length(args)) {
            arg_list$assay_name <- args[i + 1]
            i <- i + 2
        } else {
            i <- i + 1
        }
    }
    return(arg_list)
}

opt <- parse_args(args)

# Validate inputs
if (is.null(opt$h5ad) || !file.exists(opt$h5ad)) {
    stop("ERROR: Input h5ad file not found or not specified: ", opt$h5ad)
}

if (is.null(opt$output)) {
    stop("ERROR: Output path not specified")
}

if (is.null(opt$cell_type_key)) {
    opt$cell_type_key <- "scanvi_prediction"
    message("Using default cell_type_key: scanvi_prediction")
}

if (is.null(opt$assay_name)) {
    opt$assay_name <- "RNA"
}

# ============================================================================
# Convert h5ad -> SingleCellExperiment -> Seurat
# ============================================================================

message("=== h5ad -> Seurat Conversion (zellkonverter) ===")
message("Input: ", opt$h5ad)
message("Output: ", opt$output)

tryCatch({
    # FIX-65: Strip uns/log1p if present — encoding_type='null' v0.1.0 is
    # unreadable by zellkonverter's basilisk-managed anndata
    message("Pre-processing: stripping problematic uns/log1p key if present...")
    tmp_h5ad <- file.path(tempdir(), basename(opt$h5ad))
    file.copy(opt$h5ad, tmp_h5ad, overwrite = TRUE)
    # FIX-73: Use rhdf5 (available in container) instead of python3+h5py (not available)
    # FIX-80: Recursively strip ALL datasets with encoding_type='null' (not just uns/log1p)
    if (rhdf5::H5Fis_hdf5(tmp_h5ad)) {
        fid <- rhdf5::H5Fopen(tmp_h5ad)
        tryCatch({
            # Walk the HDF5 tree and collect paths of datasets with encoding_type='null'
            all_items <- rhdf5::h5ls(fid, recursive = TRUE)
            to_delete <- c()
            for (i in seq_len(nrow(all_items))) {
                item_path <- paste0(all_items$group[i], "/", all_items$name[i])
                item_path <- sub("^//", "/", item_path)
                tryCatch({
                    attrs <- rhdf5::h5readAttributes(fid, item_path)
                    if (!is.null(attrs[["encoding-type"]]) && attrs[["encoding-type"]] == "null") {
                        to_delete <- c(to_delete, item_path)
                    }
                }, error = function(e) {})
            }
            # Also always strip uns/log1p if present (may not have encoding-type attr)
            if (rhdf5::H5Lexists(fid, "uns/log1p") && !("uns/log1p" %in% to_delete)) {
                to_delete <- c(to_delete, "uns/log1p")
            }
            if (length(to_delete) > 0) {
                for (path in to_delete) {
                    # Delete parent group if it's a single dataset within a group
                    rhdf5::h5delete(fid, path)
                    message("  Stripped: ", path)
                }
            } else {
                message("  No problematic encoding_type='null' keys found")
            }
        }, finally = rhdf5::H5Fclose(fid))
    }
    strip_result <- "OK"
    message("  h5py strip result: ", paste(strip_result, collapse = " "))
    opt$h5ad <- tmp_h5ad

    # Step 1: Read h5ad via zellkonverter (Basilisk handles Python internally)
    message("Reading h5ad via zellkonverter...")
    sce <- readH5AD(opt$h5ad)
    message("  SCE dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")

    # Step 2: Convert SCE to Seurat
    message("Converting SingleCellExperiment to Seurat...")
    seurat_obj <- as.Seurat(sce, counts = "X", data = NULL)

    message("  Cells: ", ncol(seurat_obj))
    message("  Genes: ", nrow(seurat_obj))
    message("  Assays (raw): ", paste(names(seurat_obj@assays), collapse = ", "))

    # FIX-89: Rename assay to match expected name (e.g. 'originalexp' -> 'RNA')
    current_assay <- names(seurat_obj@assays)[1]
    if (current_assay != opt$assay_name) {
        message("  FIX-89: Renaming assay '", current_assay, "' -> '", opt$assay_name, "'")
        rename_args <- list(seurat_obj)
        rename_args[[current_assay]] <- opt$assay_name
        seurat_obj <- do.call(RenameAssays, rename_args)
        message("  Assays (renamed): ", paste(names(seurat_obj@assays), collapse = ", "))
    }

    # Verify cell type key exists
    if (opt$cell_type_key %in% colnames(seurat_obj@meta.data)) {
        message("  Cell type key '", opt$cell_type_key, "' found in metadata")
        cell_type_counts <- table(seurat_obj@meta.data[[opt$cell_type_key]])
        message("  Cell types (n=", length(cell_type_counts), "): ",
                paste(names(head(cell_type_counts, 5)), collapse = ", "),
                ifelse(length(cell_type_counts) > 5, "...", ""))
    } else {
        warning("Cell type key '", opt$cell_type_key, "' not found!")
        message("  Available metadata columns: ",
                paste(head(colnames(seurat_obj@meta.data), 10), collapse = ", "))
    }

    # Save
    message("Saving RDS...")
    saveRDS(seurat_obj, file = opt$output)
    message("Saved to: ", opt$output)

    # Export unique cell types for nextflow
    if (opt$cell_type_key %in% colnames(seurat_obj@meta.data)) {
        cell_types <- unique(as.character(seurat_obj@meta.data[[opt$cell_type_key]]))
        cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]
        celltypes_file <- gsub('\\.rds$', '_celltypes.txt', opt$output)
        writeLines(cell_types, celltypes_file)
        message("Exported ", length(cell_types), " cell types to: ", celltypes_file)
    } else {
        warning("Cannot export cell types - key not found in metadata")
    }

}, error = function(e) {
    stop("FAILED to convert h5ad: ", e$message)
})

message("=== Conversion Complete ===")
