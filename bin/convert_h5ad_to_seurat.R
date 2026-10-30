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
        } else if (args[i] == "--native_counts_mapping") {
            # PROBE FLAG, default OFF -- see the block at the as.Seurat call.
            arg_list$native_counts_mapping <- TRUE
            i <- i + 1
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

# Exact integral test, chunked with early exit: a normalised matrix trips on
# the first chunk, and a 909M-nnz raw matrix never materialises a full logical
# vector. Mirrors the helper in bin/run_cellchat_infer.R -- the two scripts are
# standalone by design and share no library.
is_integral_matrix <- function(m, chunk = 1e7) {
    x <- m@x
    n <- length(x)
    if (n == 0) return(TRUE)
    starts <- seq(1, n, by = chunk)
    for (s in starts) {
        e <- min(s + chunk - 1, n)
        if (any(x[s:e] != floor(x[s:e]))) return(FALSE)
    }
    TRUE
}

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
    #
    # The two producers of the h5ad we read disagree about what X means:
    #
    #   native RNA path (PLOT_POST_SCANVI, scanpy)
    #       X = log-normalised data, layers/counts = raw counts, no logcounts
    #   seurat_rds on-ramp (bin/convert_seurat_to_h5ad.R)
    #       X = raw counts, layers/logcounts = normalised data
    #
    # as.Seurat(counts = "X", data = NULL) leaves the `data` layer holding
    # whatever X held. On the native convention that is normalised values, so
    # it was right by accident. On the on-ramp convention it is RAW COUNTS,
    # and run_hdwgcna_celltype.R does SetDatExpr(..., layer = 'data').
    #
    # PRECISELY WHAT THAT BREAKS (verified in hdWGCNA source; commit 4020d42's
    # message overstates this and is corrected here). It does NOT put raw
    # counts into WGCNA's correlation matrix: MetacellsByGroups() defaults to
    # layer = "counts", and NormalizeMetacells() runs NormalizeData() on the
    # metacell object before SetDatExpr(use_metacells = TRUE) reads it, so the
    # matrix WGCNA correlates is normalised by hdWGCNA itself. What the
    # corrupted `data` layer corrupts is ScaleData -> RunPCA -> the kNN that
    # COMPOSES the metacells -- measurably: 991 vs 1005 metacells for L6 IT,
    # 440 vs 502 for Sst, which is what dragged the module calls down to
    # jaccard 0.41. A metacell-composition defect, not raw counts reaching the
    # correlation -- still wrong biology, silently, with no crash.
    #
    # The blunt "runs on raw counts" wording IS literally true of CellChat,
    # which has no metacell step: run_cellchat_infer.R hands the `data` layer
    # straight to CellChat. That path does not use this converter at all; it
    # was fixed separately in 4cf3975.
    #
    # Caught by the Tier 2 round-trip gate: the native and round-trip arms
    # differed in `data` by max 411157 (= 411165 raw max - 8.13 normalised
    # max) while `counts` was bit-identical, which is exactly why the Tier 1
    # matrix-fidelity gate could never have seen it.
    #
    # Resolution order, matching bin/run_cellchat_infer.R:
    #   1. logcounts, when the producer supplied it        (the on-ramp)
    #   2. X as-is, when X is already non-integral         (the native path)
    #   3. CPM + log1p, when X is integral and no logcounts exists
    #
    # Case 2 is the native path and is bit-for-bit what it has always been:
    # as.Seurat(data = NULL) leaves `data` holding X, which on that convention
    # is already log-normalised.
    #
    # Case 3 is an on-ramped object that was NEVER NORMALISED -- SEA-AD is
    # exactly this: its `data` slot is identical to its `counts` slot, so
    # convert_seurat_to_h5ad.R takes its "not writing logcounts" branch and no
    # logcounts layer reaches us. Without this branch `data` would hold raw
    # counts and ScaleData -> RunPCA -> the metacell kNN would run on them,
    # silently, exactly as the on-ramp convention bug did before 4020d42. The
    # CellChat readers gained this backstop in 4cf3975; this path did not, and
    # that asymmetry is what is being closed here.
    message("Converting SingleCellExperiment to Seurat...")
    sce_assays <- SummarizedExperiment::assayNames(sce)
    message("  SCE assays: ", paste(sce_assays, collapse = ", "))
    has_logcounts <- "logcounts" %in% sce_assays
    if (has_logcounts) {
        message("  logcounts present -> data layer (X stays as counts)")
    } else {
        x_mat <- SummarizedExperiment::assay(sce, "X")
        if (!inherits(x_mat, "dgCMatrix"))
            x_mat <- as(as(x_mat, "CsparseMatrix"), "dgCMatrix")
        if (is_integral_matrix(x_mat)) {
            message("  No logcounts layer and X is integral (raw counts) ",
                    "-> normalising with CPM + log1p for the data layer")
            library.size <- Matrix::colSums(x_mat)
            if (any(library.size == 0)) {
                stop("Cannot normalise: ", sum(library.size == 0),
                     " cells have zero total counts")
            }
            # Scale each column by its library size IN PLACE.
            #
            # The textbook spelling, log1p(t(t(m) / ls) * 1e4), materialises
            # two full transposes. On a 909M-non-zero object -- SEA-AD, which
            # is precisely the kind of never-normalised input that reaches
            # this branch -- that is tens of GB of temporaries and the process
            # is OOM-killed with no R traceback at all.
            #
            # A dgCMatrix is column-compressed, so every non-zero's column
            # (and therefore its library size) is recoverable from @p. That
            # turns the whole normalisation into one pass over @x. Structural
            # zeros are untouched, which is correct: log1p(0) = 0, so sparsity
            # is preserved exactly as the dense spelling preserved it.
            per_nz_libsize <- rep.int(library.size, diff(x_mat@p))
            x_mat@x <- log1p(x_mat@x / per_nz_libsize * 1e4)
            SummarizedExperiment::assay(sce, "logcounts") <- x_mat
            has_logcounts <- TRUE
            rm(x_mat, per_nz_libsize)
        } else {
            message("  No logcounts layer; X is already non-integral, so it ",
                    "supplies both counts and data (native-path convention)")
            rm(x_mat)
        }
    }
    # --- PROBE: the native convention's counts/data mapping -----------------
    # OFF BY DEFAULT. Shipped behaviour is unchanged unless the flag is passed.
    #
    # On the native convention X is log-normalised and layers/counts holds the
    # real raw counts. The shipped mapping is counts = X, data = X, so the
    # Seurat `counts` slot holds NORMALISED values and the raw counts sitting
    # in layers/counts go unused. FindVariableFeatures(selection.method="vst")
    # reads `counts`, so native hdWGCNA runs select variable features from
    # normalised data. Wrong in principle.
    #
    # Fixing it would change results on the established native path -- every
    # native hdWGCNA run FORGE has ever done -- so the decision was to QUANTIFY
    # FIRST. This flag exists to make that measurable: the probe harness runs
    # the converter twice on the same h5ad, once without the flag and once
    # with, and scores the two hdWGCNA runs against each other.
    #
    # If the measurement says fix it, this stops being a flag and becomes the
    # default for the native branch. If it says the effect is negligible, that
    # is worth knowing too and the flag documents why nothing changed.
    use_native_mapping <- isTRUE(opt$native_counts_mapping) &&
                          !has_logcounts && "counts" %in% sce_assays
    if (isTRUE(opt$native_counts_mapping)) {
        if (use_native_mapping) {
            message("  PROBE: native counts mapping ON -> counts='counts' ",
                    "(raw), data='X' (normalised)")
        } else {
            message("  PROBE: --native_counts_mapping requested but NOT ",
                    "applicable (has_logcounts=", has_logcounts,
                    ", counts assay present=", "counts" %in% sce_assays,
                    "); using the shipped mapping")
        }
    }
    seurat_obj <- if (use_native_mapping) {
        as.Seurat(sce, counts = "counts", data = "X")
    } else {
        as.Seurat(sce, counts = "X",
                  data = if (has_logcounts) "logcounts" else NULL)
    }

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
