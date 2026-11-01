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

# ------------------------------
# Resolve the NORMALISED expression matrix
# ------------------------------
# CellChat's L-R inference expects log-normalised values. Getting this wrong
# is silent: nothing crashes, the probabilities are just wrong.
#
# The two producers of the h5ad we read disagree about what X means:
#
#   native RNA path (PLOT_POST_SCANVI, scanpy)
#       X = log-normalised data, layers/counts = raw counts, no logcounts
#   seurat_rds on-ramp (bin/convert_seurat_to_h5ad.R)
#       X = RAW COUNTS, layers/logcounts = normalised data
#
# schard reads X and nothing else -- h5ad2list() has no layer support at all --
# and CreateAssayObject(counts = X) populates `data` as a copy of `counts`.
# So Layers() always reports both, `data` is always just X, and the old
# `else if ("counts" %in% ...)` fallback below was unreachable dead code.
# On the on-ramp convention that meant raw counts (max 411165 on the Allen
# reference) went straight into createCellChat(). Unlike the hdWGCNA path,
# there is no metacell step to launder them: run_cellchat_infer.R hands this
# matrix directly to the inference.
#
# Resolution order:
#   1. layers/logcounts, when the producer supplied it  (the on-ramp)
#   2. X as-is, when X is already non-integral          (the native path)
#   3. CPM + log1p, when X is integral and no logcounts exists
#      (an on-ramped object that was never normalised -- e.g. SEA-AD, whose
#      `data` is identical to `counts`, so convert_seurat_to_h5ad.R writes
#      no logcounts layer. Previously this ran raw, silently.)
#
# Case 2 is the native path and is bit-for-bit what it has always been.

# Exact integral test, chunked with early exit: a normalised matrix trips on
# the first chunk, and a 909M-nnz raw matrix never materialises a full
# logical vector.
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

h5_struct <- rhdf5::h5ls(opt$input)
has_logcounts <- any(h5_struct$group == "/layers" & h5_struct$name == "logcounts")

if (has_logcounts) {
    message("Found layers/logcounts -> using it as the normalised matrix")
    data.input <- schard::h5ad2Matrix(opt$input, "layers/logcounts")
    # h5ad2Matrix returns the matrix without dimnames; schard's Seurat object
    # carries the var/obs names already made Seurat-safe (underscores -> dashes),
    # and the layer shares X's gene and cell ordering.
    if (!identical(dim(data.input), dim(seurat_obj))) {
        stop("layers/logcounts is ", paste(dim(data.input), collapse = " x "),
             " but X is ", paste(dim(seurat_obj), collapse = " x "),
             " -- refusing to pair mismatched matrices")
    }
    dimnames(data.input) <- dimnames(seurat_obj)
} else {
    assay_obj <- seurat_obj@assays[[DefaultAssay(seurat_obj)]]
    available_layers <- Layers(assay_obj)
    if (!("data" %in% available_layers || "counts" %in% available_layers)) {
        stop("Cannot find appropriate expression matrix in h5ad")
    }
    data.input <- GetAssayData(
        seurat_obj,
        layer = if ("data" %in% available_layers) "data" else "counts"
    )
    if (is_integral_matrix(data.input)) {
        message("No logcounts layer and X is integral (raw counts) ",
                "-> normalising with CPM + log1p")
        library.size <- Matrix::colSums(data.input)
        if (any(library.size == 0)) {
            stop("Cannot normalise: ", sum(library.size == 0),
                 " cells have zero total counts")
        }
        # Scale each column by its library size IN PLACE.
        #
        # The textbook spelling, log1p(t(t(m) / ls) * 1e4), materialises two
        # full transposes of the matrix. On a 909M-non-zero object -- SEA-AD,
        # which is precisely the kind of never-normalised input that reaches
        # this branch -- that is tens of GB of temporaries and the process is
        # OOM-killed with no R traceback at all.
        #
        # A dgCMatrix is column-compressed, so every non-zero's column (and
        # therefore its library size) is recoverable from @p. That turns the
        # whole normalisation into one pass over @x with no copy of the matrix.
        # Structural zeros are untouched, which is correct: log1p(0) = 0, so
        # sparsity is preserved exactly as the dense spelling preserved it.
        if (!inherits(data.input, "dgCMatrix"))
            data.input <- as(as(data.input, "CsparseMatrix"), "dgCMatrix")
        per_nz_libsize <- rep.int(library.size, diff(data.input@p))
        data.input@x <- log1p(data.input@x / per_nz_libsize * 1e4)
    } else {
        message("X is already non-integral -> using it as the normalised ",
                "matrix (native-path convention)")
    }
}
message("  Expression matrix: ", paste(dim(data.input), collapse = " x "),
        "  max = ", signif(max(data.input), 6))

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
