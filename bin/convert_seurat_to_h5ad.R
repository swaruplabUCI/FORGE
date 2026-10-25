#!/usr/bin/env Rscript
# Seurat .rds -> h5ad conversion for the params.onramp.seurat_rds on-ramp.
#
# Writes the AnnData HDF5 layout directly with rhdf5 rather than going through
# zellkonverter. That is deliberate. zellkonverter drives a Basilisk-managed
# Python/anndata, which in this pipeline has cost us a ~79 MB Miniforge
# download per task (see the containerOptions hack on CONVERT_H5AD_TO_SEURAT)
# and three separate encoding-type workarounds in convert_h5ad_to_seurat.R
# (FIX-65/73/80). rhdf5 is already in seurat_extended.sif, needs no Python, and
# lets us emit exactly the schema FORGE consumes.
#
# The output must satisfy the same contract as PLOT_POST_SCANVI's
# annotated_with_celltype.h5ad, because it is injected at the identical point
# in main.nf (ch_integrated_rna):
#   X                      raw counts, cells x genes
#   layers/logcounts       normalised data, when distinct from counts
#   obs$cell_type          canonical label column
#   obs$cell_type_source   provenance stamp
#   obs$sample             sample identifier
#   var_names              gene symbols
#   obsm$X_umap / X_pca    embeddings, when the object carries them
#
# A note on the transpose: Seurat holds genes x cells as a dgCMatrix (CSC).
# AnnData wants cells x genes as CSR. Those two are the same bytes --
# CSC-by-cell and CSR-by-cell share data/indices/indptr -- so X is written
# without materialising a transpose.

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(rhdf5)
})

message("=== Packages loaded ===")

# ============================================================================
# AnnData HDF5 primitives
# ============================================================================

ANNDATA_VERSION   <- "0.1.0"
DATAFRAME_VERSION <- "0.2.0"
ARRAY_VERSION     <- "0.2.0"

# Attribute writer. anndata compares attrs by value, so encoding-type and
# encoding-version must be scalar variable-length UTF-8 strings, not
# 1-element arrays -- h5py surfaces the latter as ndarray and the comparison
# silently fails.
set_attr <- function(fid, path, name, value, scalar = TRUE) {
    obj_type <- tryCatch({
        o <- H5Oopen(fid, path); t <- H5Iget_type(o); H5Oclose(o); t
    }, error = function(e) NA_character_)
    is_group <- identical(obj_type, "H5I_GROUP")
    h <- if (is_group) H5Gopen(fid, path) else H5Dopen(fid, path)
    on.exit(if (is_group) H5Gclose(h) else H5Dclose(h), add = TRUE)
    if (is.character(value)) {
        h5writeAttribute(value, h, name, asScalar = scalar,
                         variableLengthString = TRUE, encoding = "UTF-8")
    } else {
        h5writeAttribute(value, h, name, asScalar = scalar)
    }
    invisible(NULL)
}

stamp <- function(fid, path, etype, eversion) {
    set_attr(fid, path, "encoding-type", etype)
    set_attr(fid, path, "encoding-version", eversion)
}

# Variable-length UTF-8 string dataset. Fixed-length storage would come back
# through h5py as bytes, and every downstream `obs['cell_type'] == 'Astrocyte'`
# comparison would quietly evaluate False.
write_str <- function(file, path, v) {
    v <- as.character(v)
    v[is.na(v)] <- ""
    h5createDataset(file, path, dims = length(v), storage.mode = "character",
                    size = NULL, encoding = "UTF-8")
    h5write(v, file, path)
}

# 1-D numeric/integer dataset, chunked + deflated. The sparse X arrays are the
# bulk of the file, so compression here is what keeps an on-ramped
# SEA-AD-sized object from ballooning on disk.
write_vec <- function(file, path, v, level = 3L) {
    mode  <- if (is.integer(v)) "integer" else "double"
    n     <- length(v)
    chunk <- max(1L, min(as.integer(n), 262144L))
    h5createDataset(file, path, dims = n, storage.mode = mode,
                    chunk = chunk, level = level)
    h5write(v, file, path)
}

# Sparse CSR group. `m` is the Seurat-native genes x cells dgCMatrix; its CSC
# slots are reused verbatim as CSR for the cells x genes AnnData orientation.
write_csr_from_csc <- function(file, fid, path, m) {
    h5createGroup(file, path)
    write_vec(file, paste0(path, "/data"),    as.numeric(m@x))
    write_vec(file, paste0(path, "/indices"), as.integer(m@i))
    write_vec(file, paste0(path, "/indptr"),  as.integer(m@p))
    stamp(fid, path, "csr_matrix", ANNDATA_VERSION)
    # shape is (n_obs, n_vars) = (cells, genes) = reversed Seurat dim
    set_attr(fid, path, "shape", as.integer(c(ncol(m), nrow(m))), scalar = FALSE)
}

# Categorical column: categories + 0-based codes, -1 for NA. Used for every
# string column so that scanpy's plotting paths (which reach for `.cat`) work
# without a conversion step.
write_categorical <- function(file, fid, path, v) {
    v     <- as.character(v)
    cats  <- sort(unique(v[!is.na(v) & nzchar(v)]))
    codes <- match(v, cats) - 1L
    codes[is.na(codes)] <- -1L
    h5createGroup(file, path)
    write_str(file, paste0(path, "/categories"), cats)
    write_vec(file, paste0(path, "/codes"), as.integer(codes))
    stamp(fid, path, "categorical", DATAFRAME_VERSION)
    set_attr(fid, path, "ordered", FALSE)
    stamp(fid, paste0(path, "/categories"), "string-array", ARRAY_VERSION)
    stamp(fid, paste0(path, "/codes"), "array", ARRAY_VERSION)
}

# A dataframe group (obs or var). column-order must be non-empty -- a
# zero-length HDF5 attribute is not writable -- so callers guarantee at least
# one column.
write_dataframe <- function(file, fid, path, df, index) {
    h5createGroup(file, path)
    write_str(file, paste0(path, "/_index"), index)
    cols <- colnames(df)
    for (cn in cols) {
        col   <- df[[cn]]
        dpath <- paste0(path, "/", cn)
        if (is.character(col) || is.factor(col)) {
            write_categorical(file, fid, dpath, col)
        } else if (is.logical(col)) {
            # Written as int rather than an HDF5 enum; anndata reads it back as
            # an integer column, which every consumer here treats the same.
            write_vec(file, dpath, as.integer(col))
            stamp(fid, dpath, "array", ARRAY_VERSION)
        } else {
            write_vec(file, dpath, as.numeric(col))
            stamp(fid, dpath, "array", ARRAY_VERSION)
        }
    }
    stamp(fid, path, "dataframe", DATAFRAME_VERSION)
    set_attr(fid, path, "_index", "_index")
    set_attr(fid, path, "column-order", as.character(cols), scalar = FALSE)
    stamp(fid, paste0(path, "/_index"), "string-array", ARRAY_VERSION)
}

# ============================================================================
# Arguments
# ============================================================================

parse_args <- function(args) {
    known <- c("--rds", "--output", "--assay", "--cell_type_key",
               "--sample_key", "--celltypes_out")
    out <- list(); i <- 1
    while (i <= length(args)) {
        if (args[i] %in% known && i < length(args)) {
            out[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2
        } else i <- i + 1
    }
    out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(opt$rds) || !file.exists(opt$rds))
    stop("ERROR: input .rds not found or not specified: ", opt$rds)
if (is.null(opt$output))
    stop("ERROR: output path not specified (--output)")
if (is.null(opt$assay))      opt$assay      <- "RNA"
if (is.null(opt$sample_key)) opt$sample_key <- "orig.ident"
if (is.null(opt$celltypes_out))
    opt$celltypes_out <- sub("\\.h5ad$", "_celltypes.txt", opt$output)

message("=== Seurat -> h5ad (rhdf5 direct) ===")
message("Input:  ", opt$rds)
message("Output: ", opt$output)

# ============================================================================
# Load and validate
# ============================================================================

obj <- readRDS(opt$rds)

# A .rds can hold anything. Report the class we actually got rather than
# letting a downstream slot access produce an inscrutable S4 error.
if (!inherits(obj, "Seurat"))
    stop("ERROR: ", opt$rds, " does not contain a Seurat object (class: ",
         paste(class(obj), collapse = ", "), "). params.onramp.seurat_rds ",
         "expects a Seurat object holding the integrated RNA data.")

# Objects saved by Seurat 3/4 lack slots that SeuratObject 5 expects (the
# Allen mouse brain reference, for one, is 3.0.0.9000 and has no `images`
# slot). Accessing any of them throws a bare "no slot of name" error long
# before we reach anything on-ramp-specific, so normalise up front. On an
# already-current object UpdateSeuratObject is a no-op, but it is not free on
# a large one, so only pay for it when the stamped version calls for it.
obj_version <- tryCatch(as.character(obj@version), error = function(e) "0")
if (utils::compareVersion(obj_version, "5.0.0") < 0) {
    message("  Object was written by Seurat ", obj_version,
            "; running UpdateSeuratObject()...")
    obj <- Seurat::UpdateSeuratObject(obj)
    message("  Updated to ", as.character(obj@version))
}

message("Loaded: ", ncol(obj), " cells x ", nrow(obj), " features")
message("  Assays: ", paste(names(obj@assays), collapse = ", "))

if (!(opt$assay %in% names(obj@assays)))
    stop("ERROR: assay '", opt$assay, "' not found. Available: ",
         paste(names(obj@assays), collapse = ", "),
         ". Set params.onramp_seurat.assay to one of these.")

assay_obj <- obj[[opt$assay]]

# A Signac ChromatinAssay is peaks x cells. It would convert cleanly and then
# fail much later inside DORC/SCENIC+ on a gene-symbol mismatch, so reject it
# here where the message can say what to do instead.
if (inherits(assay_obj, "ChromatinAssay"))
    stop("ERROR: assay '", opt$assay, "' is a Signac ChromatinAssay (peaks x ",
         "cells). params.onramp.seurat_rds is the RNA on-ramp and needs a ",
         "gene x cell assay. For ATAC use params.onramp.atac_peak_matrix_h5ad.")

DefaultAssay(obj) <- opt$assay

# Seurat v5 splits counts across layers (counts.1, counts.2, ...) when an
# object was built from multiple samples. Reading one layer would silently
# drop most of the dataset, so join first.
if (inherits(assay_obj, "Assay5")) {
    lyr <- SeuratObject::Layers(assay_obj)
    message("  v5 layers: ", paste(lyr, collapse = ", "))
    if (sum(grepl("^counts", lyr)) > 1 || sum(grepl("^data", lyr)) > 1) {
        message("  Joining split v5 layers...")
        obj <- SeuratObject::JoinLayers(obj, assay = opt$assay)
        assay_obj <- obj[[opt$assay]]
        message("  After join: ",
                paste(SeuratObject::Layers(assay_obj), collapse = ", "))
    }
}

fetch_layer <- function(a, name) {
    tryCatch(SeuratObject::LayerData(a, layer = name), error = function(e) NULL)
}

counts <- fetch_layer(assay_obj, "counts")
dat    <- fetch_layer(assay_obj, "data")

if (is.null(counts) || prod(dim(counts)) == 0) {
    # hdWGCNA, DORC and SCENIC+ all model count data; silently promoting
    # normalised values into X would corrupt them without an error anywhere.
    if (!is.null(dat) && prod(dim(dat)) > 0)
        stop("ERROR: assay '", opt$assay, "' has a data layer but no counts. ",
             "The RNA on-ramp needs raw counts (DORC and SCENIC+ model count ",
             "data). Populate the counts layer, or on-ramp an h5ad via ",
             "params.onramp.rna_integrated_h5ad instead.")
    stop("ERROR: assay '", opt$assay, "' has no counts layer.")
}

if (!inherits(counts, "dgCMatrix"))
    counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")

has_logcounts <- !is.null(dat) && prod(dim(dat)) > 0
if (has_logcounts) {
    if (!inherits(dat, "dgCMatrix"))
        dat <- as(as(dat, "CsparseMatrix"), "dgCMatrix")
    # Seurat leaves data == counts until NormalizeData runs; no point storing
    # the same matrix twice.
    if (identical(length(dat@x), length(counts@x)) &&
        isTRUE(all.equal(dat@x, counts@x))) {
        message("  data layer is identical to counts; not writing logcounts")
        has_logcounts <- FALSE
    }
}

genes <- rownames(counts)
cells <- colnames(counts)
message("  Matrix: ", length(genes), " genes x ", length(cells), " cells, ",
        length(counts@x), " non-zero")

# ============================================================================
# Cell type resolution
# ============================================================================

meta <- obj@meta.data[cells, , drop = FALSE]

# Mirror plot_postscanvi.py's precedence so an object that already passed
# through FORGE round-trips onto the column it came from.
CT_CANDIDATES <- c("cell_type", "scanvi_prediction", "celltypist_prediction",
                   "cell_type_prediction", "cell_type_marker", "seurat_annotations")

has_real_labels <- function(col) {
    v <- unique(as.character(meta[[col]]))
    v <- v[!is.na(v) & nzchar(v)]
    length(v) > 1 || (length(v) == 1 && v[1] != "Unknown")
}

requested     <- opt$cell_type_key
cell_type_col <- NULL
if (!is.null(requested) && nzchar(requested) && requested != "null") {
    if (!(requested %in% colnames(meta)))
        stop("ERROR: cell_type_key '", requested, "' not found in meta.data. ",
             "Available columns: ", paste(colnames(meta), collapse = ", "))
    if (!has_real_labels(requested))
        stop("ERROR: cell_type_key '", requested, "' exists but holds no ",
             "usable labels (all NA/empty/'Unknown').")
    cell_type_col <- requested
} else {
    for (cand in CT_CANDIDATES)
        if (cand %in% colnames(meta) && has_real_labels(cand)) {
            cell_type_col <- cand; break
        }
    if (is.null(cell_type_col))
        stop("ERROR: no cell type column found. Tried: ",
             paste(CT_CANDIDATES, collapse = ", "),
             ". Set params.onramp_seurat.cell_type_key to the meta.data column ",
             "holding your labels. Available columns: ",
             paste(colnames(meta), collapse = ", "))
    message("  Auto-detected cell type column: ", cell_type_col)
}

meta$cell_type        <- as.character(meta[[cell_type_col]])
meta$cell_type_source <- "seurat_onramp"
message("  Canonical 'cell_type' stamped from '", cell_type_col, "'")

ctt <- sort(table(meta$cell_type), decreasing = TRUE)
message("  Cell types (n=", length(ctt), "): ",
        paste(utils::head(names(ctt), 5), collapse = ", "),
        if (length(ctt) > 5) ", ..." else "")

# 'sample' is read by the QC plots and by MULTIOME_INTEGRATION's sample join.
# A single-sample object is a legitimate on-ramp, so fall back rather than error.
if (!("sample" %in% colnames(meta))) {
    if (opt$sample_key %in% colnames(meta)) {
        meta$sample <- as.character(meta[[opt$sample_key]])
        message("  'sample' stamped from '", opt$sample_key, "'")
    } else {
        meta$sample <- "sample1"
        message("  No '", opt$sample_key, "' column; 'sample' set to 'sample1'")
    }
}

# anndata cannot hold list or matrix columns; flatten what we can, drop the rest.
drop_cols <- character(0)
for (cn in colnames(meta)) {
    col <- meta[[cn]]
    if (is.list(col) || !is.null(dim(col))) {
        drop_cols <- c(drop_cols, cn)
    } else if (is.factor(col)) {
        meta[[cn]] <- as.character(col)
    } else if (is.integer(col)) {
        meta[[cn]] <- as.numeric(col)
    } else if (!is.numeric(col) && !is.character(col) && !is.logical(col)) {
        meta[[cn]] <- as.character(col)
    }
}
if (length(drop_cols)) {
    message("  Dropping non-scalar meta.data columns: ",
            paste(drop_cols, collapse = ", "))
    meta <- meta[, setdiff(colnames(meta), drop_cols), drop = FALSE]
}

# ============================================================================
# Write
# ============================================================================

if (file.exists(opt$output)) invisible(file.remove(opt$output))
h5createFile(opt$output)
fid <- H5Fopen(opt$output)

# Every write below goes through this one handle. Handing the file *path* to
# h5write/h5createGroup while `fid` is open makes rhdf5 open a second handle to
# the same file; it warns about exactly that, and concurrent handles to one
# HDF5 file risk a corrupted superblock on writes the size of a real atlas.
f <- fid

stamp(fid, "/", "anndata", ANNDATA_VERSION)

message("Writing X (counts)...")
write_csr_from_csc(f, fid, "X", counts)

h5createGroup(f, "layers"); stamp(fid, "layers", "dict", ANNDATA_VERSION)
if (has_logcounts) {
    message("Writing layers/logcounts...")
    write_csr_from_csc(f, fid, "layers/logcounts", dat)
}

message("Writing obs (", ncol(meta), " columns)...")
write_dataframe(f, fid, "obs", meta, cells)

# var carries gene_name so that column-order is never empty and so gene
# symbols survive an index rename downstream.
message("Writing var...")
write_dataframe(f, fid, "var",
                data.frame(gene_name = genes, stringsAsFactors = FALSE), genes)

# Reductions -> obsm. FORGE's Python side looks for X_umap / X_pca, so
# normalise names (umap -> X_umap) on the way out.
h5createGroup(f, "obsm"); stamp(fid, "obsm", "dict", ANNDATA_VERSION)
reds <- names(obj@reductions)
if (length(reds)) {
    for (rn in reds) {
        emb <- SeuratObject::Embeddings(obj, reduction = rn)
        if (is.null(emb) || nrow(emb) == 0) next
        emb <- emb[cells, , drop = FALSE]
        key <- if (grepl("^X_", rn)) rn else paste0("X_", tolower(rn))
        storage.mode(emb) <- "double"
        dimnames(emb) <- NULL
        # Write the transpose, not the matrix. R is column-major and HDF5 is
        # row-major, so rhdf5 lands an R (cells, k) matrix as an HDF5 (k, cells)
        # dataset and anndata rejects the obsm shape. h5write's `native=TRUE`
        # would also fix it when given a file *path*, but is silently ignored
        # when given an open handle -- the flag baked into the handle by
        # H5Fopen wins. Transposing works either way.
        h5write(t(emb), f, paste0("obsm/", key))
        stamp(fid, paste0("obsm/", key), "array", ARRAY_VERSION)
        message("  obsm/", key, ": ", nrow(emb), " x ", ncol(emb))
    }
} else {
    message("  No dimensional reductions on the object")
}

for (g in c("varm", "obsp", "varp", "uns")) {
    h5createGroup(f, g); stamp(fid, g, "dict", ANNDATA_VERSION)
}

H5Fclose(fid)
h5closeAll()
message("Saved: ", opt$output)

# Sidecar cell-type list, mirroring convert_h5ad_to_seurat.R's contract.
cts <- unique(meta$cell_type)
cts <- cts[!is.na(cts) & nzchar(cts)]
writeLines(cts, opt$celltypes_out)
message("Exported ", length(cts), " cell types to: ", opt$celltypes_out)

message("=== Conversion Complete ===")
