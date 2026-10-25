#!/usr/bin/env Rscript
# Reference side of the .rds on-ramp fidelity check.
#
# Reads the ORIGINAL Seurat object and emits a set of plain-text artifacts
# describing it. summarize_h5ad.py emits the identically-named artifacts from
# the CONVERTED h5ad, and compare.py diffs the two. Everything crosses the
# boundary as CSV/JSON so the two sides never share a runtime -- they do not
# even share a container.
#
# Deliberately recomputes every quantity from the matrices rather than reading
# cached meta.data columns: the point is to catch a converter that moved the
# numbers, and a cached column would be copied through unchanged and hide it.

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(jsonlite)
})

parse_args <- function(args) {
    known <- c("--rds", "--outdir", "--assay", "--cell_type_key", "--markers")
    out <- list(); i <- 1
    while (i <= length(args)) {
        if (args[i] %in% known && i < length(args)) {
            out[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2
        } else i <- i + 1
    }
    out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$rds))    stop("--rds required")
if (is.null(opt$outdir)) stop("--outdir required")
if (is.null(opt$assay))  opt$assay <- "RNA"
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(opt$rds)

# Match the converter's handling of Seurat 3/4 objects, so the reference side
# is measuring the same object the converter saw.
obj_version <- tryCatch(as.character(obj@version), error = function(e) "0")
if (utils::compareVersion(obj_version, "5.0.0") < 0) {
    obj <- Seurat::UpdateSeuratObject(obj)
}

DefaultAssay(obj) <- opt$assay
a <- obj[[opt$assay]]

if (inherits(a, "Assay5")) {
    lyr <- SeuratObject::Layers(a)
    if (sum(grepl("^counts", lyr)) > 1 || sum(grepl("^data", lyr)) > 1) {
        obj <- SeuratObject::JoinLayers(obj, assay = opt$assay)
        a <- obj[[opt$assay]]
    }
}

counts <- SeuratObject::LayerData(a, layer = "counts")
if (!inherits(counts, "dgCMatrix")) counts <- as(as(counts, "CsparseMatrix"), "dgCMatrix")
lognorm <- tryCatch(SeuratObject::LayerData(a, layer = "data"), error = function(e) NULL)
if (!is.null(lognorm) && !inherits(lognorm, "dgCMatrix"))
    lognorm <- as(as(lognorm, "CsparseMatrix"), "dgCMatrix")

cells <- colnames(counts)
genes <- rownames(counts)
meta  <- obj@meta.data[cells, , drop = FALSE]

# Resolve the cell type column the same way the converter does, so both sides
# label by the same thing.
ct_key <- opt$cell_type_key
if (is.null(ct_key) || !nzchar(ct_key) || !(ct_key %in% colnames(meta))) {
    for (cand in c("cell_type", "scanvi_prediction", "celltypist_prediction",
                   "cell_type_prediction", "cell_type_marker", "seurat_annotations")) {
        if (cand %in% colnames(meta)) { ct_key <- cand; break }
    }
}
cell_type <- as.character(meta[[ct_key]])

# ---- per-cell QC, recomputed from the counts matrix ----
per_cell <- data.frame(
    cell_id      = cells,
    total_counts = Matrix::colSums(counts),
    n_genes      = Matrix::colSums(counts > 0),
    cell_type    = cell_type,
    stringsAsFactors = FALSE
)
write.csv(per_cell, file.path(opt$outdir, "per_cell.csv"), row.names = FALSE)

# ---- per-gene totals ----
per_gene <- data.frame(
    gene         = genes,
    total_counts = Matrix::rowSums(counts),
    n_cells      = Matrix::rowSums(counts > 0),
    stringsAsFactors = FALSE
)
write.csv(per_gene, file.path(opt$outdir, "per_gene.csv"), row.names = FALSE)

# ---- embeddings ----
reds <- names(obj@reductions)
emb_rows <- list()
for (rn in reds) {
    e <- SeuratObject::Embeddings(obj, reduction = rn)
    if (is.null(e) || nrow(e) == 0) next
    e <- e[cells, , drop = FALSE]
    key <- if (grepl("^X_", rn)) rn else paste0("X_", tolower(rn))
    k <- min(2L, ncol(e))
    emb_rows[[key]] <- data.frame(
        cell_id = cells, embedding = key,
        dim1 = as.numeric(e[, 1]),
        dim2 = if (k > 1) as.numeric(e[, 2]) else NA_real_,
        stringsAsFactors = FALSE)
}
emb_df <- if (length(emb_rows)) do.call(rbind, emb_rows) else
    data.frame(cell_id = character(0), embedding = character(0),
               dim1 = numeric(0), dim2 = numeric(0))
write.csv(emb_df, file.path(opt$outdir, "embeddings.csv"), row.names = FALSE)

# ---- dot-plot quantities ----
# Mean lognorm expression and percent-expressing per (cell_type, gene).
# Computed by hand on both sides rather than via DotPlot()/sc.pl.dotplot() so
# the comparison measures the conversion, not two plotting libraries'
# differing scaling conventions.
expr <- if (!is.null(lognorm) && prod(dim(lognorm)) > 0) lognorm else counts
if (!is.null(opt$markers) && nzchar(opt$markers)) {
    marker_genes <- strsplit(opt$markers, ",")[[1]]
    marker_genes <- trimws(marker_genes)
    marker_genes <- marker_genes[marker_genes %in% genes]
} else {
    # Top 20 by total expression -- deterministic and present on both sides.
    marker_genes <- genes[order(Matrix::rowSums(expr), decreasing = TRUE)][1:min(20L, length(genes))]
}

cts <- sort(unique(cell_type[!is.na(cell_type) & nzchar(cell_type)]))
rows <- list()
for (ct in cts) {
    idx <- which(cell_type == ct)
    sub <- expr[marker_genes, idx, drop = FALSE]
    rows[[ct]] <- data.frame(
        cell_type     = ct,
        gene          = marker_genes,
        mean_expr     = as.numeric(Matrix::rowMeans(sub)),
        pct_expressed = as.numeric(Matrix::rowSums(sub > 0)) / length(idx) * 100,
        stringsAsFactors = FALSE)
}
dot_df <- do.call(rbind, rows)
write.csv(dot_df, file.path(opt$outdir, "dotplot.csv"), row.names = FALSE)

# ---- scalar summary ----
num_cols <- names(meta)[vapply(meta, is.numeric, logical(1))]
summary_obj <- list(
    source          = "seurat_rds",
    n_cells         = ncol(counts),
    n_genes         = nrow(counts),
    x_sum           = sum(counts@x),
    x_nnz           = length(counts@x),
    has_lognorm     = !is.null(lognorm) && prod(dim(lognorm)) > 0,
    lognorm_sum     = if (!is.null(lognorm) && prod(dim(lognorm)) > 0) sum(lognorm@x) else NA_real_,
    cell_type_key   = ct_key,
    cell_type_counts = as.list(table(cell_type)),
    embeddings      = as.list(vapply(emb_rows, function(d) nrow(d), numeric(1))),
    obs_numeric_means = as.list(vapply(num_cols, function(cn) mean(meta[[cn]], na.rm = TRUE), numeric(1))),
    first_cells     = head(cells, 5),
    first_genes     = head(genes, 5),
    marker_genes    = marker_genes
)
write(toJSON(summary_obj, auto_unbox = TRUE, digits = 10, na = "null"),
      file.path(opt$outdir, "summary.json"))

cat("Wrote Seurat-side artifacts to ", opt$outdir, "\n", sep = "")
cat("  cells=", ncol(counts), " genes=", nrow(counts),
    " x_sum=", format(sum(counts@x), scientific = FALSE),
    " cell_type_key=", ct_key, "\n", sep = "")
