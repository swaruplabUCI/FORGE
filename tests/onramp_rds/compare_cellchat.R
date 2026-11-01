#!/usr/bin/env Rscript
# ============================================================================
# CellChat round-trip equivalence gate for the Seurat .rds on-ramp.
#
# WHY THIS EXISTS, SEPARATELY FROM compare_hdwgcna.R
#
# hdWGCNA and CellChat were wired onto on-ramped RNA objects together
# (c61611b), but they do not share a reader. hdWGCNA goes through
# CONVERT_H5AD_TO_SEURAT (zellkonverter); CellChat reads the h5ad itself with
# schard, inside bin/run_cellchat_infer.R. So the logcounts fix in 4020d42
# covered hdWGCNA and left CellChat untouched.
#
# The failure it guards against was real and silent. The two producers of the
# integrated h5ad disagree about what X means:
#
#   native RNA path (PLOT_POST_SCANVI, scanpy)
#       X = log-normalised data, layers/counts = raw counts, no logcounts
#   seurat_rds on-ramp (bin/convert_seurat_to_h5ad.R)
#       X = RAW COUNTS, layers/logcounts = normalised data
#
# schard reads X and nothing else, and CreateAssayObject(counts = X) copies it
# into `data`. run_cellchat_infer.R took `data` as normalised, so on an
# on-ramped object CellChat received raw counts (max 411165 on the Allen
# reference). Unlike hdWGCNA there is no metacell step in between to launder
# them -- the matrix goes straight into createCellChat().
#
# THE TWO ARMS
#
# Both arms run the SAME production script, bin/run_cellchat_infer.R. They
# differ ONLY in which h5ad convention they are handed, which is exactly the
# variable under test:
#
#   NATIVE   : an h5ad written in the native convention (X = normalised,
#              no logcounts layer) -> the "X is already non-integral" branch
#   ROUNDTRIP: an h5ad from convert_seurat_to_h5ad.R (X = counts +
#              layers/logcounts) -> the "found layers/logcounts" branch
#
# Both are built from the SAME source object and the same normalised values,
# so a correct reader must land in the same place. That makes this gate do
# double duty: it proves the on-ramp branch is right AND that the native
# branch is unchanged, because a regression in either one breaks agreement.
#
# Exit 0 = within tolerance. Exit 1 = a metric drifted past its threshold.
#
# NOTE ON THRESHOLDS: the hdWGCNA determinism control (3ab607e) established a
# ZERO noise floor for that pipeline. CellChat's inference is also
# deterministic given a fixed matrix -- computeCommunProb()'s permutation test
# seeds off the data, not the clock -- so the honest expectation here is also
# exact agreement. The thresholds below are deliberately a little looser than
# that so a genuine numerical-precision effect is distinguishable from a
# control-flow bug, but anything short of ~1.0 is worth reading as signal, not
# as noise. Do NOT loosen them to make a run pass.
# ============================================================================
suppressPackageStartupMessages({
    library(CellChat)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) && i < length(args)) args[i + 1] else default
}
native_rds <- get_arg("--native")
rt_rds     <- get_arg("--roundtrip")
label      <- get_arg("--label", "unnamed")
out_tsv    <- get_arg("--out", NULL)

if (is.null(native_rds) || is.null(rt_rds))
    stop("usage: compare_cellchat.R --native <a.rds> --roundtrip <b.rds> [--label L] [--out tsv]")

cat("================================================================\n")
cat("CellChat ROUND-TRIP EQUIVALENCE --", label, "\n")
cat("================================================================\n")

a <- readRDS(native_rds)
b <- readRDS(rt_rds)

fails <- character(0)
rows  <- list()
record <- function(metric, native, roundtrip, verdict) {
    rows[[length(rows) + 1]] <<- data.frame(
        label = label, metric = metric,
        native = as.character(native), roundtrip = as.character(roundtrip),
        verdict = verdict, stringsAsFactors = FALSE)
    cat(sprintf("  %-28s native=%-18s rt=%-18s %s\n",
                metric, native, roundtrip, verdict))
}
check <- function(metric, nv, rv, ok) {
    v <- if (ok) "PASS" else "FAIL"
    if (!ok) fails <<- c(fails, metric)
    record(metric, nv, rv, v)
}

# ---- 1. cell groups --------------------------------------------------------
# If the two arms disagree about the grouping, nothing downstream is
# comparable and every later metric would be measuring the wrong thing.
ga <- levels(a@idents); gb <- levels(b@idents)
check("n_cell_groups", length(ga), length(gb), length(ga) == length(gb))
check("cell_groups_identical", paste(length(ga), "groups"),
      paste(length(gb), "groups"), identical(ga, gb))
check("n_cells", length(a@idents), length(b@idents),
      length(a@idents) == length(b@idents))

# ---- 2. over-expressed signalling genes ------------------------------------
# The first point where the expression VALUES (not just the gene set) steer
# the result: identifyOverExpressedGenes() does a per-group one-vs-rest test.
# Feeding it raw counts instead of normalised data changes this set, and
# everything downstream inherits that.
#
# identifyOverExpressedGenes() stores the significant features as a plain
# character vector under var.features[[features.name]] (default "features"),
# and the accompanying data.frame under "<name>.info". Pull the vector.
#
# Guard against the empty case explicitly: `$features` on a list with no such
# element returns NULL, as.character(NULL) is character(0) WITHOUT raising, so
# a tryCatch alone would let both arms resolve to empty and score jaccard = 1.
# That is the same shape of inert-metric bug as the soft_power NA/NA PASS that
# eae3fff fixed in the hdWGCNA gate -- an unreadable metric must FAIL, not pass.
pull_features <- function(o) {
    vf <- o@var.features
    v <- tryCatch(vf[["features"]], error = function(e) NULL)
    if (is.null(v) || length(v) == 0) {
        plain <- names(vf)[!grepl("\\.info$", names(vf))]
        if (length(plain)) v <- tryCatch(vf[[plain[1]]], error = function(e) NULL)
    }
    v <- as.character(v)
    v[!is.na(v)]
}
oga <- pull_features(a); ogb <- pull_features(b)
if (length(oga) == 0 || length(ogb) == 0) {
    check("overexpressed_genes_readable",
          length(oga), length(ogb), FALSE)
}
jac <- function(x, y) {
    x <- unique(x); y <- unique(y)
    if (length(x) == 0 && length(y) == 0) return(1)
    length(intersect(x, y)) / length(union(x, y))
}
j_og <- jac(oga, ogb)
check("n_overexpressed_genes", length(unique(oga)), length(unique(ogb)),
      length(unique(oga)) == length(unique(ogb)))
check("overexpressed_gene_jaccard", sprintf("%.4f", j_og), sprintf("%.4f", j_og),
      j_og >= 0.95)

# ---- 3. inferred communications --------------------------------------------
# subsetCommunication() is what the pipeline actually exports as the CSV, so
# this is the user-visible result.
na <- subsetCommunication(a)
nb <- subsetCommunication(b)
check("n_interactions", nrow(na), nrow(nb), nrow(na) == nrow(nb))

key <- function(df) paste(df$source, df$target, df$interaction_name, sep = "|")
ka <- key(na); kb <- key(nb)
j_int <- jac(ka, kb)
check("interaction_jaccard", sprintf("%.4f", j_int), sprintf("%.4f", j_int),
      j_int >= 0.95)

# ---- 4. communication probabilities on shared interactions -----------------
# A set-level match can still hide shifted probabilities, which is what
# determines the strength of every reported interaction.
shared_k <- intersect(ka, kb)
check("shared_interactions", length(shared_k), length(shared_k),
      length(shared_k) > 0)
if (length(shared_k) > 1) {
    pa <- na$prob[match(shared_k, ka)]
    pb <- nb$prob[match(shared_k, kb)]
    cor_p <- suppressWarnings(cor(pa, pb))
    if (!is.finite(cor_p)) cor_p <- if (isTRUE(all.equal(pa, pb))) 1 else NA_real_
    check("cor_prob", sprintf("%.6f", cor_p), sprintf("%.6f", cor_p),
          is.finite(cor_p) && cor_p >= 0.99)
    denom <- max(abs(c(pa, pb)), na.rm = TRUE)
    rel <- if (denom > 0) max(abs(pa - pb)) / denom else 0
    check("max_rel_prob_diff", sprintf("%.3e", rel), sprintf("%.3e", rel),
          is.finite(rel) && rel <= 0.01)
} else {
    check("cor_prob", "too few shared", "too few shared", FALSE)
}

# ---- 5. pathways -----------------------------------------------------------
pa_set <- unique(as.character(na$pathway_name))
pb_set <- unique(as.character(nb$pathway_name))
check("n_pathways", length(pa_set), length(pb_set),
      length(pa_set) == length(pb_set))
j_pw <- jac(pa_set, pb_set)
check("pathway_jaccard", sprintf("%.4f", j_pw), sprintf("%.4f", j_pw),
      j_pw >= 0.95)

# ---- 6. aggregated network -------------------------------------------------
# aggregateNet() is what the downstream plots and the centrality analysis
# consume, so an agreement here is the closest thing to "the figures match".
agg_ok <- TRUE
ca <- tryCatch(a@net$count,  error = function(e) NULL)
cb <- tryCatch(b@net$count,  error = function(e) NULL)
wa <- tryCatch(a@net$weight, error = function(e) NULL)
wb <- tryCatch(b@net$weight, error = function(e) NULL)
if (is.null(ca) || is.null(cb) || !identical(dim(ca), dim(cb))) {
    check("net_count_identical", "unreadable", "unreadable", FALSE)
    agg_ok <- FALSE
} else {
    check("net_count_identical", sum(ca), sum(cb), identical(ca, cb))
}
if (agg_ok && !is.null(wa) && !is.null(wb) && identical(dim(wa), dim(wb))) {
    cw <- suppressWarnings(cor(as.vector(wa), as.vector(wb)))
    if (!is.finite(cw)) cw <- if (isTRUE(all.equal(as.vector(wa), as.vector(wb)))) 1 else NA_real_
    check("cor_net_weight", sprintf("%.6f", cw), sprintf("%.6f", cw),
          is.finite(cw) && cw >= 0.99)
} else {
    check("cor_net_weight", "unreadable", "unreadable", FALSE)
}

# ---- verdict ---------------------------------------------------------------
res <- do.call(rbind, rows)
if (!is.null(out_tsv)) {
    write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\nwrote", out_tsv, "\n")
}
cat("----------------------------------------------------------------\n")
if (length(fails)) {
    cat("RESULT:", label, "FAIL (", length(fails), "of", nrow(res), ") --",
        paste(fails, collapse = ", "), "\n")
    quit(status = 1)
}
cat("RESULT:", label, "PASS (", nrow(res), "of", nrow(res), ")\n")
quit(status = 0)
