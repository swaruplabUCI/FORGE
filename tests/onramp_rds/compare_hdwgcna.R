#!/usr/bin/env Rscript
# ============================================================================
# Round-trip hdWGCNA equivalence gate for the Seurat .rds on-ramp.
#
# Tier 1 (compare.py) proves the CONVERTED MATRIX is bit-exact. It cannot
# prove that an analysis run on top of that matrix lands in the same place --
# hdWGCNA soft-power selection and the WGCNA tree cut are both threshold
# decisions, so a 1e-15 drift could in principle tip a module split.
#
# This compares two hdWGCNA runs on the SAME cell type:
#   NATIVE    : run_hdwgcna_celltype.R straight off the original .rds
#   ROUNDTRIP : original .rds -> convert_seurat_to_h5ad.R (new, rhdf5)
#                            -> convert_h5ad_to_seurat.R (existing, zellkonverter)
#                            -> run_hdwgcna_celltype.R
#
# It therefore exercises BOTH converters and would catch a round-trip
# asymmetry that Tier 1 alone is blind to.
#
# Exit 0 = within tolerance. Exit 1 = a metric drifted past its threshold.
# ============================================================================
suppressPackageStartupMessages({
    library(Seurat); library(hdWGCNA)
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
    stop("usage: compare_hdwgcna.R --native <a.rds> --roundtrip <b.rds> [--label L] [--out tsv]")

cat("================================================================\n")
cat("hdWGCNA ROUND-TRIP EQUIVALENCE --", label, "\n")
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

# ---- 1. soft power ---------------------------------------------------------
# @misc[[1]] is "active_wgcna" -- an atomic character vector naming the active
# slot, NOT the slot itself. Indexing it with $ raises "$ operator is invalid
# for atomic vectors", which the old tryCatch swallowed into NA; NA vs NA then
# scored PASS. That silently neutralised the most diagnostic metric in the
# gate: job 58546541 reported soft_power NA/NA PASS while both arms had in
# fact selected power 9. Read it through the canonical accessor instead, and
# treat an unreadable power as a FAIL rather than a pass -- a metric we cannot
# evaluate is not a metric that agrees.
sp <- function(o) {
    p <- tryCatch(GetWGCNAParams(o)[["power"]], error = function(e) NULL)
    if (is.null(p) || length(p) != 1) return(NA_real_)
    p <- suppressWarnings(as.numeric(p))
    if (!is.finite(p)) return(NA_real_)
    p
}
spa <- sp(a); spb <- sp(b)
check("soft_power", spa, spb,
      is.finite(spa) && is.finite(spb) && isTRUE(all.equal(spa, spb)))

# ---- 2. modules ------------------------------------------------------------
ma <- GetModules(a); mb <- GetModules(b)
mods_a <- setdiff(unique(as.character(ma$module)), "grey")
mods_b <- setdiff(unique(as.character(mb$module)), "grey")
check("n_modules", length(mods_a), length(mods_b), length(mods_a) == length(mods_b))

# ---- 3. gene-level module assignment agreement -----------------------------
# Modules are colour-named per run; colours are not stable across runs. Match
# module A->B by best Jaccard, then score the fraction of genes that land in
# the matched module.
shared <- intersect(rownames(ma), rownames(mb))
check("genes_in_common", length(shared), length(shared), length(shared) > 0)

if (length(shared) > 0 && length(mods_a) > 0 && length(mods_b) > 0) {
    ga <- split(shared, as.character(ma[shared, "module"]))
    gb <- split(shared, as.character(mb[shared, "module"]))
    ga <- ga[setdiff(names(ga), "grey")]; gb <- gb[setdiff(names(gb), "grey")]

    jac <- function(x, y) length(intersect(x, y)) / length(union(x, y))
    best <- sapply(names(ga), function(m) {
        if (!length(gb)) return(0)
        max(sapply(gb, function(g) jac(ga[[m]], g)))
    })
    mean_j <- if (length(best)) mean(best) else 0
    check("mean_best_module_jaccard", sprintf("%.4f", mean_j),
          sprintf("%.4f", mean_j), mean_j >= 0.90)

    # weighted agreement: fraction of all module genes recovered in the match
    agree <- sum(sapply(names(ga), function(m) {
        if (!length(gb)) return(0)
        j <- sapply(gb, function(g) jac(ga[[m]], g))
        length(intersect(ga[[m]], gb[[names(which.max(j))]]))
    }))
    tot <- length(unlist(ga))
    frac <- if (tot) agree / tot else 0
    check("frac_genes_same_module", sprintf("%.4f", frac),
          sprintf("%.4f", frac), frac >= 0.90)
}

# ---- 4. kME correlation ----------------------------------------------------
kme_cols_a <- grep("^kME_", colnames(ma), value = TRUE)
kme_cols_b <- grep("^kME_", colnames(mb), value = TRUE)
if (length(kme_cols_a) && length(kme_cols_b) && length(shared)) {
    # compare the per-gene MAX kME, which is colour-label independent
    maxa <- apply(ma[shared, kme_cols_a, drop = FALSE], 1, max, na.rm = TRUE)
    maxb <- apply(mb[shared, kme_cols_b, drop = FALSE], 1, max, na.rm = TRUE)
    ok <- is.finite(maxa) & is.finite(maxb)
    r <- if (sum(ok) > 2) cor(maxa[ok], maxb[ok]) else NA_real_
    check("cor_max_kME", sprintf("%.6f", r), sprintf("%.6f", r),
          !is.na(r) && r >= 0.95)
}

# ---- 5. module eigengenes --------------------------------------------------
mea <- tryCatch(GetMEs(a, harmonized = FALSE), error = function(e) NULL)
meb <- tryCatch(GetMEs(b, harmonized = FALSE), error = function(e) NULL)
if (!is.null(mea) && !is.null(meb)) {
    cells <- intersect(rownames(mea), rownames(meb))
    check("cells_in_common", length(cells), length(cells), length(cells) > 0)
    if (length(cells) > 2) {
        # best |correlation| per native ME against any round-trip ME
        cm <- suppressWarnings(abs(cor(mea[cells, , drop = FALSE],
                                       meb[cells, , drop = FALSE])))
        cm[!is.finite(cm)] <- 0
        bestr <- apply(cm, 1, max)
        mb_r <- mean(bestr)
        check("mean_best_ME_cor", sprintf("%.6f", mb_r), sprintf("%.6f", mb_r),
              mb_r >= 0.95)
    }
}

# ---- 6. hub genes ----------------------------------------------------------
hub <- function(o, n = 25) {
    h <- tryCatch(GetHubGenes(o, n_hubs = n), error = function(e) NULL)
    if (is.null(h)) return(character(0))
    unique(as.character(h$gene_name))
}
ha <- hub(a); hb <- hub(b)
if (length(ha) && length(hb)) {
    ov <- length(intersect(ha, hb)) / length(union(ha, hb))
    check("hub_gene_jaccard", sprintf("%.4f", ov), sprintf("%.4f", ov),
          ov >= 0.80)
}

# ---- verdict ---------------------------------------------------------------
res <- do.call(rbind, rows)
if (!is.null(out_tsv)) {
    write.table(res, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\nWrote", out_tsv, "\n")
}
cat("================================================================\n")
if (length(fails)) {
    cat("RESULT:", label, "-- FAIL (", length(fails), "metric(s) ):",
        paste(fails, collapse = ", "), "\n")
    quit(status = 1)
}
cat("RESULT:", label, "-- PASS (", nrow(res), "metrics )\n")
quit(status = 0)
