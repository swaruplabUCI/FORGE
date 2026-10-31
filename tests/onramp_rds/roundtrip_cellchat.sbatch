#!/bin/bash
#SBATCH --job-name=forge_rt_cellchat
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL
# ============================================================================
# CellChat round-trip equivalence for the Seurat .rds on-ramp.
#
# The companion to roundtrip_hdwgcna.sbatch, and NOT redundant with it:
# hdWGCNA and CellChat were wired onto on-ramps together (c61611b) but they do
# not share a reader. hdWGCNA goes through CONVERT_H5AD_TO_SEURAT
# (zellkonverter); CellChat reads the h5ad itself with schard, inside
# bin/run_cellchat_infer.R. The logcounts fix in 4020d42 therefore covered
# hdWGCNA and left CellChat reading raw counts.
#
# THE TRICK THAT MAKES THIS A FAIR TEST
#
# Both arms run the SAME production script on h5ads built by the SAME
# production converter from the SAME prepped object. The only difference is
# which slot the normalised values live in:
#
#   ONRAMP arm : convert_seurat_to_h5ad.R on the prepped object
#                -> X = raw counts, layers/logcounts = normalised
#   NATIVE arm : same converter, on a copy whose `counts` slot has been
#                replaced by its `data` slot. The converter skips logcounts
#                when data == counts (bin/convert_seurat_to_h5ad.R:270), so
#                this emits X = normalised, no logcounts -- which is exactly
#                the native scanpy convention that PLOT_POST_SCANVI produces.
#
# Both h5ads therefore carry NUMERICALLY IDENTICAL normalised values, plus
# identical obs/var/obsm, and a correct reader must land in the same place.
# That makes this gate do double duty: it proves the on-ramp branch is right
# AND that the native branch is unchanged, since a regression in either one
# breaks agreement.
#
# Usage:
#   sbatch --export=ALL,RDS=/path/x.rds,CT_KEY=subclass,SPECIES=mouse,\
#   N_CT=6,OUT=/path/out tests/onramp_rds/roundtrip_cellchat.sbatch
# ============================================================================
set -uo pipefail

REPO=${REPO:-/dfs7/swaruplab/lesolano/src_FORGE}
SING=${SING:-/opt/apps/singularity/3.11.3/bin/singularity}
SIF=$REPO/singularity_cache/seurat_extended.sif
RDS=${RDS:?set RDS=/path/to/input.rds}
CT_KEY=${CT_KEY:?set CT_KEY=<meta.data column>}
SPECIES=${SPECIES:?set SPECIES=human|mouse}
# CellChat runs across ALL cell groups at once rather than per cell type, so
# N_CT caps the number of groups instead of looping. 0 = use every group.
N_CT=${N_CT:-6}
SAMPLE_KEY=${SAMPLE_KEY:-orig.ident}
OUT=${OUT:?set OUT=/path/to/outdir}
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT"/{native,onramp,logs}
export HOME=$(mktemp -d /tmp/rt_cc_home_XXXXXX)
R="$SING exec -B /dfs7 -B /tmp $SIF Rscript"

echo "=============================================================="
echo "CellChat ROUND-TRIP EQUIVALENCE"
echo "  rds      : $RDS"
echo "  ct_key   : $CT_KEY"
echo "  species  : $SPECIES"
echo "  n_ct     : $N_CT  (0 = all groups)"
echo "  out      : $OUT"
echo "  threads  : $THREADS"
echo "  started  : $(date)"
echo "=============================================================="

# --- 0. prep the shared source object -------------------------------------
# Mirrors roundtrip_hdwgcna.sbatch step 0b: migrate pre-v5 objects (the Allen
# reference is 3.0.0.9000 and even basic accessors throw without this) and
# stamp `sample` so both arms agree on what a sample is.
#
# The normalisation ASSERTION below is the part specific to this gate. If the
# source object's `data` is identical to its `counts` -- i.e. it was never
# normalised, which is exactly the SEA-AD case -- then convert_seurat_to_h5ad.R
# writes no logcounts layer, BOTH arms come out byte-identical, and the gate
# passes while testing nothing at all. That is a false green, so refuse it
# loudly rather than report a meaningless PASS.
if [[ ! -s "$OUT/prepped.rds" ]]; then
    $R -e "
    suppressPackageStartupMessages({library(Seurat); library(Matrix)})
    o <- readRDS('$RDS')
    v <- package_version(as.character(o@version))
    cat('Input Seurat version:', as.character(v), '\n')
    if (v < '5.0.0') { cat('Pre-v5 object -- running UpdateSeuratObject()\n')
                       o <- UpdateSeuratObject(o) }
    stopifnot('$CT_KEY' %in% colnames(o@meta.data))
    md <- o@meta.data
    if (!('sample' %in% colnames(md))) {
        if ('$SAMPLE_KEY' %in% colnames(md)) {
            md[['sample']] <- as.character(md[['$SAMPLE_KEY']])
        } else md[['sample']] <- 'sample1'
        o@meta.data <- md
    }
    # keep only the N most abundant groups, so CellChat is not asked to infer
    # across groups too small for filterCommunication(min.cells = 10)
    n_ct <- $N_CT
    if (n_ct > 0) {
        t <- sort(table(as.character(o@meta.data[['$CT_KEY']])), decreasing = TRUE)
        t <- t[names(t) != 'NA' & nzchar(names(t))]
        keep <- head(names(t), n_ct)
        cat('Keeping', length(keep), 'groups:', paste(keep, collapse=' | '), '\n')
        cat('Counts :', paste(head(as.integer(t), n_ct), collapse=' | '), '\n')
        o <- o[, as.character(o@meta.data[['$CT_KEY']]) %in% keep]
    }
    cat('Prepped object:', nrow(o), 'genes x', ncol(o), 'cells\n')

    a <- o[[DefaultAssay(o)]]
    cnt <- SeuratObject::LayerData(a, layer='counts')
    dat <- tryCatch(SeuratObject::LayerData(a, layer='data'), error=function(e) NULL)
    if (is.null(dat) || prod(dim(dat)) == 0)
        stop('FATAL: object has no data layer -- nothing to compare against')
    same <- identical(length(dat@x), length(cnt@x)) && isTRUE(all.equal(dat@x, cnt@x))
    cat('counts max:', max(cnt), ' data max:', max(dat), ' identical:', same, '\n')
    if (same)
        stop('FATAL: this object was never normalised (data == counts). Both ',
             'arms would be byte-identical and the gate would pass without ',
             'testing anything. Use a normalised object.')
    cat('Normalisation assertion OK -- data is distinct from counts\n')
    saveRDS(o, '$OUT/prepped.rds')
    " 2>&1 | tee "$OUT/logs/00_prep.log"
fi
[[ -s "$OUT/prepped.rds" ]] || { echo "FATAL: prep failed"; exit 1; }

# --- 1. ONRAMP-convention h5ad --------------------------------------------
if [[ -s "$OUT/onramp.h5ad" ]]; then
    echo ">>> [1/5] SKIP onramp h5ad (present)"
else
    echo ">>> [1/5] convert_seurat_to_h5ad.R -> onramp.h5ad (X=counts + logcounts)"
    $R "$REPO/bin/convert_seurat_to_h5ad.R" \
        --rds "$OUT/prepped.rds" --output "$OUT/onramp.h5ad" \
        --cell_type_key "$CT_KEY" --sample_key "$SAMPLE_KEY" \
        2>&1 | tee "$OUT/logs/01_onramp_h5ad.log"
fi
[[ -s "$OUT/onramp.h5ad" ]] || { echo "FATAL: onramp h5ad not produced"; exit 1; }

# --- 2. NATIVE-convention h5ad --------------------------------------------
# Replace `counts` with `data` so the converter's data==counts branch fires
# and no logcounts layer is written. X then holds the normalised values --
# the native scanpy convention -- with identical obs/var to arm 1.
if [[ -s "$OUT/native.h5ad" ]]; then
    echo ">>> [2/5] SKIP native h5ad (present)"
else
    echo ">>> [2/5] building native-convention object (counts <- data)"
    $R -e "
    suppressPackageStartupMessages({library(Seurat); library(Matrix)})
    o <- readRDS('$OUT/prepped.rds')
    a <- o[[DefaultAssay(o)]]
    dat <- SeuratObject::LayerData(a, layer='data')
    SeuratObject::LayerData(o[[DefaultAssay(o)]], layer='counts') <- dat
    chk <- o[[DefaultAssay(o)]]
    cat('counts max now:', max(SeuratObject::LayerData(chk, layer='counts')), '\n')
    saveRDS(o, '$OUT/native_conv.rds')
    " 2>&1 | tee "$OUT/logs/02a_native_obj.log"
    [[ -s "$OUT/native_conv.rds" ]] || { echo "FATAL: native obj failed"; exit 1; }

    echo ">>> [2/5] convert_seurat_to_h5ad.R -> native.h5ad (X=normalised)"
    $R "$REPO/bin/convert_seurat_to_h5ad.R" \
        --rds "$OUT/native_conv.rds" --output "$OUT/native.h5ad" \
        --cell_type_key "$CT_KEY" --sample_key "$SAMPLE_KEY" \
        2>&1 | tee "$OUT/logs/02b_native_h5ad.log"
fi
[[ -s "$OUT/native.h5ad" ]] || { echo "FATAL: native h5ad not produced"; exit 1; }

# Assert the two conventions actually came out different, or the gate is inert.
echo ">>> [2/5] verifying the two h5ads really differ in convention"
$R -e "
st_o <- rhdf5::h5ls('$OUT/onramp.h5ad')
st_n <- rhdf5::h5ls('$OUT/native.h5ad')
lo <- any(st_o\$group=='/layers' & st_o\$name=='logcounts')
ln <- any(st_n\$group=='/layers' & st_n\$name=='logcounts')
cat('onramp has logcounts:', lo, '  native has logcounts:', ln, '\n')
if (!lo) stop('FATAL: onramp h5ad has no logcounts layer -- gate would be inert')
if (ln)  stop('FATAL: native h5ad HAS a logcounts layer -- not the native convention')
cat('Convention assertion OK\n')
" 2>&1 | tee "$OUT/logs/02c_convention_check.log"
[[ ${PIPESTATUS[0]} -ne 0 ]] && { echo "FATAL: convention check failed"; exit 1; }

# --- 3+4. both CellChat arms ----------------------------------------------
rc_all=0
for ARM in native onramp; do
    if [[ -s "$OUT/$ARM/${ARM}_cellchat.rds" ]]; then
        echo ">>> [3/5] SKIP CellChat $ARM (result present)"
        continue
    fi
    echo ">>> [3/5] run_cellchat_infer.R -- $ARM arm"
    $R "$REPO/bin/run_cellchat_infer.R" \
        --input "$OUT/$ARM.h5ad" \
        --output_prefix "$OUT/$ARM/$ARM" \
        --cell_type_key cell_type \
        --species "$SPECIES" \
        --threads "$THREADS" \
        2>&1 | tee "$OUT/logs/03_cellchat_$ARM.log"
    [[ ${PIPESTATUS[0]} -ne 0 ]] && rc_all=1
done

NA_RDS="$OUT/native/native_cellchat.rds"
RT_RDS="$OUT/onramp/onramp_cellchat.rds"
if [[ ! -s "$NA_RDS" || ! -s "$RT_RDS" ]]; then
    echo "FATAL: one or both arms produced no cellchat object"
    echo "  native=$NA_RDS onramp=$RT_RDS"
    exit 1
fi

# --- 5. compare -----------------------------------------------------------
echo ">>> [5/5] compare"
$R "$REPO/tests/onramp_rds/compare_cellchat.R" \
    --native "$NA_RDS" --roundtrip "$RT_RDS" \
    --label "$(basename "$OUT")" \
    --out "$OUT/compare_cellchat.tsv" 2>&1 | tee "$OUT/logs/04_compare.log"
[[ ${PIPESTATUS[0]} -ne 0 ]] && rc_all=1

echo "=============================================================="
echo "finished : $(date)"
cat "$OUT/compare_cellchat.tsv" 2>/dev/null
[[ $rc_all -eq 0 ]] && echo "OVERALL: PASS" || echo "OVERALL: FAIL (see logs)"
exit $rc_all
