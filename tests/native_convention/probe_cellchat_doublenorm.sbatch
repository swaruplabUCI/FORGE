#!/bin/bash
#SBATCH --job-name=forge_probe_ccdn
#SBATCH -A vswarup_lab
#SBATCH --partition=highmem
#SBATCH --time=24:00:00
#SBATCH --mem=250GB
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL
# ============================================================================
# QUANTIFY-FIRST PROBE: does the native-path CellChat double normalisation
# actually change the inferred ligand-receptor communication?
#
# A MEASUREMENT, NOT A GATE. The arms are deliberately different; a FAIL is
# the expected shape of the answer and the MAGNITUDE is the result.
#
# THE DEFECT
#   bin/run_cellchat_per_condition.R:113 (live, via CELLCHAT_PER_CONDITION)
#   and bin/run_cellchat.R:52 (no nextflow reference; dead but wrong) take
#       data.input <- log1p(t(t(X) / colSums(X)) * 1e4)
#   unconditionally. On the NATIVE convention X is ALREADY log-normalised, so
#   this CPM+log1p's data that has already been normalised. Note colSums(X) of
#   logged values is not a library size at all -- it is the sum of logs -- so
#   the per-cell scale factor is a quantity with no meaning.
#
#   run_cellchat_infer.R does NOT have this defect: 4cf3975 gave it an
#   integrality test, so it uses a non-integral X as-is. The fix here is to
#   make the other two consumers agree with it.
#
# THE TWO ARMS, on ONE native-convention h5ad, same cells, same genes:
#   LEGACY    : X -> log1p(t(t(X)/colSums(X))*1e4)   [what ships today]
#   CORRECTED : X as-is                              [what infer.R already does]
#   Both then go through run_cellchat_infer.R, which takes the identical
#   non-integral-X branch for each, so CellChat itself is held constant and
#   the matrix is the only variable.
#
# Usage:
#   sbatch --export=ALL,H5AD=/path/native.h5ad,CT_KEY=scanvi_prediction,\
#N_CT=6,SPECIES=human,OUT=/path/outdir \
#     tests/native_convention/probe_cellchat_doublenorm.sbatch
# ============================================================================
set -uo pipefail

REPO=${REPO:-/dfs7/swaruplab/lesolano/src_FORGE}
SING=${SING:-/opt/apps/singularity/3.11.3/bin/singularity}
SIF=$REPO/singularity_cache/seurat_extended.sif
H5AD=${H5AD:?set H5AD=/path/to/native_integrated.h5ad}
CT_KEY=${CT_KEY:?set CT_KEY=<obs column>}
N_CT=${N_CT:-6}
SPECIES=${SPECIES:-human}
OUT=${OUT:?set OUT=/path/to/outdir}
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUT"/{legacy,corrected,logs}
export HOME=$(mktemp -d /tmp/ccdn_home_XXXXXX)
R="$SING exec -B /dfs7 -B /tmp $SIF Rscript"

echo "=============================================================="
echo "PROBE: CellChat native double normalisation"
echo "  h5ad   : $H5AD"
echo "  ct_key : $CT_KEY   n_ct: $N_CT   species: $SPECIES"
echo "  out    : $OUT"
echo "  started: $(date)"
echo "=============================================================="

# --- 1. build both matrices, quantify the distortion, write two h5ads -----
# Stage 1 prints the matrix-level answer BEFORE the expensive CellChat runs,
# so a decision is possible early if the distortion turns out to be trivial.
if [[ -s "$OUT/corrected.h5ad" && -s "$OUT/legacy.h5ad" ]]; then
    echo ">>> [1/3] SKIP matrix build (both h5ads present)"
else
    echo ">>> [1/3] build both matrices + quantify"
    $R -e "
suppressPackageStartupMessages({
    library(zellkonverter); library(SingleCellExperiment); library(Matrix)
})
sce <- readH5AD('$H5AD')
X <- assay(sce, 'X')
if (!inherits(X, 'dgCMatrix')) X <- as(as(X, 'CsparseMatrix'), 'dgCMatrix')
stopifnot(!all(X@x == floor(X@x)))   # must be the native convention
ct <- as.character(colData(sce)[['$CT_KEY']])
keep_ct <- names(sort(table(ct), decreasing = TRUE))[1:$N_CT]
sel <- ct %in% keep_ct
cat('cell types kept:', paste(keep_ct, collapse=' | '), '\n')
cat('cells kept     :', sum(sel), 'of', length(ct), '\n\n')
sce <- sce[, sel]; X <- X[, sel]; ct <- ct[sel]

# the two matrices
corrected <- X
ls_ <- Matrix::colSums(X)
legacy <- X
per_nz <- rep.int(ls_, diff(legacy@p))
legacy@x <- log1p(legacy@x / per_nz * 1e4)

cat('LEGACY    max', signif(max(legacy@x),6), ' CORRECTED max', signif(max(corrected@x),6), '\n')
cat('per-cell scale factor 1e4/colSums(X): min', signif(min(1e4/ls_),4),
    ' median', signif(median(1e4/ls_),4), ' max', signif(max(1e4/ls_),4), '\n\n')

# Element-wise agreement is nearly meaningless here because the transform is
# monotone WITHIN a cell. What CellChat actually consumes is a per-cell-group
# summary, so score THAT: per (cell type, gene) mean, both ways.
gm <- function(M) {
    t(sapply(keep_ct, function(k) Matrix::rowMeans(M[, ct == k, drop = FALSE])))
}
ga <- gm(legacy); gb <- gm(corrected)
cat('*** PER-(CELL TYPE, GENE) MEAN EXPRESSION -- what CellChat consumes ***\n')
cat(sprintf('  pearson  %.6f\n', cor(as.vector(ga), as.vector(gb))))
cat(sprintf('  spearman %.6f\n', cor(as.vector(ga), as.vector(gb), method='spearman')))
for (k in keep_ct) {
    a <- ga[k, ]; b <- gb[k, ]
    cat(sprintf('  %-22s within-group gene rank cor %.6f\n', k,
                cor(a, b, method='spearman')))
}
# the cross-cell-type ordering per gene is what drives differential signalling
pg <- sapply(seq_len(ncol(ga)), function(j) {
    if (sd(ga[,j])==0 || sd(gb[,j])==0) NA_real_ else cor(ga[,j], gb[,j])
})
cat(sprintf('\n  per-gene across-cell-type correlation: median %.6f, frac < 0.95: %.4f\n',
            median(pg, na.rm=TRUE), mean(pg < 0.95, na.rm=TRUE)))

assay(sce, 'X') <- corrected
writeH5AD(sce, '$OUT/corrected.h5ad', X_name = 'X')
assay(sce, 'X') <- legacy
writeH5AD(sce, '$OUT/legacy.h5ad', X_name = 'X')
cat('\nwrote both h5ads\n')
" 2>&1 | tee "$OUT/logs/01_build_quantify.log"
    [[ ${PIPESTATUS[0]} -ne 0 ]] && { echo "FATAL: matrix build failed"; exit 1; }
fi
[[ -s "$OUT/legacy.h5ad" && -s "$OUT/corrected.h5ad" ]] || { echo "FATAL: h5ads missing"; exit 1; }

# --- 2. CellChat through the SAME validated reader for both arms ----------
# run_cellchat_infer.R sees a non-integral X in both cases and uses it as-is,
# so it applies no normalisation of its own and CellChat is held constant.
for ARM in legacy corrected; do
    if [[ -s "$OUT/$ARM/${ARM}_cellchat.rds" ]]; then
        echo ">>> [2/3] SKIP CellChat $ARM (present)"; continue
    fi
    echo ">>> [2/3] CellChat -- $ARM"
    $R "$REPO/bin/run_cellchat_infer.R" \
        --input "$OUT/$ARM.h5ad" \
        --output_prefix "$OUT/$ARM/$ARM" \
        --cell_type_key "$CT_KEY" \
        --species "$SPECIES" \
        --threads "$THREADS" \
        2>&1 | tee "$OUT/logs/02_cellchat_$ARM.log"
done
A="$OUT/legacy/legacy_cellchat.rds"; B="$OUT/corrected/corrected_cellchat.rds"
[[ -s "$A" && -s "$B" ]] || { echo "FATAL: an arm produced no cellchat object"; exit 1; }

# --- 3. score -------------------------------------------------------------
echo ">>> [3/3] score  (native column = LEGACY, roundtrip column = CORRECTED)"
$R "$REPO/tests/onramp_rds/compare_cellchat.R" \
    --native "$A" --roundtrip "$B" --label "double_norm" \
    --out "$OUT/probe_cellchat.tsv" 2>&1 | tee "$OUT/logs/03_compare.log"

echo "=============================================================="
echo "finished : $(date)"
cat "$OUT/probe_cellchat.tsv" 2>/dev/null
echo
echo "REMINDER: a FAIL here is the EXPECTED shape. Read the MAGNITUDE."
grep -A12 'PER-(CELL TYPE, GENE)' "$OUT/logs/01_build_quantify.log" 2>/dev/null
echo "=============================================================="
exit 0
