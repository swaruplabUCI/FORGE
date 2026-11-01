#!/bin/bash
#SBATCH --job-name=forge_probe_vst
#SBATCH -A vswarup_lab
#SBATCH --partition=highmem
#SBATCH --time=24:00:00
#SBATCH --mem=250GB
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL
# ============================================================================
# QUANTIFY-FIRST PROBE: does the native convention's counts/data mapping
# actually change hdWGCNA results?
#
# THIS IS A MEASUREMENT, NOT A GATE. It is expected to FAIL the equivalence
# thresholds -- the two arms are deliberately different. The number that
# matters is HOW different, because that is what decides whether correcting
# the native path is worth changing established results over.
#
# THE DEFECT BEING MEASURED
#   On the native convention X is log-normalised and layers/counts holds the
#   real raw counts. convert_h5ad_to_seurat.R maps counts = X, data = X, so
#   the Seurat `counts` slot holds NORMALISED values and the raw counts go
#   unused. FindVariableFeatures(selection.method = "vst") reads `counts`, so
#   every native hdWGCNA run FORGE has ever done selected variable features
#   from normalised data. Wrong in principle; unknown in magnitude.
#
# THE TWO ARMS, on ONE native-convention h5ad:
#   CURRENT   : convert_h5ad_to_seurat.R --legacy_counts_mapping
#               -> counts = X (normalised), data = X   [pre-2026-10-30]
#   CORRECTED : convert_h5ad_to_seurat.R as shipped
#               -> counts = layers/counts (raw), data = X (normalised)
#   then run_hdwgcna_celltype.R on each, scored with compare_hdwgcna.R.
#
# NOTE: the arms swapped roles on 2026-10-30. The measurement below said the
# legacy mapping changed results materially (Oligodendrocyte 6 vs 8 modules),
# so the corrected mapping became the default and the flag now restores the
# old behaviour. This harness is kept as the regression instrument: re-running
# it must keep reproducing the same delta.
#
# Both arms come from the SAME converter invocation path, differing only in
# that flag, so nothing else can drift between them.
#
# HOW TO READ THE RESULT
#   metrics at ~1.0  -> the mapping barely matters. Correcting it is cheap and
#                       safe; the published record is not in question.
#   metrics far off  -> correcting it would move established results. That is
#                       a finding about the published record, not a licence to
#                       quietly change the default. Report it; decide then.
#   Note soft_power and n_modules especially: those are threshold decisions
#   and the most likely place for a real divergence to show.
#
# Usage:
#   sbatch --export=ALL,H5AD=/path/native.h5ad,CT_KEY=scanvi_prediction,\
#N_CT=2,OUT=/path/outdir tests/native_convention/probe_vst_counts.sbatch
# ============================================================================
set -uo pipefail

REPO=${REPO:-/dfs7/swaruplab/lesolano/src_FORGE}
SING=${SING:-/opt/apps/singularity/3.11.3/bin/singularity}
SIF=$REPO/singularity_cache/seurat_extended.sif
H5AD=${H5AD:?set H5AD=/path/to/native_integrated.h5ad}
CT_KEY=${CT_KEY:?set CT_KEY=<obs column>}
N_CT=${N_CT:-2}
OUT=${OUT:?set OUT=/path/to/outdir}
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUT"/{current,corrected,logs}
export HOME=$(mktemp -d /tmp/probe_home_XXXXXX)
R="$SING exec -B /dfs7 -B /tmp $SIF Rscript"

echo "=============================================================="
echo "PROBE: native counts/data mapping -> hdWGCNA"
echo "  h5ad    : $H5AD"
echo "  ct_key  : $CT_KEY"
echo "  n_ct    : $N_CT"
echo "  out     : $OUT"
echo "  started : $(date)"
echo "=============================================================="

# --- 0. assert the input really is the native convention ------------------
# The whole probe is meaningless on an on-ramped h5ad: there X is already raw
# and there is no defect to measure. Fail loudly rather than produce a
# confidently wrong number.
$R -e "
suppressPackageStartupMessages(library(rhdf5))
f <- '$H5AD'
x <- as.numeric(h5read(f, '/X/data'))
cat('X: max', max(x), ' integral', all(x == floor(x)), '\n')
has_counts <- tryCatch({ h5read(f, '/layers/counts/data'); TRUE },
                       error = function(e) FALSE)
cat('layers/counts present:', has_counts, '\n')
if (all(x == floor(x)))
    stop('X is integral -- this is NOT a native-convention h5ad. The probe ',
         'only means something where X is log-normalised and layers/counts ',
         'holds the raw counts.')
if (!has_counts)
    stop('No layers/counts -- nothing to map the corrected arm onto.')
cat('Native convention confirmed.\n')
" 2>&1 | tee "$OUT/logs/00_assert.log"
[[ ${PIPESTATUS[0]} -ne 0 ]] && { echo "FATAL: input is not native convention"; exit 1; }

# --- 1. build both objects ------------------------------------------------
# Cached on their outputs so a resubmit after a harness fix does not redo the
# conversions. Delete $OUT to force a rebuild.
if [[ -s "$OUT/current.rds" ]]; then
    echo ">>> [1/4] SKIP current.rds (present)"
else
    echo ">>> [1/4] convert -- LEGACY mapping (pre-2026-10-30 behaviour)"
    $R "$REPO/bin/convert_h5ad_to_seurat.R" \
        --h5ad "$H5AD" --output "$OUT/current.rds" --cell_type_key "$CT_KEY" \
        --legacy_counts_mapping \
        2>&1 | tee "$OUT/logs/01_current.log"
fi
[[ -s "$OUT/current.rds" ]] || { echo "FATAL: legacy arm not built"; exit 1; }

if [[ -s "$OUT/corrected.rds" ]]; then
    echo ">>> [2/4] SKIP corrected.rds (present)"
else
    echo ">>> [2/4] convert -- CORRECTED mapping (the shipped default)"
    $R "$REPO/bin/convert_h5ad_to_seurat.R" \
        --h5ad "$H5AD" --output "$OUT/corrected.rds" --cell_type_key "$CT_KEY" \
        2>&1 | tee "$OUT/logs/02_corrected.log"
fi
[[ -s "$OUT/corrected.rds" ]] || { echo "FATAL: corrected arm not built"; exit 1; }

# --- 2. confirm the arms actually differ, and in the expected direction ----
# A probe whose two arms are accidentally identical reports a reassuring 1.0
# and means nothing. Prove the difference exists before spending the hours.
$R -e "
suppressPackageStartupMessages(library(Seurat))
a <- readRDS('$OUT/current.rds'); b <- readRDS('$OUT/corrected.rds')
ca <- GetAssayData(a, layer='counts'); cb <- GetAssayData(b, layer='counts')
da <- GetAssayData(a, layer='data');   db <- GetAssayData(b, layer='data')
isint <- function(m) all(m@x == floor(m@x))
cat('CURRENT   counts: max', max(ca@x), ' integral', isint(ca), '\n')
cat('CORRECTED counts: max', max(cb@x), ' integral', isint(cb), '\n')
cat('CURRENT   data  : max', max(da@x), ' integral', isint(da), '\n')
cat('CORRECTED data  : max', max(db@x), ' integral', isint(db), '\n')
if (isint(ca)) stop('CURRENT counts is integral -- the shipped mapping did ',
                    'not put normalised values in counts. Probe assumption broken.')
if (!isint(cb)) stop('CORRECTED counts is NOT integral -- the flag did not ',
                     'pick up layers/counts. Probe assumption broken.')
cat('Arms differ as expected.\n')
" 2>&1 | tee "$OUT/logs/03_arms_differ.log"
[[ ${PIPESTATUS[0]} -ne 0 ]] && { echo "FATAL: probe assumptions broken"; exit 1; }

# --- 3. pick cell types and run hdWGCNA on both arms ----------------------
$R -e "
suppressPackageStartupMessages(library(Seurat))
o <- readRDS('$OUT/current.rds')
t <- sort(table(as.character(o@meta.data[['$CT_KEY']])), decreasing=TRUE)
t <- t[names(t) != 'NA' & nzchar(names(t))]
writeLines(head(names(t), $N_CT), '$OUT/celltypes.txt')
cat('Selected:', paste(head(names(t), $N_CT), collapse=' | '), '\n')
cat('Counts  :', paste(head(as.integer(t), $N_CT), collapse=' | '), '\n')
" 2>&1 | tee "$OUT/logs/04_select.log"
[[ -s "$OUT/celltypes.txt" ]] || { echo "FATAL: no cell types"; exit 1; }

# --- 3b. the mechanism, measured the way the pipeline actually does it -----
# run_hdwgcna_celltype.R subsets to ONE cell type, then calls
# FindVariableFeatures(selection.method="vst", nfeatures=12000). An earlier
# version of this probe measured the overlap on the WHOLE object at
# nfeatures=2000 and got a number that looked right for the wrong reasons.
# Measure what the pipeline measures, or do not measure it.
#
# Report the rank correlation too: RunPCA consumes VariableFeatures in order,
# so a reshuffle of the shared genes matters even when membership does not.
$R -e "
suppressPackageStartupMessages(library(Seurat))
a <- readRDS('$OUT/current.rds'); b <- readRDS('$OUT/corrected.rds')
N <- nrow(a); K <- 12000
cat(sprintf('genes=%d  nfeatures=%d  jaccard expected from two INDEPENDENT draws: %.4f\n\n',
            N, K, (K*K/N)/(2*K - K*K/N)))
for (ct in readLines('$OUT/celltypes.txt')) {
  sa <- a[, a@meta.data[['$CT_KEY']] == ct]
  sb <- b[, b@meta.data[['$CT_KEY']] == ct]
  va <- VariableFeatures(FindVariableFeatures(sa, selection.method='vst',
                         nfeatures=K, verbose=FALSE))
  vb <- VariableFeatures(FindVariableFeatures(sb, selection.method='vst',
                         nfeatures=K, verbose=FALSE))
  i <- length(intersect(va, vb))
  ra <- match(intersect(va, vb), va); rb <- match(intersect(va, vb), vb)
  cat(sprintf('%-20s n=%6d  shared=%5d/%d  jaccard=%.4f  rank_cor=%.4f\n',
              ct, ncol(sa), i, K, i/length(union(va, vb)),
              cor(ra, rb, method='spearman')))
}
" 2>&1 | tee "$OUT/logs/07_variable_features.log"

rc_all=0
while IFS= read -r CT; do
    [[ -z "$CT" ]] && continue
    SAFE=$(echo "$CT" | tr -c '[:alnum:]' '_')
    for ARM in current corrected; do
        if ls "$OUT/$ARM/$SAFE"/*_seurat_with_wgcna.rds >/dev/null 2>&1; then
            echo ">>> [3/4] SKIP hdWGCNA $ARM -- $CT"
        else
            echo ">>> [3/4] hdWGCNA $ARM -- $CT"
            mkdir -p "$OUT/$ARM/$SAFE"
            ( cd "$OUT/$ARM/$SAFE" && mkdir -p FiguresTables && \
              $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$OUT/$ARM.rds" \
                 --cell_type "$CT" --cell_type_key "$CT_KEY" --output_dir . \
                 --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/05_${ARM}_$SAFE.log"
        fi
    done

    A=$(ls "$OUT/current/$SAFE"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
    B=$(ls "$OUT/corrected/$SAFE"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
    if [[ -z "$A" || -z "$B" ]]; then
        echo "SKIP $CT: an arm produced no wgcna object (current=$A corrected=$B)"
        rc_all=1; continue
    fi

    echo ">>> [4/4] score -- $CT  (current = native, corrected = roundtrip column)"
    $R "$REPO/tests/onramp_rds/compare_hdwgcna.R" \
        --native "$A" --roundtrip "$B" --label "$CT" \
        --out "$OUT/probe_$SAFE.tsv" 2>&1 | tee "$OUT/logs/06_score_$SAFE.log"
done < "$OUT/celltypes.txt"

echo "=============================================================="
echo "finished : $(date)"
cat "$OUT"/probe_*.tsv 2>/dev/null | awk 'NR==1 || $1!="label"'
echo
echo "REMINDER: a FAIL here is the EXPECTED shape of the answer. The arms are"
echo "deliberately different. Read the MAGNITUDE, not the verdict column."
cat "$OUT/logs/07_variable_features.log" 2>/dev/null | tail -5
echo "=============================================================="
exit 0
