#!/bin/bash
#SBATCH --job-name=forge_rt_hdwgcna
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=24:00:00
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL
# ============================================================================
# Tier 2: hdWGCNA round-trip equivalence for the Seurat .rds on-ramp.
#
# Tier 1 (compare.py) proves the converted MATRIX is bit-exact. It cannot
# prove an ANALYSIS on top of that matrix lands in the same place -- hdWGCNA
# soft-power selection and the WGCNA tree cut are threshold decisions, so a
# 1e-15 drift could in principle tip a module split.
#
# Two arms on the same cell type:
#   NATIVE    : run_hdwgcna_celltype.R straight off the original .rds
#   ROUNDTRIP : .rds -> convert_seurat_to_h5ad.R (new, rhdf5)
#                    -> convert_h5ad_to_seurat.R (existing, zellkonverter)
#                    -> run_hdwgcna_celltype.R
#
# This exercises BOTH converters, so it catches a round-trip asymmetry that
# Tier 1 is blind to.
#
# Usage:
#   sbatch --export=ALL,RDS=/path/x.rds,CT_KEY=subclass,N_CT=2,OUT=/path/out \
#          tests/onramp_rds/roundtrip_hdwgcna.sbatch
# ============================================================================
set -uo pipefail

REPO=${REPO:-/dfs7/swaruplab/lesolano/src_FORGE}
SING=${SING:-/opt/apps/singularity/3.11.3/bin/singularity}
SIF=$REPO/singularity_cache/seurat_extended.sif
RDS=${RDS:?set RDS=/path/to/input.rds}
CT_KEY=${CT_KEY:?set CT_KEY=<meta.data column>}
N_CT=${N_CT:-2}
OUT=${OUT:?set OUT=/path/to/outdir}
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT"/{native,roundtrip,logs}
export HOME=$(mktemp -d /tmp/rt_home_XXXXXX)
R="$SING exec -B /dfs7 -B /tmp $SIF Rscript"

echo "=============================================================="
echo "ROUND-TRIP hdWGCNA EQUIVALENCE"
echo "  rds      : $RDS"
echo "  ct_key   : $CT_KEY"
echo "  n_ct     : $N_CT  (most abundant cell types)"
echo "  out      : $OUT"
echo "  threads  : $THREADS"
echo "  started  : $(date)"
echo "=============================================================="

# --- 0. pick the N most abundant cell types -------------------------------
# Chosen by abundance rather than hardcoded, so this works on any object and
# we never compare networks built from too few cells to be meaningful.
$R -e "
suppressPackageStartupMessages(library(Seurat))
o <- readRDS('$RDS')
if (package_version(as.character(o@version)) < '5.0.0') o <- UpdateSeuratObject(o)
stopifnot('$CT_KEY' %in% colnames(o@meta.data))
t <- sort(table(as.character(o@meta.data[['$CT_KEY']])), decreasing=TRUE)
t <- t[names(t) != 'NA' & nzchar(names(t))]
writeLines(head(names(t), $N_CT), '$OUT/celltypes.txt')
cat('Selected:', paste(head(names(t), $N_CT), collapse=' | '), '\n')
cat('Counts  :', paste(head(as.integer(t), $N_CT), collapse=' | '), '\n')
" 2>&1 | tee "$OUT/logs/00_select.log"
[[ -s "$OUT/celltypes.txt" ]] || { echo "FATAL: no cell types selected"; exit 1; }

# --- 0b. Seurat version migration for the NATIVE arm ----------------------
# run_hdwgcna_celltype.R cannot read a pre-v5 Seurat object: the Allen
# reference is 3.0.0.9000 and has no `images` slot, so even basic accessors
# throw "no slot of name images" (handoff trap (f)).
#
# The ROUND-TRIP arm never hits this because convert_seurat_to_h5ad.R calls
# UpdateSeuratObject() internally. Applying the same migration to the native
# arm is what makes the two arms COMPARABLE: it isolates the variable actually
# under test (our converter round trip) from an incidental Seurat-version
# difference that the other arm already got for free.
#
# This is a test-harness concern only. In FORGE, run_hdwgcna_celltype.R is only
# ever fed CONVERT_H5AD_TO_SEURAT output -- a freshly built v5 object. It never
# receives a raw user .rds.
NATIVE_RDS="$RDS"
if [[ ! -s "$OUT/native_updated.rds" ]]; then
    $R -e "
    suppressPackageStartupMessages(library(Seurat))
    o <- readRDS('$RDS')
    v <- package_version(as.character(o@version))
    cat('Input Seurat version:', as.character(v), '\n')
    if (v < '5.0.0') {
        cat('Pre-v5 object -- running UpdateSeuratObject()\n')
        o <- UpdateSeuratObject(o)
        saveRDS(o, '$OUT/native_updated.rds')
        cat('Wrote $OUT/native_updated.rds\n')
    } else {
        cat('Already v5+ -- native arm will use the original .rds\n')
    }
    " 2>&1 | tee "$OUT/logs/00b_update.log"
fi
[[ -s "$OUT/native_updated.rds" ]] && NATIVE_RDS="$OUT/native_updated.rds"
echo "NATIVE arm will read: $NATIVE_RDS"

# --- 1. round trip through both converters --------------------------------
# Steps are cached on their outputs so a resubmit after a harness fix does not
# redo work that already succeeded.
if [[ -s "$OUT/rt.h5ad" ]]; then
    echo ">>> [1/4] SKIP convert_seurat_to_h5ad.R (rt.h5ad present)"
else
    echo ">>> [1/4] convert_seurat_to_h5ad.R (new, rhdf5)"
    $R "$REPO/bin/convert_seurat_to_h5ad.R" \
        --rds "$RDS" --output "$OUT/rt.h5ad" --cell_type_key "$CT_KEY" \
        2>&1 | tee "$OUT/logs/01_seurat_to_h5ad.log"
fi
[[ -s "$OUT/rt.h5ad" ]] || { echo "FATAL: h5ad not produced"; exit 1; }

if [[ -s "$OUT/rt.rds" ]]; then
    echo ">>> [2/4] SKIP convert_h5ad_to_seurat.R (rt.rds present)"
else
    echo ">>> [2/4] convert_h5ad_to_seurat.R (existing, zellkonverter)"
    $R "$REPO/bin/convert_h5ad_to_seurat.R" \
        --h5ad "$OUT/rt.h5ad" --output "$OUT/rt.rds" --cell_type_key cell_type \
        2>&1 | tee "$OUT/logs/02_h5ad_to_seurat.log"
fi
[[ -s "$OUT/rt.rds" ]] || { echo "FATAL: round-trip rds not produced"; exit 1; }

# --- 2+3. both hdWGCNA arms, per cell type --------------------------------
rc_all=0
while IFS= read -r CT; do
    [[ -z "$CT" ]] && continue
    SAFE=$(echo "$CT" | tr -c '[:alnum:]' '_')
    if ls "$OUT/native/$SAFE"/*_seurat_with_wgcna.rds >/dev/null 2>&1; then
        echo ">>> [3/4] SKIP hdWGCNA NATIVE    -- $CT (result present)"
    else
        echo ">>> [3/4] hdWGCNA NATIVE    -- $CT"
        mkdir -p "$OUT/native/$SAFE"
        ( cd "$OUT/native/$SAFE" && mkdir -p FiguresTables && \
          $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$NATIVE_RDS" \
             --cell_type "$CT" --cell_type_key "$CT_KEY" --output_dir . \
             --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/03_native_$SAFE.log"
    fi

    if ls "$OUT/roundtrip/$SAFE"/*_seurat_with_wgcna.rds >/dev/null 2>&1; then
        echo ">>> [3/4] SKIP hdWGCNA ROUNDTRIP -- $CT (result present)"
    else
        echo ">>> [3/4] hdWGCNA ROUNDTRIP -- $CT"
        mkdir -p "$OUT/roundtrip/$SAFE"
        ( cd "$OUT/roundtrip/$SAFE" && mkdir -p FiguresTables && \
          $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$OUT/rt.rds" \
             --cell_type "$CT" --cell_type_key cell_type --output_dir . \
             --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/03_rt_$SAFE.log"
    fi

    NA_RDS=$(ls "$OUT/native/$SAFE"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
    RT_RDS=$(ls "$OUT/roundtrip/$SAFE"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
    if [[ -z "$NA_RDS" || -z "$RT_RDS" ]]; then
        echo "SKIP $CT: one or both arms produced no wgcna object"
        echo "  native=$NA_RDS roundtrip=$RT_RDS"
        rc_all=1; continue
    fi

    echo ">>> [4/4] compare -- $CT"
    $R "$REPO/tests/onramp_rds/compare_hdwgcna.R" \
        --native "$NA_RDS" --roundtrip "$RT_RDS" --label "$CT" \
        --out "$OUT/compare_$SAFE.tsv" 2>&1 | tee "$OUT/logs/04_compare_$SAFE.log"
    [[ ${PIPESTATUS[0]} -ne 0 ]] && rc_all=1
done < "$OUT/celltypes.txt"

echo "=============================================================="
echo "finished : $(date)"
cat "$OUT"/compare_*.tsv 2>/dev/null | awk 'NR==1 || $1!="label"'
[[ $rc_all -eq 0 ]] && echo "OVERALL: PASS" || echo "OVERALL: FAIL (see logs)"
echo "=============================================================="
exit $rc_all
