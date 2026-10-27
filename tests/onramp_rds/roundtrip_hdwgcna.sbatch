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

# --- 1. round trip through both converters --------------------------------
echo ">>> [1/4] convert_seurat_to_h5ad.R (new, rhdf5)"
$R "$REPO/bin/convert_seurat_to_h5ad.R" \
    --rds "$RDS" --output "$OUT/rt.h5ad" --cell_type_key "$CT_KEY" \
    2>&1 | tee "$OUT/logs/01_seurat_to_h5ad.log"
[[ -s "$OUT/rt.h5ad" ]] || { echo "FATAL: h5ad not produced"; exit 1; }

echo ">>> [2/4] convert_h5ad_to_seurat.R (existing, zellkonverter)"
$R "$REPO/bin/convert_h5ad_to_seurat.R" \
    --h5ad "$OUT/rt.h5ad" --output "$OUT/rt.rds" --cell_type_key cell_type \
    2>&1 | tee "$OUT/logs/02_h5ad_to_seurat.log"
[[ -s "$OUT/rt.rds" ]] || { echo "FATAL: round-trip rds not produced"; exit 1; }

# --- 2+3. both hdWGCNA arms, per cell type --------------------------------
rc_all=0
while IFS= read -r CT; do
    [[ -z "$CT" ]] && continue
    SAFE=$(echo "$CT" | tr -c '[:alnum:]' '_')
    echo ">>> [3/4] hdWGCNA NATIVE    -- $CT"
    mkdir -p "$OUT/native/$SAFE"
    ( cd "$OUT/native/$SAFE" && mkdir -p FiguresTables && \
      $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$RDS" \
         --cell_type "$CT" --cell_type_key "$CT_KEY" --output_dir . \
         --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/03_native_$SAFE.log"

    echo ">>> [3/4] hdWGCNA ROUNDTRIP -- $CT"
    mkdir -p "$OUT/roundtrip/$SAFE"
    ( cd "$OUT/roundtrip/$SAFE" && mkdir -p FiguresTables && \
      $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$OUT/rt.rds" \
         --cell_type "$CT" --cell_type_key cell_type --output_dir . \
         --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/03_rt_$SAFE.log"

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
