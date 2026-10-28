#!/bin/bash
#SBATCH --job-name=forge_determ_ctrl
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --mem=120GB
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL
# ============================================================================
# DETERMINISM CONTROL for the Tier 2 round-trip gate.
#
# roundtrip_hdwgcna.sbatch scores NATIVE vs ROUNDTRIP against fixed
# thresholds. That is only meaningful if two runs of the SAME arm on the SAME
# input agree: if hdWGCNA has a run-to-run noise floor, the thresholds are
# scoring noise and no converter work can ever turn the gate green.
#
# This runs run_hdwgcna_celltype.R TWICE on one identical object and scores
# the pair with the gate's own comparator.
#
#   exit 0 -> deterministic; any Tier 2 drift is REAL and attributable to the
#             converters, with no noise budget to hide in
#   exit 1 -> the gate is scoring noise; the thresholds, and probably the
#             gate's design, need rethinking before any converter is blamed
#
# Measured 2026-10-27 on the Allen reference (L6 IT, 1872 cells): all 9
# metrics exactly 1.0000, identical metacell counts (1005), identical module
# calls. The noise floor is ZERO. That is what licenses reading a
# native-vs-roundtrip score below 1.0 as converter signal -- and it is why
# the current thresholds (0.90 / 0.95 / 0.80) are far more permissive than
# the pipeline actually requires.
#
# Re-run this whenever hdWGCNA, WGCNA, Seurat or the container changes: the
# zero noise floor is a property of those versions, not a law.
#
# Usage:
#   sbatch --export=ALL,RDS=/path/obj.rds,CT_KEY=subclass,CT='L6 IT',OUT=/path/out \
#          tests/onramp_rds/determinism_control.sbatch
#
# RDS should be an object run_hdwgcna_celltype.R can read directly -- e.g.
# the native_prepped.rds that roundtrip_hdwgcna.sbatch step 0b materialises.
# ============================================================================
set -uo pipefail

REPO=${REPO:-/dfs7/swaruplab/lesolano/src_FORGE}
SING=${SING:-/opt/apps/singularity/3.11.3/bin/singularity}
SIF=$REPO/singularity_cache/seurat_extended.sif
RDS=${RDS:?set RDS=/path/to/object.rds}
CT=${CT:?set CT=<cell type>}
CT_KEY=${CT_KEY:?set CT_KEY=<meta.data column>}
OUT=${OUT:?set OUT=/path/to/outdir}
THREADS=${SLURM_CPUS_PER_TASK:-20}

mkdir -p "$OUT"/{repA,repB,logs}
export HOME=$(mktemp -d /tmp/determ_home_XXXXXX)
R="$SING exec -B /dfs7 -B /tmp $SIF Rscript"

echo "=============================================================="
echo "DETERMINISM CONTROL"
echo "  rds     : $RDS"
echo "  cell    : $CT  (key: $CT_KEY)"
echo "  out     : $OUT"
echo "  threads : $THREADS"
echo "=============================================================="

# Both replicates must see an identical thread count: WGCNA block decomposition
# is thread-count dependent, so varying it would confound the control.
for REP in repA repB; do
    if ls "$OUT/$REP"/*_seurat_with_wgcna.rds >/dev/null 2>&1; then
        echo ">>> SKIP $REP (result present)"; continue
    fi
    echo ">>> hdWGCNA replicate $REP -- $CT"
    ( cd "$OUT/$REP" && mkdir -p FiguresTables && \
      $R "$REPO/bin/run_hdwgcna_celltype.R" --seurat_rds "$RDS" \
         --cell_type "$CT" --cell_type_key "$CT_KEY" --output_dir . \
         --threads "$THREADS" ) 2>&1 | tee "$OUT/logs/$REP.log"
done

A=$(ls "$OUT/repA"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
B=$(ls "$OUT/repB"/*_seurat_with_wgcna.rds 2>/dev/null | head -1)
if [[ -z "$A" || -z "$B" ]]; then
    echo "FATAL: a replicate produced no wgcna object (A=$A B=$B)"; exit 1
fi

echo ">>> compare replicate A vs replicate B"
$R "$REPO/tests/onramp_rds/compare_hdwgcna.R" \
   --native "$A" --roundtrip "$B" --label "$CT (native vs native)" \
   --out "$OUT/compare_determinism.tsv" 2>&1 | tee "$OUT/logs/compare.log"
rc=${PIPESTATUS[0]}

echo "=============================================================="
cat "$OUT/compare_determinism.tsv" 2>/dev/null
if [[ $rc -eq 0 ]]; then
    echo "DETERMINISM CONTROL: PASS -- deterministic; Tier 2 drift is REAL"
else
    echo "DETERMINISM CONTROL: FAIL -- the gate is scoring run-to-run noise."
    echo "  Do NOT blame the converters on a Tier 2 failure until this passes."
fi
echo "=============================================================="
exit $rc
