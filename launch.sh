#!/bin/bash
#SBATCH --job-name=refactor5
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/out.refactor5_%j.log
#SBATCH --error=logs/err.refactor5_%j.log
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =========================================================================
# REFACTOR5: Unified Multiomics Pipeline v3.0.0
#
# Usage:
#   sbatch launch.sh -c configs/datasets/pbmc_10x_10k.config
#   sbatch launch.sh -c configs/datasets/bd_90plus_brain.config
#   sbatch launch.sh -c configs/datasets/my_dataset.config --dry-run
#
# Arguments:
#   -c <config>     Dataset config file (REQUIRED)
#   --dry-run       Preview mode (validates config + channels only)
#   --resume        Resume from previous run (default: enabled)
#   --no-resume     Start fresh (no cache)
# =========================================================================

# Parse arguments
DATASET_CONFIG=""
DRY_RUN=""
RESUME="-resume"
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -c) DATASET_CONFIG="$2"; shift 2 ;;
        --dry-run) DRY_RUN="true"; shift ;;
        --no-resume) RESUME=""; shift ;;
        *) EXTRA_ARGS="$EXTRA_ARGS $1"; shift ;;
    esac
done

if [[ -z "$DATASET_CONFIG" ]]; then
    echo "ERROR: Dataset config required. Usage: sbatch launch.sh -c configs/datasets/<dataset>.config"
    exit 1
fi

# =========================================================================
# Environment
# =========================================================================
PROJECT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "${PROJECT_DIR}" || exit 1
mkdir -p logs work results

module load singularity 2>/dev/null || true
export PATH="/dfs7/swaruplab/lesolano/tools:$PATH"

export NXF_WORK="${PROJECT_DIR}/work"
export TMPDIR=/tmp
export NXF_TEMP=/tmp
export HDF5_USE_FILE_LOCKING="FALSE"
export SINGULARITY_BINDPATH="/dfs7,/tmp"

# =========================================================================
# Auto-detect resource tier from manifest
# =========================================================================
MANIFEST=$(grep 'metadata_file' "$DATASET_CONFIG" | grep -oP "'[^']+'" | tr -d "'" | head -1)
if [[ -n "$MANIFEST" && -f "$MANIFEST" ]]; then
    # FIX-F10: Exclude blank lines from sample count
    N_SAMPLES=$(tail -n +2 "$MANIFEST" | grep -c '[^[:space:]]')
    if [[ $N_SAMPLES -le 5 ]]; then
        RESOURCE_TIER="small"
    elif [[ $N_SAMPLES -le 50 ]]; then
        RESOURCE_TIER="medium"
    else
        RESOURCE_TIER="large"
    fi
    echo "Auto-detected resource tier: ${RESOURCE_TIER} (${N_SAMPLES} samples)"
else
    RESOURCE_TIER="small"
    echo "Could not read manifest — defaulting to resource_tier=small"
fi

# Check if dataset config overrides resource_tier
CONFIG_TIER=$(grep 'resource_tier' "$DATASET_CONFIG" | grep -oP "'[^']+'" | tr -d "'" | head -1)
if [[ -n "$CONFIG_TIER" && "$CONFIG_TIER" != "auto" ]]; then
    RESOURCE_TIER="$CONFIG_TIER"
    echo "Dataset config overrides resource tier: ${RESOURCE_TIER}"
fi

# =========================================================================
# Pre-flight checks
# =========================================================================
echo ""
echo "========================================="
echo "REFACTOR5: Unified Multiomics Pipeline v3.0.0"
echo "  Dataset config: ${DATASET_CONFIG}"
echo "  Resource tier:  ${RESOURCE_TIER}"
echo "  Project dir:    ${PROJECT_DIR}"
echo "  Resume:         ${RESUME:-disabled}"
if [[ -n "$DRY_RUN" ]]; then
    echo "  *** DRY-RUN MODE ***"
fi
echo "========================================="
echo ""

echo "Project structure check:"
for item in main.nf nextflow.config "$DATASET_CONFIG"; do
    if [ -f "${item}" ]; then
        echo "  OK: ${item}"
    else
        echo "  MISSING: ${item}"
        exit 1
    fi
done
for dir in modules bin; do
    if [ -d "${dir}" ]; then
        echo "  OK: ${dir}/"
    else
        echo "  MISSING: ${dir}/"
        exit 1
    fi
done
echo ""

# FIX-R1-7: Derive container list from nextflow.config params.containers map
# instead of hardcoding. Nextflow also validates at startup (validateStartupParams).
echo "Container check:"
CONTAINER_SIFS=$(grep -A20 'containers = \[' nextflow.config \
    | grep '\.sif"' \
    | grep -oP 'singularity_cache/[^"]+' \
    | sort -u)
if [[ -z "$CONTAINER_SIFS" ]]; then
    echo "  WARNING: Could not parse container list from nextflow.config; falling back to known list"
    CONTAINER_SIFS="singularity_cache/scgpu_extended.sif singularity_cache/snapatac_extended.sif singularity_cache/seurat_extended.sif singularity_cache/cicero.sif singularity_cache/scenicplus.sif"
fi
CONTAINER_OK=true
for sif in $CONTAINER_SIFS; do
    if [ -f "$sif" ]; then
        echo "  OK: $(basename "$sif")"
    else
        echo "  MISSING: $sif"
        CONTAINER_OK=false
    fi
done
if [[ "$CONTAINER_OK" != "true" ]]; then
    echo "  ERROR: Missing containers. Pull/build before launching."
    exit 1
fi
echo ""

# =========================================================================
# Nextflow command
# =========================================================================
OUTDIR="${PROJECT_DIR}/results"

NF_CMD="nextflow run main.nf \
  -c ${DATASET_CONFIG} \
  -profile cluster,gpu,singularity \
  ${RESUME} \
  --resource_tier ${RESOURCE_TIER} \
  --outdir ${OUTDIR} \
  -with-report   ${OUTDIR}/pipeline_info/nextflow_report.html \
  -with-timeline ${OUTDIR}/pipeline_info/nextflow_timeline.html \
  -with-trace    ${OUTDIR}/pipeline_info/trace.tsv \
  -with-dag      ${OUTDIR}/pipeline_info/nextflow_dag.pdf \
  ${EXTRA_ARGS}"

if [[ -n "$DRY_RUN" ]]; then
    echo "========================================="
    echo "DRY-RUN PREVIEW"
    echo "========================================="
    echo "Command: $NF_CMD"
    echo ""

    nextflow run main.nf \
      -c "${DATASET_CONFIG}" \
      -profile cluster,singularity \
      -preview \
      --resource_tier "${RESOURCE_TIER}" \
      --outdir "${OUTDIR}" \
      2>&1 | tee logs/dry_run.log

    echo ""
    echo "Dry-run complete. Exit code: $?"
    echo "If passed, launch for real: sbatch launch.sh -c ${DATASET_CONFIG}"

else
    echo "========================================="
    echo "LAUNCHING NEXTFLOW (PRODUCTION)"
    echo "========================================="
    echo "Start: $(date)"
    echo ""

    eval "$NF_CMD"
    exit_code=$?

    echo ""
    echo "========================================="
    echo "Pipeline completed at $(date)"
    echo "Exit code: ${exit_code}"
    echo "Results:   ${OUTDIR}/"
    echo "========================================="
    exit ${exit_code}
fi
