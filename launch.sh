#!/bin/bash
#SBATCH --job-name=forge_ad_mm_10x_r5
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --output=/dfs7/swaruplab/lesolano/FORGE/AD_Mm_10X_r5/logs/out.forge_ad_mm_10x_%j.log
#SBATCH --error=/dfs7/swaruplab/lesolano/FORGE/AD_Mm_10X_r5/logs/err.forge_ad_mm_10x_%j.log
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =========================================================================
# FORGE: 10x Multiome Alzheimer's Mouse Model (CRND8)
#
# SELF-CONTAINED INSTANCE: All pipeline code, containers, and data
# are in this directory.
#
# KEY SETTINGS:
#   - 12 samples (6 TG CRND8 + 6 WT), true multiome (RNA + ATAC)
#   - RNA Path B (CellTypist direct, no scANVI reference atlas)
#   - Developing_Mouse_Brain.pkl CellTypist model
#   - Brain tissue type for ATAC annotation
#   - Differential enabled: WT vs TG
#   - Resource tier: medium
#
# Usage:
#   sbatch launch.sh              # full production run
#   sbatch launch.sh --dry-run    # dry-run preview only
# =========================================================================

DRY_RUN=""
RESUME="-resume"

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run) DRY_RUN="true"; shift ;;
        --no-resume) RESUME=""; shift ;;
        *) shift ;;
    esac
done

# =========================================================================
# Environment
# =========================================================================
PROJECT_DIR="/dfs7/swaruplab/lesolano/FORGE/AD_Mm_10X_r5"
DATASET_CONFIG="configs/datasets/ad_mm_10x.config"

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
# Resource tier
# =========================================================================
RESOURCE_TIER="medium"
CONFIG_TIER=$(grep 'resource_tier' "$DATASET_CONFIG" | grep -oP "'[^']+'" | tr -d "'" | head -1)
if [[ -n "$CONFIG_TIER" && "$CONFIG_TIER" != "auto" ]]; then
    RESOURCE_TIER="$CONFIG_TIER"
fi

# =========================================================================
# Pre-flight checks
# =========================================================================
echo ""
echo "========================================="
echo "FORGE: 10x Multiome AD Mouse (CRND8)"
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
for item in main.nf nextflow.config "$DATASET_CONFIG" ad_mm_10x_manifest.csv; do
    if [ -f "${item}" ]; then
        echo "  OK: ${item}"
    else
        echo "  MISSING: ${item}"
        exit 1
    fi
done
for dir in modules bin configs singularity_cache data; do
    if [ -d "${dir}" ]; then
        echo "  OK: ${dir}/"
    else
        echo "  MISSING: ${dir}/"
        exit 1
    fi
done
echo ""

# Container check
echo "Container check:"
CONTAINER_SIFS=$(grep -A20 'containers = \[' nextflow.config \
    | grep '\.sif"' \
    | grep -oP 'singularity_cache/[^"]+' \
    | sort -u)
if [[ -z "$CONTAINER_SIFS" ]]; then
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
    echo "  ERROR: Missing containers — aborting"
    exit 1
fi
echo ""

# Data file check — verify all 12 samples have RNA + ATAC + .tbi
echo "Data file check:"
DATA_DIR="${PROJECT_DIR}/data"
DATA_OK=true
for sample in AD_17p9_rep4 AD_17p9_rep5 AD_2p5_rep2 AD_2p5_rep3 AD_5p7_rep2 AD_5p7_rep6 \
              WT_13p4_rep2 WT_13p4_rep5 WT_2p5_rep2 WT_2p5_rep7 WT_5p7_rep2 WT_5p7_rep3; do
    for f in "${DATA_DIR}/${sample}_raw_feature_bc_matrix.h5" \
             "${DATA_DIR}/${sample}_atac_fragments.tsv.gz" \
             "${DATA_DIR}/${sample}_atac_fragments.tsv.gz.tbi"; do
        if [ -f "$f" ]; then
            echo "  OK: $(basename $f)"
        else
            echo "  MISSING: $f"
            DATA_OK=false
        fi
    done
done
if [[ "$DATA_OK" != "true" ]]; then
    echo "  ERROR: Missing data files — aborting"
    exit 1
fi
echo ""

# Manifest summary
echo "Manifest summary (ad_mm_10x_manifest.csv):"
echo "  Total rows:  $(tail -n +2 ad_mm_10x_manifest.csv | wc -l)"
echo "  Lane rows:   $(grep ',lane,' ad_mm_10x_manifest.csv | wc -l)"
echo "  WT samples:  $(grep ',WT,' ad_mm_10x_manifest.csv | wc -l)"
echo "  TG samples:  $(grep ',TG,' ad_mm_10x_manifest.csv | wc -l)"
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
  -with-dag      ${OUTDIR}/pipeline_info/nextflow_dag.pdf"

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

    dry_exit=$?
    echo ""
    echo "Dry-run complete. Exit code: ${dry_exit}"
    echo "Log: ${PROJECT_DIR}/logs/dry_run.log"
    echo ""
    echo "If passed, launch for real: sbatch launch.sh"
    exit ${dry_exit}

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
