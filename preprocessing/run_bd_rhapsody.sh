#!/bin/bash
#SBATCH --job-name=bd_preproc
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=16
#SBATCH --output=preprocessing/logs/out.bd_preproc_%j.log
#SBATCH --error=preprocessing/logs/err.bd_preproc_%j.log
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =========================================================================
# BD Rhapsody WTA+ATAC Preprocessing
#
# Runs the BD Rhapsody Docker-free install bundle to demultiplex and
# align BD multiome FASTQs. Produces:
#   - MEX matrix (RNA counts)
#   - fragments.tsv.gz (ATAC)
#   - Cell metrics and QC reports
#
# These outputs should be staged into the pipeline's data/ directory.
#
# Usage:
#   1. Copy this script to your project directory
#   2. Edit PREPROC_DIR, BUNDLE, CONDA_ENV, and CONFIG below
#   3. Create a CWL input YAML (see examples/brain_mm_bd.yml)
#   4. sbatch run_bd_rhapsody.sh
# =========================================================================

# ---- EDIT THESE ----
PREPROC_DIR="$(pwd)/preprocessing"
BUNDLE="/dfs7/swaruplab/lesolano/tools/bd_rhapsody_bundle"
CONDA_ENV="/dfs7/swaruplab/lesolano/tools/bd_rhapsody_cwl"
CONFIG="examples/brain_mm_bd.yml"   # Your CWL input YAML
# --------------------

cd "${PREPROC_DIR}" || exit 1

echo ""
echo "========================================="
echo "BD Rhapsody Preprocessing"
echo "  Started:  $(date)"
echo "  Workdir:  ${PREPROC_DIR}"
echo "  Bundle:   ${BUNDLE}"
echo "  Config:   ${CONFIG}"
echo "========================================="
echo ""

# =========================================================================
# Pre-flight checks
# =========================================================================
echo "Pre-flight checks:"

if [ ! -f "${CONFIG}" ]; then
    echo "  MISSING: ${CONFIG}"; exit 1
fi
echo "  OK: ${CONFIG}"

if [ ! -x "${BUNDLE}/rhapsody" ]; then
    echo "  MISSING: ${BUNDLE}/rhapsody"; exit 1
fi
echo "  OK: ${BUNDLE}/rhapsody"

echo ""

# =========================================================================
# Activate conda environment (cwl-runner)
# =========================================================================
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate "${CONDA_ENV}"
echo "  Python:     $(which python3)"
echo "  cwl-runner: $(cwl-runner --version 2>&1)"
echo ""

# =========================================================================
# Run BD Rhapsody pipeline
# =========================================================================
mkdir -p logs results

echo "========================================="
echo "LAUNCHING BD RHAPSODY PIPELINE"
echo "  $(date)"
echo "========================================="
echo ""

"${BUNDLE}/rhapsody" pipeline \
    --outdir "${PREPROC_DIR}/results" \
    "${CONFIG}"

exit_code=$?

echo ""
echo "========================================="
echo "BD Rhapsody pipeline completed"
echo "  Exit code: ${exit_code}"
echo "  Finished:  $(date)"
echo "  Results:   ${PREPROC_DIR}/results/"
echo "========================================="

if [ ${exit_code} -eq 0 ]; then
    echo ""
    echo "Output summary:"
    ls -lh "${PREPROC_DIR}/results/" 2>/dev/null | head -30
fi

exit ${exit_code}
