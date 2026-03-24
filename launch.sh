#!/bin/bash
#SBATCH --job-name=v4_10x_pbmc
#SBATCH -A vswarup_lab
#SBATCH --partition=standard
#SBATCH --time=8:00:00  # capped from 72h for maintenance
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
#SBATCH --output=/dfs7/swaruplab/lesolano/nf_refactor4/PBMC_10x_nf/logs/out.v4_10x_pbmc_%j.log
#SBATCH --error=/dfs7/swaruplab/lesolano/nf_refactor4/PBMC_10x_nf/logs/err.v4_10x_pbmc_%j.log
#SBATCH --mail-user=lesolano@uci.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# =========================================================================
# V4 UNIFIED MULTIOMICS PIPELINE v2.2.0 — 10x Genomics 10k Human PBMC
#
# SELF-CONTAINED INSTANCE: This directory is an isolated Nextflow project.
# main.nf, modules/, bin/, and singularity_cache/ are copies from the
# parent pipeline. References are accessed via absolute paths in
# nextflow.config — no symlinks needed.
#
# This isolation means:
#   - .nextflow.log stays in PBMC_10x_nf/ (no cross-contamination)
#   - NXF_WORK is local to this project
#   - Can run simultaneously with BD pipeline in nf_refactor4/
#
# KEY DIFFERENCES FROM BD:
#   - Path B (CellTypist direct, no scANVI reference atlas)
#   - Single sample, no demultiplexing
#   - Immune_All_Low.pkl CellTypist model (vs Pan_Fetal_Human for brain)
#   - PBMC tissue type for ATAC annotation
#   - Scaled-down resources (~10k cells vs ~141 BD samples)
#
# Usage:
#   sbatch Launch_10x_pbmc.sh              # full production run
#   sbatch Launch_10x_pbmc.sh --dry-run    # dry-run preview only
# =========================================================================

DRY_RUN="${1:-}"

echo "========================================="
echo "V4 UNIFIED MULTIOMICS PIPELINE v2.2.0"
echo "10x Genomics 10k Human PBMC Multiome"
echo "  RNA Path:    B (CellTypist direct)"
echo "  CellTypist:  Immune_All_Low.pkl"
echo "  ATAC tissue: PBMC"
echo "  Modules:     Full pipeline (no differential)"
if [[ "$DRY_RUN" == "--dry-run" ]]; then
    echo "  *** DRY-RUN MODE ***"
fi
echo "========================================="
echo "Start time: $(date)"
echo ""

# ------------------------------
# Environment
# ------------------------------
source /pub/lesolano/miniconda3/etc/profile.d/conda.sh
module load singularity 2>/dev/null || true
export PATH="/dfs7/swaruplab/lesolano/tools:$PATH"

# ------------------------------
# Working directory
# This IS the project root — main.nf, nextflow.config, modules/, and
# bin/ all live here. Nextflow's projectDir resolves to this directory.
# Containers and references are accessed via absolute paths in the config.
# ------------------------------
PROJECT_DIR="/dfs7/swaruplab/lesolano/nf_refactor4/PBMC_10x_nf"
PARENT_DIR="/dfs7/swaruplab/lesolano/nf_refactor4"
cd "${PROJECT_DIR}" || exit 1
mkdir -p logs work data

# ------------------------------
# Nextflow runtime env
# ------------------------------
export NXF_WORK="${PROJECT_DIR}/work"
export TMPDIR=/tmp
export NXF_TEMP=/tmp
export HDF5_USE_FILE_LOCKING="FALSE"
export SINGULARITY_BINDPATH="/dfs7,/tmp"

# ------------------------------
# Inputs / outputs
# ------------------------------
PIPELINE_NF="${PROJECT_DIR}/main.nf"
OUTDIR="${PROJECT_DIR}/results_10x_pbmc"

echo "Configuration:"
echo "  Project dir:  ${PROJECT_DIR}"
echo "  Parent dir:   ${PARENT_DIR}"
echo "  Pipeline:     ${PIPELINE_NF}"
echo "  Config:       ${PROJECT_DIR}/nextflow.config (fully self-contained)"
echo "  Manifest:     ${PROJECT_DIR}/10k_pbmc_manifest.csv"
echo "  Mode:         full multiome (RNA + ATAC + integration)"
echo "  Outdir:       ${OUTDIR}"
echo "  NXF_WORK:     ${NXF_WORK}"
echo ""

# ------------------------------
# Verify local project assets
# ------------------------------
echo "Project structure check:"
for item in main.nf nextflow.config 10k_pbmc_manifest.csv; do
    if [ -f "${PROJECT_DIR}/${item}" ]; then
        echo "  OK: ${item}"
    else
        echo "  MISSING: ${item}"
        exit 1
    fi
done
for dir in modules bin; do
    if [ -d "${PROJECT_DIR}/${dir}" ]; then
        echo "  OK: ${dir}/"
    else
        echo "  MISSING: ${dir}/ — copy from parent: cp -r ${PARENT_DIR}/${dir} ."
        exit 1
    fi
done
echo ""

# NOTE: No parent config dependency — nextflow.config is fully self-contained.
# The includeConfig '../nextflow.config' was removed in v1.3.0.
echo ""

# Verify test data files
echo "Test data check:"
DATA_DIR="${PROJECT_DIR}/data"
for f in \
    "${DATA_DIR}/10k_PBMC_Multiome_nextgem_Chromium_X_raw_feature_bc_matrix.h5" \
    "${DATA_DIR}/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz" \
    "${DATA_DIR}/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz.tbi"
do
    if [ -f "$f" ]; then
        echo "  OK: $(basename $f)"
    else
        echo "  MISSING: $f — aborting"
        exit 1
    fi
done
echo ""

# Verify containers (local singularity_cache/ — copied from parent)
echo "Container check:"
CACHE_DIR="${PROJECT_DIR}/singularity_cache"
for c in scgpu_extended snapatac_extended seurat_extended cicero scenicplus; do
    if [ -f "${CACHE_DIR}/${c}.sif" ]; then
        echo "  OK: ${c}.sif"
    else
        echo "  MISSING: ${CACHE_DIR}/${c}.sif — aborting"
        exit 1
    fi
done
echo ""

# ------------------------------
# Quick manifest stats
# ------------------------------
echo "Manifest summary (10k_pbmc_manifest.csv):"
echo "  Total rows:  $(tail -n +2 ${PROJECT_DIR}/10k_pbmc_manifest.csv | wc -l)"
echo "  Lane rows:   $(grep ',lane,' ${PROJECT_DIR}/10k_pbmc_manifest.csv | wc -l)"
echo "  Demux rows:  $(grep ',demux,' ${PROJECT_DIR}/10k_pbmc_manifest.csv | wc -l)"
echo "  Batches:     $(tail -n +2 ${PROJECT_DIR}/10k_pbmc_manifest.csv | cut -d',' -f2 | sort -u | tr '\n' ' ')"
echo ""

# =========================================================================
# NEXTFLOW COMMAND
# Running from PROJECT_DIR with local main.nf — no -c overlay needed.
# nextflow.config in this directory is auto-detected by Nextflow.
# projectDir = PROJECT_DIR, so ${projectDir}/modules/ and ${projectDir}/bin/
# resolve to the local copies. Container and reference paths are absolute
# in the config, pointing to the parent's singularity_cache/ and references/.
# .nextflow.log stays in this directory (no cross-contamination with BD).
# =========================================================================

NF_CMD="nextflow run ${PIPELINE_NF} \
  -profile cluster,gpu,singularity \
  -resume ccb23db0-f705-4636-8645-e043f8d50af9 \
  --outdir ${OUTDIR} \
  -with-report   ${OUTDIR}/pipeline_info/nextflow_report.html \
  -with-timeline ${OUTDIR}/pipeline_info/nextflow_timeline.html \
  -with-trace    ${OUTDIR}/pipeline_info/trace.tsv \
  -with-dag      ${OUTDIR}/pipeline_info/nextflow_dag.pdf"

if [[ "$DRY_RUN" == "--dry-run" ]]; then
    echo "========================================="
    echo "DRY-RUN PREVIEW"
    echo "========================================="
    echo ""
    echo "Command that would be executed:"
    echo "  $NF_CMD"
    echo ""

    # Run Nextflow with -preview to validate config + channels without executing
    echo "Running Nextflow config validation (-preview)..."
    echo ""

    nextflow run "${PIPELINE_NF}" \
      -profile cluster,singularity \
      -preview \
      --outdir "${OUTDIR}" \
      2>&1 | tee "${PROJECT_DIR}/logs/dry_run_10x_pbmc.log"

    dry_exit=$?

    echo ""
    echo "========================================="
    echo "DRY-RUN COMPLETE"
    echo "  Exit code: ${dry_exit}"
    echo "  Log: ${PROJECT_DIR}/logs/dry_run_10x_pbmc.log"
    echo ""
    echo "If the dry run passed, launch the real run:"
    echo "  sbatch Launch_10x_pbmc.sh"
    echo "========================================="

    exit ${dry_exit}

else
    echo "========================================="
    echo "LAUNCHING NEXTFLOW (PRODUCTION)"
    echo "========================================="

    eval "$NF_CMD"

    exit_code=$?

    echo ""
    echo "========================================="
    echo "Pipeline completed at $(date)"
    echo "Exit code: ${exit_code}"
    echo "Outdir: ${OUTDIR}"
    echo ""
    echo "Check outputs:"
    echo "  RNA:         ${OUTDIR}/rna/"
    echo "  ATAC:        ${OUTDIR}/atac/"
    echo "  Multiome:    ${OUTDIR}/multiome/"
    echo "  Regulatory:  ${OUTDIR}/cicero/"
    echo "  GRN:         ${OUTDIR}/pycistopic/ + ${OUTDIR}/scenicplus/"
    echo "  Logs:        ${OUTDIR}/pipeline_info/"
    echo "========================================="

    exit ${exit_code}
fi
