#!/bin/bash
#SBATCH --job-name=tut-onramp-chromvar
#SBATCH -A vswarup_lab_gpu
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A30:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/dfs7/swaruplab/lesolano/src_FORGE/dev_notes/phase3/onramp_chromvar_%j.log
#
# =========================================================================
# Build the tutorial's ChromVAR on-ramp bundle.
#
# STATUS: NOT SHIPPED (decided 2026-08-07). Kept because it works and because
# it records two findings that cost real time to establish. Do not wire this
# into the tutorial without re-reading the second one.
#
# FINDING 1 — the bundle is useless under the tutorial's own settings.
# It builds correctly and is barcode-compatible (817/817 verified below), but
# EVERY consumer of ChromVAR output is gated behind something the tutorial
# disables:
#     VIS_CHROMVAR                  needs has_da_peaks (differential conditions)
#     EXTRACT_CHROMVAR_MOTIFS       needs params.scprinter.run
#     MAP_TF_TO_TARGET_GENES        needs params.scprinter.run
#     DIFFERENTIAL_TF_ACCESSIBILITY needs params.differential_tf.run + comparisons
#     ENHANCER_FOOTPRINTING_RECIPES needs the scPRINTER printer object
# The tutorial is a single sample with one condition_group and scprinter.run =
# false, so on-ramping this produces zero additional work. Making it useful
# would mean enabling scPRINTER — 54% of all compute across the four published
# datasets — or inventing conditions for a single-condition dataset. Neither is
# justified for a wiring demo.
#
# Resurrect this only if the tutorial ever enables scPRINTER or gains a second
# condition. Runtime is ~2 minutes on an A30.
#
# FINDING 2 — WHY IT IS NOT "TRIMMED FROM THE PRODUCTION PBMC RUN"
#
# PHASE3_HANDOFF.md Step 4 proposed shipping precomputed ChromVAR output taken
# from the real PBMC_Hs_10X_r5 run. That does not work. ChromVAR's outputs are
# a cells x motifs matrix keyed by the ATAC barcodes of the run that produced
# them, and every downstream consumer joins them back against that same run:
#
#   EXTRACT_CHROMVAR_MOTIFS(dev, ..., atac_cell_type_key)  -> per-cell-type,
#       so it needs the tutorial's cells and its cell_type column
#   MAP_TF_TO_TARGET_GENES(raw, motifs, cicero_ccan, gtf)  -> joins against the
#       tutorial's Cicero CCANs, which are chr21/chr22 peaks only
#   DIFFERENTIAL_TF_ACCESSIBILITY(dev, peak_matrix, tasks) -> takes the
#       tutorial's peak matrix alongside the deviations
#
# The production run has different barcodes (its own QC), a genome-wide peak
# universe, and CCANs on all chromosomes. Joining it to the tutorial's objects
# would overlap on almost nothing and would fail quietly — empty or near-empty
# results that still look like a successful run. That is the worst possible
# failure mode for a tutorial.
#
# So the bundle is generated FROM THE TUTORIAL DATASET ITSELF, by running the
# one stage the tutorial cannot run (ChromVAR is a hard cupy/rmm GPU dependency
# with no CPU fallback) once, here, on a GPU node. Users then on-ramp the
# result and everything joins correctly with no GPU.
#
# Usage:
#   sbatch dev_notes/phase3/build_onramp_chromvar.sh [PEAK_MATRIX]
#
# PEAK_MATRIX defaults to the tutorial run's own peak matrix. It MUST come from
# a tutorial run whose ATAC parameters match what ships, or the barcodes will
# not line up with what users produce. Verify before publishing:
#   compare the peak matrix used here against a fresh tutorial run's
#   results_tutorial/atac/final/peak_matrix.h5ad (same n_obs/n_vars, same
#   obs_names) — see the CHECK step at the end of this script.
# =========================================================================

set -euo pipefail

REPO=/dfs7/swaruplab/lesolano/src_FORGE
PEAK_MATRIX="${1:-${REPO}/results_tutorial/atac/final/peak_matrix.h5ad}"
OUTDIR="${REPO}/tutorial_onramp/chromvar"

CACHE_DIR=/dfs7/swaruplab/lesolano/ref/scprinter
PFMS="${CACHE_DIR}/JASPAR2022_core_nonredundant.jaspar"
GENOME=hg38
CHUNK=30000

if [[ ! -f "$PEAK_MATRIX" ]]; then
    echo "ERROR: peak matrix not found: $PEAK_MATRIX" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=========================================="
echo "Tutorial ChromVAR on-ramp"
echo "  peak matrix: $PEAK_MATRIX"
echo "  outdir:      $OUTDIR"
echo "  genome:      $GENOME"
echo "=========================================="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || {
    echo "ERROR: no GPU visible. ChromVAR has a hard cupy/rmm dependency." >&2
    exit 1
}

module load singularity 2>/dev/null || true

# Mirrors the container flags nextflow.config applies to GPU_CHROMVAR
# (label 'process_gpu' + containerOptions '--nv' + the snapatac container).
singularity exec --nv --contain --home /tmp \
    --bind /dfs7 --bind /tmp --bind /dev/shm \
    --env PYTHONNOUSERSITE=1 \
    --env HDF5_USE_FILE_LOCKING=FALSE \
    --env NUMBA_CACHE_DIR=/tmp/numba_cache \
    --env MPLCONFIGDIR=/tmp/matplotlib \
    --env XDG_CACHE_HOME=/tmp/cache \
    --env CUPY_CACHE_DIR=/tmp/cupy_cache \
    --env PATH="/opt/conda/bin:/usr/local/bin:/usr/bin:/bin:${REPO}/bin" \
    "${REPO}/singularity_cache/snapatac_extended.sif" \
    "${REPO}/bin/gpu_chromvar_nf.py" \
      --peak-matrix "$PEAK_MATRIX" \
      --cache-dir "$CACHE_DIR" \
      --pfms "$PFMS" \
      --genome "$GENOME" \
      --out-prefix "${OUTDIR}/chromvar_matrix" \
      --chunk-size "$CHUNK" \
      --da-peaks ''
# NOTE: --out-prefix MUST be absolute. `--contain --home /tmp` makes the
# container's CWD /tmp, so a relative prefix writes both h5ads inside the
# container and they vanish when it exits — the run still exits 0 and logs
# "Wrote ...", which is why this looked like a success the first time.

echo ""
echo "=== Produced ==="
ls -la "${OUTDIR}/chromvar_matrix_deviations.h5ad" "${OUTDIR}/chromvar_matrix.h5ad"

# --- CHECK: the bundle must line up with the tutorial's own ATAC cells -------
singularity exec --contain --home /tmp --bind /dfs7 --bind /tmp \
    --env PYTHONNOUSERSITE=1 --env HDF5_USE_FILE_LOCKING=FALSE \
    --env NUMBA_CACHE_DIR=/tmp/nb --env MPLCONFIGDIR=/tmp/mpl --env XDG_CACHE_HOME=/tmp/c \
    "${REPO}/singularity_cache/snapatac_extended.sif" python3 - <<PYEOF
import anndata as ad
pm  = ad.read_h5ad("$PEAK_MATRIX", backed="r")
dev = ad.read_h5ad("${OUTDIR}/chromvar_matrix_deviations.h5ad", backed='r')
raw = ad.read_h5ad("${OUTDIR}/chromvar_matrix.h5ad", backed='r')
print(f"peak matrix : {pm.n_obs} cells x {pm.n_vars} peaks")
print(f"deviations  : {dev.n_obs} cells x {dev.n_vars} motifs")
print(f"raw         : {raw.n_obs} cells x {raw.n_vars}")
shared = len(set(pm.obs_names) & set(dev.obs_names))
print(f"barcodes shared with peak matrix: {shared} / {pm.n_obs}")
assert shared == pm.n_obs, "BARCODE MISMATCH — this bundle is not usable as an on-ramp"
print("OK — on-ramp is barcode-compatible with the tutorial dataset")
PYEOF

echo ""
echo "Ship these two together (all-or-none pair):"
echo "  onramp.chromvar_deviations = chromvar_matrix_deviations.h5ad"
echo "  onramp.chromvar_raw        = chromvar_matrix.h5ad"
