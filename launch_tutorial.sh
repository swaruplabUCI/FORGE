#!/bin/bash
#SBATCH --job-name=forge-tutorial
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=forge-tutorial-%j.log
#
# NOTE: --account and --partition are deliberately NOT set here. They are
# site-specific and there is no portable default. Pass them to sbatch:
#
#     sbatch -A <account> -p <partition> launch_tutorial.sh
#
# =========================================================================
# FORGE tier-2 tutorial launcher
#
# Runs the ~1,000-cell tutorial dataset end to end. Works two ways:
#
#   WITH SLURM      sbatch -A <account> -p <partition> launch_tutorial.sh
#   WITHOUT SLURM   ./launch_tutorial.sh
#
# Both paths run the SAME pipeline with the SAME `tutorial` profile. The only
# difference is that the SLURM path additionally loads
# configs/tutorial_slurm.config, which sizes Nextflow's local executor to the
# allocation (see the header of that file for why that is not optional).
#
# The whole pipeline runs inside ONE allocation using Nextflow's local
# executor. It does NOT submit a SLURM job per process — that would require
# site-specific account/partition/QOS settings in every withName block.
#
# Sizing: the heaviest single task measured 8.2 GB peak RSS and the tier
# requests at most 12 GB, so --mem=48G allows ~4 heavy tasks at once. Raise
# --cpus-per-task/--mem for more parallelism; the 45-way hdWGCNA and 25-way
# Cicero fan-outs are what benefit.
#
# Options:
#   --outdir DIR   publish results to DIR       (default: results_tutorial)
#   --tutorial_data DIR
#                  where the unpacked dataset lives (default: <repo>/tutorial_data)
#   --no-resume    ignore the cache, run cold   (default: -resume)
#   --preview      build the DAG and exit, ~15s (no compute, no containers run)
# =========================================================================

set -euo pipefail

OUTDIR=""
RESUME="-resume"
PREVIEW=""
TUTORIAL_DATA_ARG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --outdir)    OUTDIR="$2"; shift 2 ;;
        --tutorial_data) TUTORIAL_DATA_ARG="$2"; shift 2 ;;
        --no-resume) RESUME="";   shift ;;
        --preview)   PREVIEW="-preview"; RESUME=""; shift ;;
        # Print the header block up to its closing rule. Anchored on the rule
        # rather than a line number, which silently truncates the help the next
        # time a line is added above it.
        -h|--help)   awk 'NR>=14 { if (/^# ={10,}/) exit; print }' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

# --- Locate the repo -----------------------------------------------------
# $0 is unreliable under sbatch (SLURM copies the script), so prefer the
# submit directory when SLURM tells us what it was.
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$PROJECT_DIR"

if [[ ! -f main.nf || ! -f nextflow.config ]]; then
    echo "ERROR: main.nf / nextflow.config not found in $PROJECT_DIR" >&2
    echo "Run this from the FORGE repository root." >&2
    exit 1
fi

# --- Environment ---------------------------------------------------------
module load singularity 2>/dev/null || true

# Nextflow must be on PATH. The bundled launcher is preferred when present.
if [[ -x /dfs7/swaruplab/lesolano/tools/nextflow ]]; then
    export PATH="/dfs7/swaruplab/lesolano/tools:$PATH"
fi
if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: nextflow not found on PATH." >&2
    echo "Install it (https://nextflow.io) or add its directory to PATH." >&2
    exit 1
fi

export TMPDIR="${TMPDIR:-/tmp}"
export NXF_TEMP="$TMPDIR"
export HDF5_USE_FILE_LOCKING=FALSE
export SINGULARITY_BINDPATH="${SINGULARITY_BINDPATH:-$PROJECT_DIR,$TMPDIR}"

# --- Input check ---------------------------------------------------------
TUTORIAL_DATA="${TUTORIAL_DATA_ARG:-${PROJECT_DIR}/tutorial_data}"
if [[ ! -f "${TUTORIAL_DATA}/manifest.csv" ]]; then
    echo "ERROR: tutorial dataset not found at ${TUTORIAL_DATA}" >&2
    echo "Download and unpack the tutorial data release asset first;" >&2
    echo "see docs/tutorial.md. Override the location with --tutorial_data." >&2
    exit 1
fi

# --- Assemble the command ------------------------------------------------
CONFIGS=(-c configs/datasets/tutorial_pbmc.config)

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    # Inside an allocation: cap the local executor to what we were given.
    CONFIGS+=(-c configs/tutorial_slurm.config)
    MODE="SLURM job ${SLURM_JOB_ID} (${SLURM_CPUS_PER_TASK:-?} CPUs, ${SLURM_MEM_PER_NODE:-?} MB)"
else
    MODE="local (no SLURM; Nextflow will use this machine's full CPU/RAM)"
fi

OUTDIR_ARG=()
[[ -n "$OUTDIR" ]] && OUTDIR_ARG=(--outdir "$OUTDIR")

DATA_ARG=()
[[ -n "$TUTORIAL_DATA_ARG" ]] && DATA_ARG=(--tutorial_data "$TUTORIAL_DATA")

echo "========================================="
echo "FORGE tier-2 tutorial"
echo "  Project dir: ${PROJECT_DIR}"
echo "  Mode:        ${MODE}"
echo "  Outdir:      ${OUTDIR:-results_tutorial (profile default)}"
echo "  Dataset:     ${TUTORIAL_DATA}"
echo "  Resume:      ${RESUME:-disabled}"
[[ -n "$PREVIEW" ]] && echo "  *** PREVIEW ONLY — no compute ***"
echo "========================================="

set -x
nextflow run main.nf \
    -profile tutorial,singularity \
    "${CONFIGS[@]}" \
    ${PREVIEW} \
    ${RESUME} \
    "${OUTDIR_ARG[@]}" \
    "${DATA_ARG[@]}"
