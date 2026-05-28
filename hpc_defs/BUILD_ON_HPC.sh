#!/bin/bash
# ============================================================================
# BUILD_ON_HPC.sh — build FORGE Singularity containers from .def files on HPC.
#
# Source recipes: ../docs/defs/<container>.def (the canonical, version-pinned
# Singularity definition files).
#
# Usage:
#   # Interactive (recommended for first build / debugging):
#   srun --partition=free --time=04:00:00 --mem=24G --cpus-per-task=8 --pty bash
#   module load singularity   # or: module load apptainer
#   cd /path/to/forge
#   bash hpc_defs/BUILD_ON_HPC.sh all
#
#   # Or a single container:
#   bash hpc_defs/BUILD_ON_HPC.sh seurat_extended
#
#   # Or via sbatch (uses the SBATCH directives below):
#   sbatch hpc_defs/BUILD_ON_HPC.sh all
#
# Options:
#   --rebuild   rebuild even if the .sif already exists (default: skip)
#   --no-test   skip the %test section in each .def (faster, less safe)
#
# If --fakeroot is not available on your HPC, the script falls back to an
# unprivileged build; some apt-get steps may need a remote build then.
# ============================================================================
#SBATCH --job-name=forge-containers
#SBATCH --partition=free
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --output=container_build_%j.log

set -uo pipefail

# ---- Argument parsing ----
TARGETS=()
REBUILD=0
NO_TEST=0
for arg in "$@"; do
    case "$arg" in
        --rebuild)  REBUILD=1 ;;
        --no-test)  NO_TEST=1 ;;
        -h|--help)
            sed -n '2,28p' "$0"
            exit 0 ;;
        *)          TARGETS+=("$arg") ;;
    esac
done
[ ${#TARGETS[@]} -eq 0 ] && TARGETS=("all")

# ---- Paths ----
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DEFS_DIR="$REPO_DIR/docs/defs"
OUTPUT_DIR="${FORGE_SIF_DIR:-$REPO_DIR/singularity_cache}"
LOG_DIR="$OUTPUT_DIR/build_logs"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# ---- Container list (build order = increasing complexity) ----
ALL_CONTAINERS=(
    cicero
    scgpu_extended
    snapatac_extended
    seurat_extended
    scenicplus
)

# Expand "all" → full list
EXPANDED=()
for t in "${TARGETS[@]}"; do
    if [ "$t" = "all" ]; then
        EXPANDED+=("${ALL_CONTAINERS[@]}")
    else
        EXPANDED+=("$t")
    fi
done

# ---- Module load (best-effort) ----
module load singularity 2>/dev/null || module load apptainer 2>/dev/null || true
BUILDER="$(command -v singularity || command -v apptainer || echo '')"
if [ -z "$BUILDER" ]; then
    echo "ERROR: neither singularity nor apptainer is on PATH after 'module load'."
    echo "       Try 'module avail singularity' and load the right module."
    exit 1
fi
echo "Builder: $BUILDER ($($BUILDER --version 2>&1 | head -1))"
echo "Output:  $OUTPUT_DIR"
echo "Logs:    $LOG_DIR"
echo "Targets: ${EXPANDED[*]}"
echo ""

# ---- Build one container ----
build_one() {
    local name="$1"
    local def="$DEFS_DIR/${name}.def"
    local sif="$OUTPUT_DIR/${name}.sif"
    local log="$LOG_DIR/${name}_build_$(date +%Y%m%d_%H%M%S).log"

    if [ ! -f "$def" ]; then
        echo "[$name] SKIP — definition file missing: $def"
        return 2
    fi

    if [ -f "$sif" ] && [ "$REBUILD" -eq 0 ]; then
        echo "[$name] SKIP — $sif already exists. Pass --rebuild to force."
        return 0
    fi
    [ -f "$sif" ] && rm -f "$sif"

    echo "[$name] Building from $def"
    echo "         logging to     $log"
    local t0
    t0=$(date +%s)

    # Try --fakeroot first, fall back to unprivileged.
    local build_args=()
    [ "$NO_TEST" -eq 1 ] && build_args+=("--notest")

    if $BUILDER build --fakeroot "${build_args[@]}" "$sif" "$def" >"$log" 2>&1; then
        echo "[$name] OK (--fakeroot)"
    elif $BUILDER build "${build_args[@]}" "$sif" "$def" >>"$log" 2>&1; then
        echo "[$name] OK (unprivileged)"
    else
        echo "[$name] FAIL — see $log"
        return 1
    fi

    local t1 elapsed size
    t1=$(date +%s)
    elapsed=$(( t1 - t0 ))
    size=$(du -h "$sif" | cut -f1)
    echo "[$name] $size in ${elapsed}s"
    return 0
}

# ---- Run + summarise ----
declare -A RESULT
declare -A SIZE
declare -A ELAPSED
START=$(date +%s)

for name in "${EXPANDED[@]}"; do
    t0=$(date +%s)
    if build_one "$name"; then
        RESULT[$name]="OK"
        if [ -f "$OUTPUT_DIR/${name}.sif" ]; then
            SIZE[$name]=$(du -h "$OUTPUT_DIR/${name}.sif" | cut -f1)
        else
            SIZE[$name]="-"
        fi
    else
        rc=$?
        RESULT[$name]=$([ "$rc" = "2" ] && echo "MISSING_DEF" || echo "FAIL")
        SIZE[$name]="-"
    fi
    ELAPSED[$name]=$(( $(date +%s) - t0 ))
done
END=$(date +%s)

# ---- Summary table ----
echo ""
echo "============================================================"
echo "  Build summary"
echo "  Total elapsed: $(( END - START ))s"
echo "============================================================"
printf "  %-22s  %-12s  %8s  %8s\n" "Container" "Result" "Size" "Time"
printf "  %-22s  %-12s  %8s  %8s\n" "----------------------" "------------" "--------" "--------"
EXIT_RC=0
for name in "${EXPANDED[@]}"; do
    printf "  %-22s  %-12s  %8s  %8s\n" \
        "$name" "${RESULT[$name]}" "${SIZE[$name]}" "${ELAPSED[$name]}s"
    [ "${RESULT[$name]}" != "OK" ] && EXIT_RC=1
done
echo ""

# ---- Hash manifest for diffing against expected ----
echo "============================================================"
echo "  SHA256 manifest (paste this back when comparing builds)"
echo "============================================================"
for name in "${EXPANDED[@]}"; do
    sif="$OUTPUT_DIR/${name}.sif"
    if [ -f "$sif" ]; then
        printf "  %-26s  %12s  %s\n" "${name}.sif" "$(stat -c %s "$sif")" "$(sha256sum "$sif" | awk '{print $1}')"
    fi
done

exit $EXIT_RC
