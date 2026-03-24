#!/bin/bash
# atac_pipeline_wrapper.sh
# Wrapper to fix LD_LIBRARY_PATH before running ATAC pipeline

set -euo pipefail

# ========================================
# DYNAMIC PATH RESOLUTION
# ========================================
# Get the directory where THIS script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The Python script is in the same bin/ directory
PYTHON_SCRIPT="${SCRIPT_DIR}/atac_initial_qc.py"

# ========================================
# SET UP ENVIRONMENT
# ========================================
# Set up conda environment path (assumes pysnATAC_nf is always at this location)
# export LD_LIBRARY_PATH=/dfs7/swaruplab/lesolano/tools/pysnATAC_nf/lib:$LD_LIBRARY_PATH

# ========================================
# VERIFY SCRIPT EXISTS
# ========================================
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "ERROR: Python script not found at: $PYTHON_SCRIPT" >&2
    echo "Script directory: $SCRIPT_DIR" >&2
    exit 1
fi

# ========================================
# RUN THE PIPELINE
# ========================================
echo "Running: python3 $PYTHON_SCRIPT with $(($# )) arguments"
python3 "$PYTHON_SCRIPT" "$@"
