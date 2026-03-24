// modules/rna/cellbender.nf
// CellBender ambient RNA removal
//
// HISTORY:
//   FIX-19: Resource allocation restored from legacy (128 GB / 16 cpus / 72h).
//           Previous 32 GB / 4 cpus / 6h caused repeated SLURM timeouts (exit 140).
//   FIX-22: CellBender 0.3.0 + PyTorch 2.4 broke torch.save() checkpoint
//           serialisation ("cannot pickle 'weakref.ReferenceType' object").
//           PyTorch 2.4 changed nn.Module hook dicts to use weakref-backed
//           containers.  Without checkpoint, posterior computation crashes,
//           so *_filtered.h5 and *_report.html are never generated.
//           Symptom: jobs exit 0, .h5 present, but missing report.html
//           triggers infinite Nextflow retries.
//     FIX-22a (attempted): cellbender_wrapper.py monkey-patched torch.save
//           to strip hook dicts.  Partial fix — weakrefs lived deeper in
//           Pyro's variational inference internals.
//     FIX-22b (attempted): copyreg handler to pickle weakref.ref as None,
//           plus hook stripping.  Worked but fragile across Pyro versions.
//     FIX-22c (CURRENT — container-level fix): Upgraded CellBender in the
//           scgpu_extended.sif container using the upstream fix from Broadinstitute:
//             pip install --no-cache-dir -U \
//               git+https://github.com/broadinstitute/CellBender.git@4334e89
//           This commit properly handles serialisation under PyTorch 2.x.
//           Confirmed working by multiple users with torch 2.5.1.
//           The wrapper script (bin/cellbender_wrapper.py) is no longer needed
//           but kept in the repo for reference.

process CELLBENDER {
    tag "${sample_id}"
    publishDir "${params.outdir}/cellbender", mode: 'copy'

    input:
    tuple val(sample_id), path(raw_h5)

    output:
    tuple val(sample_id), path("${sample_id}_cellbender_filtered.h5"), emit: filtered_h5
    tuple val(sample_id), path("${sample_id}_cellbender.h5"), emit: full_h5
    path("${sample_id}_cellbender_report.html"), emit: report

    script:
    """
    set -euo pipefail
    export LD_LIBRARY_PATH="\${CONDA_PREFIX:-/usr}/lib:\${LD_LIBRARY_PATH:-}"

    # ------------------------------------------------------------------
    # 1. Resolve input format (BD zip / 10x h5 / directory)
    # ------------------------------------------------------------------
    INPUT_PATH=""

    if [[ ${raw_h5} == *.zip ]]; then
        unzip -q ${raw_h5}

        if [ -f "matrix.mtx.gz" ] && [ -f "barcodes.tsv.gz" ] && [ -f "features.tsv.gz" ]; then
            # BD Rhapsody format — files at root, reorganise
            mkdir -p raw_feature_bc_matrix
            mv matrix.mtx.gz barcodes.tsv.gz features.tsv.gz raw_feature_bc_matrix/
            INPUT_PATH="raw_feature_bc_matrix"
        else
            # 10x format — look for standard directory
            INPUT_PATH=\$(find . -type d \\( -name "raw_feature_bc_matrix" -o -name "raw_gene_bc_matrices*" \\) | head -n 1)
        fi
    elif [[ ${raw_h5} == *.h5 ]] || [[ ${raw_h5} == *.h5ad ]]; then
        INPUT_PATH="${raw_h5}"
    elif [ -d "${raw_h5}" ]; then
        INPUT_PATH="${raw_h5}"
    fi

    if [ -z "\$INPUT_PATH" ]; then
        echo "ERROR: Could not resolve CellBender input from ${raw_h5}"
        ls -la
        exit 1
    fi

    # ------------------------------------------------------------------
    # 2. Cap total_droplets_included to barcode count (prevents CB crash)
    # ------------------------------------------------------------------
    TOTAL_DROPLETS=${params.cellbender.total_droplets}

    if [ -d "\$INPUT_PATH" ]; then
        N_BARCODES=0
        if [ -f "\$INPUT_PATH/barcodes.tsv.gz" ]; then
            N_BARCODES=\$(zcat "\$INPUT_PATH/barcodes.tsv.gz" | wc -l)
        elif [ -f "\$INPUT_PATH/barcodes.tsv" ]; then
            N_BARCODES=\$(wc -l < "\$INPUT_PATH/barcodes.tsv")
        fi
        if [ "\$N_BARCODES" -gt 0 ] && [ "\$N_BARCODES" -le "\$TOTAL_DROPLETS" ]; then
            TOTAL_DROPLETS=\$((N_BARCODES - 1))
            if [ "\$TOTAL_DROPLETS" -lt 1000 ]; then
                TOTAL_DROPLETS=1000
            fi
        fi
        echo "CellBender: sample ${sample_id} — directory input, \$N_BARCODES barcodes, total_droplets=\$TOTAL_DROPLETS"
    else
        echo "CellBender: sample ${sample_id} — .h5 file input, total_droplets=\$TOTAL_DROPLETS"
    fi

    # ------------------------------------------------------------------
    # 3. Detect GPU availability
    # ------------------------------------------------------------------
    CUDA_FLAG=""
    python -c "import torch; exit(0 if torch.cuda.is_available() else 1)" 2>/dev/null && CUDA_FLAG="--cuda"

    # ------------------------------------------------------------------
    # 4. Run CellBender
    #    FIX-22c: Native cellbender CLI — the upstream commit (4334e89)
    #    fixes PyTorch 2.x checkpoint serialisation.  No wrapper needed.
    #    If this regresses, see bin/cellbender_wrapper.py (FIX-22b) as
    #    a fallback that patches torch.save via copyreg.
    # ------------------------------------------------------------------
    cellbender remove-background \\
        --input \$INPUT_PATH \\
        --output ${sample_id}_cellbender.h5 \\
        --total-droplets-included \$TOTAL_DROPLETS \\
        --expected-cells ${params.cellbender.expected_cells} \\
        --fpr ${params.cellbender.fpr} \\
        --epochs ${params.cellbender.epochs} \\
        \$CUDA_FLAG \\
        --low-count-threshold ${params.cellbender.low_count_threshold}

    # CellBender should produce *_filtered.h5 automatically after
    # posterior computation.  If missing, fall back to copying the
    # main output (this was the pre-FIX-22c behaviour that masked
    # the checkpoint failure).
    if [ ! -f "${sample_id}_cellbender_filtered.h5" ]; then
        echo "WARNING: _filtered.h5 not created by CellBender — copying main output as fallback"
        cp ${sample_id}_cellbender.h5 ${sample_id}_cellbender_filtered.h5
    fi
    """
}
